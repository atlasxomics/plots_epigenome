import anndata
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns 
import scanpy as sc
import snapatac2 as snap
import squidpy as sq

from functools import lru_cache
from pathlib import Path
from plotly.graph_objs.layout import Title
from plotly.subplots import make_subplots
from typing import List

from lplots import submit_widget_state
from lplots.widgets.checkbox import w_checkbox
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.multiselect import w_multi_select
from lplots.widgets.row import w_row
from lplots.widgets.select import w_select
from lplots.widgets.text import w_text_input, w_text_output
from latch.ldata.path import LPath


w_text_output(content="# **ATX Spatial Epigenomics Report**")

w_text_output(content="""

<details>
<summary><i>details</i></summary>

To load in data, click the 'Select File' icon and choose a directorycontaining an AnnData objects from Latch Data.
The selected directory should contain at least one of the following files:
- adata_ge.h5ad: a file containing a SnapATAC2 AnnData Object with .X as a gene accessibility matrix,
- adata_motifs.h5ad: a file containing a SnapATAC2 AnnData Object with .X as a motif deviation matrix.
Loading data into memory may take a couple minutes for large datasets.

> For AtlasXomics data, the standard location in Data for compatible files is `latch:///snap_outs/project_name/`.

</details>

""")

# Globals ------------------------------------------------------------------

obsm_keys = ("X_umap", "spatial")
na_keys = ['barcode', 'on_off', 'row', 'col', 'xcor', 'ycor', 'score']
all_colors = (
    'Paired', 'Paired_r',
    'Set1', 'Set1_r',
    'Set2', 'Set2_r',
    'tab10', 'tab10_r',
    'tab20', 'tab20_r',
    'tab20b', 'tab20b_r',
    'tab20c', 'tab20c_r',
    'deep', 'muted', 'pastel', 'bright', 'dark', 'colorblind',
    "Alphabet", "Dark24", "Light24"
)


# Functions ----------------------------------------------------------------
def convert_feature_expression(anndata, feature_name):
    """Create a new column in .obs with feature value from .X.
    """
    try:
        feature_value = anndata[:, feature_name].X.flatten()
        anndata.obs[feature_name] = list(feature_value)
    except KeyError as e:
        print(f"Error {e}: No gene {feature_name} found in .var_names")


def create_proportion_dataframe(
    adata, group_by, stack_by, return_type="proportion"
):
    """Create a DataFrame for a proportion or raw count plot (stacked bar
    graph) from an AnnData object.
    """

    # Create a cross-tabulation (contingency table) to get the counts
    count_df = pd.crosstab(adata.obs[group_by], adata.obs[stack_by])

    if return_type == 'proportion':
        # Calculate proportions for each group
        result_df = count_df.div(count_df.sum(axis=1), axis=0)
    elif return_type == 'counts':
        # Return raw counts
        result_df = count_df
    else:
        raise ValueError("Invalid return_type. Use 'proportion' or 'counts'.")

    # Reshape the DataFrame to long format for easier plotting (stacked format)
    long_df = result_df.reset_index().melt(
      id_vars=group_by,
      value_name="value",
      var_name=stack_by
    )
    long_df.columns = ["group_by", "stack_by", "value"]

    return long_df


def create_violin_data(adata, group_by, plot_data, data_type="obs"):
    """Create a DataFrame from an AnnData object to be used for violin plots;
    returns pandas DataFrame with columns: 'group', 'value', and 'type'
    (either 'obs' or 'gene').
    """

    # Check if the data to plot is from .obs
    if data_type == "obs":

        # Extract the data from .obs
        values = adata.obs[plot_data]
        df = pd.DataFrame({
            'group': adata.obs[group_by],
            'value': values
        })

    # Check if the data to plot is gene expression from .X
    elif data_type == "gene":

        # Extract gene expression values for the gene
        values = adata[:, plot_data].X.flatten()

        df = pd.DataFrame({
            "group": adata.obs[group_by],
            "value": values
        })

    else:
        raise ValueError("data_type must be either 'obs' or 'gene'.")

    return df


def custom_plotly(
    snap_fig, color_scheme='bright', width=400, height=500, hide_axes=True
):

    orig_data = snap_fig.data
    unique_groups = sorted([d.name for d in orig_data if d.name is not None])
    n_groups = len(unique_groups)

    # Generate new colors
    colors = generate_color_palette(n_groups, color_scheme)

    group_color_map = {unique_groups[i]: colors[i] for i in range(n_groups)}

    # Create new figure
    new_fig = go.Figure()

    new_fig.update_layout(
        snap_fig.layout
    )
    # Add traces with new colors
    for i, trace in enumerate(orig_data):
        if trace.name is not None:  # Skip traces without names
            new_trace = go.Scatter(
                x=trace.x,
                y=trace.y,
                mode=trace.mode,
                name=trace.name,
                marker=dict(
                    size=trace.marker.size,
                    color=group_color_map[trace.name],
                    opacity=trace.marker.opacity if trace.marker.opacity else 0.7
                ),
                showlegend=trace.showlegend
            )
            new_fig.add_trace(new_trace)

    # Update only width and height
    new_fig.update_layout(
        width=width,
        height=height
    )

    if hide_axes:
        new_fig.update_layout(
            xaxis=dict(showticklabels=False, showline=False, zeroline=False, showgrid=False),
            yaxis=dict(showticklabels=False, showline=False, zeroline=False, showgrid=False)
        )

    return new_fig


def filter_adata_by_groups(adata, group, group_a, group_b="All"):
    """Filter adata to two values in obs."""
    assert group_a != group_b, "Groups must be different."
    return adata[adata.obs[group].isin([group_a, group_b])]


def generate_color_palette(length, scheme="bright"):
    """
    Generate a list of hex color codes  with the specified number of colors.
    
    length : int
        Number of colors in the palette
    scheme : str, default="bright"
        Name of the color scheme to use. Can be one of:
        - Matplotlib palettes: 'Paired', 'Set1', 'Set2', 'tab10', 'tab20', etc.
        - Seaborn palettes: 'deep', 'muted', 'pastel', 'bright', 'dark', 'colorblind'
        - Plotly palettes: 'Alphabet', 'Dark24', 'Light24'
        
    """

    
    # Helper function to convert RGB to hex
    def rgb_to_hex(rgb):
        return "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        )
    
    matplotlib_palettes = (
        "Paired", "Paired_r",
        "Set1", "Set1_r",
        "Set2", "Set2_r",
        "tab10", "tab10_r",
        "tab20", "tab20_r",
        "tab20b", "tab20b_r",
        "tab20c", "tab20c_r",
    )
    
    seaborn_palettes = (
        "deep", "muted", "pastel", "bright", "dark", "colorblind"
    )
    
    plotly_palettes = {
        "Alphabet": [
            "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
            "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD", 
            "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E", 
            "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", 
            "#FBE426", "#FA0087"
        ],
        "Dark24": [
            "#2E91E5", "#E15F99", "#1CA71C", "#FB0D0D", "#DA16FF", "#222A2A", 
            "#B68100", "#750D86", "#EB663B", "#511CFB", "#00A08B", "#FB00D1", 
            "#FC0080", "#B2828D", "#6C7C32", "#778AAE", "#862A16", "#A777F1", 
            "#620042", "#1616A7", "#DA60CA", "#6C4516", "#0D2A63", "#AF0038"
        ],
        "Light24": [
            "#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", 
            "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", "#D626FF", 
            "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", 
            "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72"
        ]
    }
    
    if scheme in matplotlib_palettes:
        if scheme == 'rainbow':
            cm = plt.cm.rainbow
        else:
            cm = plt.get_cmap(scheme)
        colors = [rgb_to_hex(cm(i)[:3]) for i in np.linspace(0, 1, length)]
    elif scheme in seaborn_palettes:
        colors = sns.color_palette(scheme, length)
        colors = [rgb_to_hex(c) for c in colors]
    elif scheme in plotly_palettes:
        # Get Plotly palette and cycle if needed
        plotly_colors = plotly_palettes[scheme]
        colors = [plotly_colors[i % len(plotly_colors)] for i in range(length)]
    else:
        cm = plt.get_cmap('viridis')
        colors = [rgb_to_hex(cm(i)[:3]) for i in np.linspace(0, 1, length)]
    
    return colors


def get_feature_heatmap(df, features, rank_by="logfoldchanges"):
    """Returns a heatmap DataFrame by selecting the selected features per group
    based on a ranking column.
    """
    groups = df["group"].unique()

    # Extract gene scores from each group
    feats_dfs = []
    for group in groups:

        # Filter the main DataFrame to only include the selected top genes
        feats_df = df[df["names"].isin(features)]

        feats_dfs.append(feats_df)

    # Combine all top n genes into a single DataFrame
    combined_df = pd.concat(feats_dfs, ignore_index=True)

    # Prepare the DataFrame for heatmap plotting
    pivot_df = combined_df[["group", "names", rank_by]]
    pivot_df["names"] = pd.Categorical(
      pivot_df["names"], categories=pivot_df["names"].unique(), ordered=True
    )

    heatmap_df = pivot_df.pivot_table(
      index='group', columns='names', values=rank_by, aggfunc='mean'
    )

    return heatmap_df


def get_groups(adata: anndata.AnnData) -> List[str]:
    """Set 'groups' list for differential analysis."""

    groups = []
    for group in ["cluster", "sample", "condition"]:
        length = 1
        try:
            length = len(adata.obs[group].unique())
        except KeyError as e:
            print(f"{e}")

        if length > 1:
            groups.append(group)

    return groups


def get_top_n_heatmap(df, rank_by="logfoldchanges", n_top=5):
    """Returns a heatmap DataFrame by selecting the top n features per group
    based on a ranking column.
    """

    groups = df["group"].unique()

    # Loop through each group and select the top n genes
    topn_gene_dfs = []
    for group in groups:

        group_df = df[df["group"] == group]

        if rank_by in ["logfoldchanges", "scores"]:
            topn = group_df.nlargest(n_top, columns=rank_by)["names"]
        elif rank_by in ["pvals", "pvals_adj"]:
            topn = group_df.nsmallest(n_top, columns=rank_by)["names"]

        # Filter the main DataFrame to only include the selected top genes
        gene_df = df[df["names"].isin(topn)]

        topn_gene_dfs.append(gene_df)

    # Combine all top n genes into a single DataFrame
    combined_df = pd.concat(topn_gene_dfs, ignore_index=True)

    # Prepare the DataFrame for heatmap plotting
    pivot_df = combined_df[["group", "names", rank_by]]
    pivot_df["names"] = pd.Categorical(
      pivot_df["names"], categories=pivot_df["names"].unique(), ordered=True
    )

    heatmap_df = pivot_df.pivot_table(
      index='group', columns='names', values=rank_by, aggfunc='mean'
    )

    return heatmap_df


def make_volcano_df(
    adata,
    group,
    group_a,
    group_b,
    feature,
    threshold=0.01,
    display_pval=True,
    display_gm=True
):
    """Using sc.get.rank_genes_groups_df, make dataframe for volcano plot.
    Replace p-values that are 0 or NaN with a small number.
    Optionally filter out genes with the "Gm" prefix.
    """

    assert group_a != group_b, "Groups must be different."
    assert group in adata.obs.columns, f"No group {group} for in AnnData."

    subgroups = adata.obs[group].unique()
    assert group_a in subgroups, f"Group A {group_a} not found in subgroups."

    key = f"{group_a}_{group_b}_genes"

    if group_b == "All":
        df = sc.get.rank_genes_groups_df(
            adata, group=group_a, key=f"{group}_{feature}"
        )
    else:
        assert group_b in subgroups, f"Group B {group_b} not found in subgroups."
        adata = filter_adata_by_groups(adata, group, group_a, group_b)
        sc.tl.rank_genes_groups(
            adata,
            groupby=group,
            method="t-test",
            key_added=key,
            use_raw=False
        )
        df = sc.get.rank_genes_groups_df(adata, group=group_a, key=key)

    # Filter out rows with NaN logfoldchanges
    df = df[~df["logfoldchanges"].isna()]

    if not display_gm:
        df = df[~df["names"].str.startswith("Gm")]

    if display_pval:
        # Ensure pvals_adj column exists
        if "pvals_adj" not in df.columns:
            raise ValueError("pvals_adj column is missing from the dataframe.")

        # Replace NaN with 0 (so we handle both 0s and NaNs consistently)
        df["pvals_adj"].fillna(0, inplace=True)

        # Find the smallest nonzero p-value
        min_nonzero_pval = df.loc[df['pvals_adj'] > 0, 'pvals_adj'].min()

        # If all p-values are zero or NaN, use the smallest possible float
        if min_nonzero_pval is None or np.isnan(min_nonzero_pval):
            min_replacement = np.finfo(float).eps
        elif min_nonzero_pval > threshold:
            min_replacement = np.finfo(float).eps
        else:
            min_replacement = min_nonzero_pval / 10  # Or use an even smaller fraction

        # Replace both 0 and NaN values with min_replacement
        df.loc[df['pvals_adj'] == 0, 'pvals_adj'] = min_replacement

    else:
        df = df[df['pvals_adj'] != 0]

    return df


def plot_umap_for_samples(
  adata,
  samples,
  color_by='cluster',
  pt_size=3,
  coords="spatial",
  flipY=True,
  color_scheme='bright'
):
    # Check if color_by is discrete or continuous
    obs_groups = sorted(adata.obs[color_by].unique())
    print(obs_groups)
    is_discrete = len(obs_groups) < 30

    # Generate color mapping
    colors = generate_color_palette(len(obs_groups), color_scheme)
    group_color_map = {
        obs_groups[i]: colors[i] for i in range(len(obs_groups))
    }

    # Calculate number of rows needed
    num_rows = (len(samples) - 1) // 3 + 1
    num_cols = min(len(samples), 3)

    # Create subplots with calculated rows and columns
    combined_fig = make_subplots(
      rows=num_rows,
      cols=num_cols,
      subplot_titles=samples,
      horizontal_spacing=0,
      vertical_spacing=0.05
    )

    # Track clusters that have already been added to the legend for discrete
    shown_clusters = set() if is_discrete else None

    # We have to flip the y axis to match scanpy plots
    if flipY:
        flipped_y = adata.obsm['spatial'].copy()
        flipped_y[:, 1] = -flipped_y[:, 1]
        adata.obsm['spatial_flippedY'] = flipped_y
        coords = 'spatial_flippedY'

    # Loop through each sample
    for i, sample in enumerate(samples):
        row = i // 3 + 1
        col = i % 3 + 1

        # Subset the data by sample
        sample_data = adata[adata.obs['sample'] == sample]

        # Create UMAP plot for the subset
        fig = snap.pl.umap(
            sample_data,
            color=color_by,
            use_rep=coords,
            marker_size=pt_size,
            show=False
        )

        # Add each trace to the combined figure
        for trace in fig.data:
            if is_discrete:
                trace['marker']['color'] = group_color_map[trace['name']]
                # Handle duplicates for discrete data
                if trace.name not in shown_clusters:
                    shown_clusters.add(trace.name)
                    combined_fig.add_trace(trace, row=row, col=col)
                else:
                    trace.showlegend = False  # Hide duplicate legend entries
                    combined_fig.add_trace(trace, row=row, col=col)
            else:
                # For continuous data, always show the legend
                combined_fig.add_trace(trace, row=row, col=col)

        # original_width = fig.layout.width
        # original_height = fig.layout.height

    subplot_width = 500
    subplot_height = 500

    combined_fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        autosize=False,
        width=subplot_width * num_cols,
        height=subplot_height * num_rows,
        margin=dict(
            l=0,   # left margin
            r=0,    # right margin
            t=25,   # top margin (just enough for titles)
            b=0,    # bottom margin
            pad=0   # padding between plots
        ),
        legend=dict(
          font=dict(size=18),
          itemsizing='constant',
          itemwidth=30,
          xanchor='right',
          yanchor='top',
          x=1.1
        )
      )

    combined_fig.update_coloraxes(
      colorscale='Spectral_r', colorbar_title=color_by
    )

    for i in range(1, num_rows + 1):
        for j in range(1, num_cols + 1):
            combined_fig.update_xaxes(
              # scaleanchor="y",
              scaleratio=1,
              showticklabels=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              row=i,
              col=j
            )
            combined_fig.update_yaxes(
              showticklabels=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              row=i,
              col=j
            )

    return combined_fig


def plot_volcano(
  vol_df,
  pvals_adj_threshold,
  log2fc_threshold,
  group_a,
  group_b,
  plot_width=1000,
  plot_height=425,
  top_n=2
):
    """Creates a volcano plot using Plotly with labels for the top n points by
    p-value, highest log2 fold change, and lowest log2 fold change, minimizing
    label overlap and alternating label positions.
    """
    fig_volcano_plot = go.Figure()

    # Add scatter plot
    fig_volcano_plot.add_trace(go.Scattergl(
        x=vol_df['logfoldchanges'],
        y=-np.log10(vol_df['pvals_adj']),
        mode='markers',
        marker=dict(
            size=5,
            color=np.where(
              (vol_df['pvals_adj'] < pvals_adj_threshold)
              & (abs(vol_df['logfoldchanges']) > log2fc_threshold),
              np.where(vol_df['logfoldchanges'] > 0, 'red', 'blue'),
              'grey'
            )
        ),
        text=vol_df['names'],
        hoverinfo='text',
        hovertemplate='<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(Adj P-value): %{y:.2f}<br>Adj P-value: %{customdata:.2e}<extra></extra>',
        customdata=vol_df['pvals_adj']
    ))

    # Add labels for top n points by p-value, highest and lowest log2 fold change
    top_pvals = vol_df.nsmallest(top_n, 'pvals_adj')
    top_logfc = vol_df.nlargest(top_n, 'logfoldchanges')
    top_logfc_neg = vol_df.nsmallest(top_n, 'logfoldchanges')
    lowest_logfc = vol_df.nlargest(top_n, 'logfoldchanges')
    top_points = pd.concat([top_pvals, top_logfc, top_logfc_neg, lowest_logfc]).drop_duplicates()

    annotations = []
    for i, row in enumerate(top_points.itertuples()):
        annotations.append(
            dict(
                x=row.logfoldchanges,
                y=-np.log10(row.pvals_adj),
                text=row.names,
                showarrow=True,
                arrowhead=2,
                ax=0,
                ay=(-40 if i % 2 == 0 else 40),  # Alternate label placement
                xanchor='auto',
                yanchor='auto',
                textangle=0,
                align='center'
            )
        )

    # Adjust annotations to minimize overlap
    fig_volcano_plot.update_layout(annotations=annotations)

    # Add horizontal line for p-value threshold
    fig_volcano_plot.add_hline(
        y=-np.log10(pvals_adj_threshold), line_dash="dash", line_color="grey"
    )

    # Add vertical lines for log2 fold change thresholds
    fig_volcano_plot.add_vline(x=-log2fc_threshold, line_dash="dash", line_color="grey")
    fig_volcano_plot.add_vline(x=log2fc_threshold, line_dash="dash", line_color="grey")

    # Update layout
    fig_volcano_plot.update_layout(
        title=f"Volcano Plot: {group_a} vs {group_b}",
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10(Adjusted P-value)",
        showlegend=False,
        width=plot_width,
        height=plot_height
    )

    return fig_volcano_plot


def rgb_to_hex(rgb):
    """Convert RGB tuple to hex color code"""
    return '#{:02x}{:02x}{:02x}'.format(
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255)
    )


# Select input data -----------------------------------------------------------

data_path = w_ldata_picker(
  label="atx_snap output folder",
  appearance={
    "placeholder": "Placeholder…"
  }
)

if data_path.value is None:
    adata_g = None
    adata_m = None
    adata = None
    exit()

if not data_path.value.is_dir():
    w_text_output(
        content="Selected resource must be a directory...",
        appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()

# Get .h5ad files -------------------------------------------------------------

adata_g = [f for f in data_path.value.iterdir() if "sm_ge.h5ad" in f.name()]
adata_m = [f for f in data_path.value.iterdir() if "sm_motifs.h5ad" in f.name()]

if len(adata_g) == 1:
    adata_g = adata_g[0]
elif len(adata_g) == 0:
    adata_g = None
    w_text_output(
        content="No file with suffix 'sm_ge.h5ad' (gene data) found in selected folder; please ensure the output folder contains a file ending in '_ge.h5ad'",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
elif len(adata_g) > 1:
    adata_g = None
    w_text_output(
        content="Multiple files with suffix 'sm_ge.h5ad' (gene data) found in selected folder; please ensure the output folder contains only one file ending in '_ge.h5ad'",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()

if len(adata_m) == 1:
    adata_m = adata_m[0]
elif len(adata_m) == 0:
    adata_m = None
    w_text_output(
        content="No file with suffix 'sm_motifs.h5ad' (motif data) found in selected folder; please ensure the output folder contains a file ending in '_motifs.h5ad'",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
elif len(adata_m) > 1:
    adata_m = None
    w_text_output(
        content="Multiple files with suffix 'sm_motifs.h5ad' (motif data) found in selected folder; please ensure the output folder contains only one file ending in '_motifs.h5ad'",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()

if adata_g is None and adata_m is None:
    exit()

# Download files --------------------------------------------------------------

w_text_output(
  content="Downloading files...",
  appearance={"message_box": "info"}
)
submit_widget_state()

for data in [adata_g, adata_m]:
    if data is not None:
        data.download(Path(data.name()), cache=True)

# Load files ------------------------------------------------------------------

w_text_output(
  content="Loading data into memory; this may take a few minutes...",
  appearance={"message_box": "info"}
)
submit_widget_state()

if adata_g is not None:
    adata_g = sc.read(Path(adata_g.name()))
    available_genes = list(adata_g.var_names)

    # Convert n_fragment to float for plotting
    adata_g.obs["n_fragment"] = adata_g.obs["n_fragment"].astype(float)

    w_text_output(
        content=f"Successfully loaded data with {adata_g.n_obs} cells and {adata_g.n_vars} genomic features.",
        appearance={"message_box": "success"}
    )
    submit_widget_state()

if adata_m is not None:
    adata_m = sc.read(Path(adata_m.name()))
    available_motifs = list(adata_m.var_names)

    w_text_output(
      content=f"Successfully loaded data with {adata_m.n_obs} cells and {adata_m.n_vars} motifs.",
      appearance={"message_box": "success"}
    )
    submit_widget_state()

# Set default values ----------------------------------------------------------

if adata_g is not None:
    adata = adata_g
elif adata_m is not None:
    adata = adata_m
else:
    adata = None
    exit()

# Downsample adata_g to only highly variable genes for volcano plot.
if "highly_variable" not in adata_g.var.keys():
    w_text_output(
      content="Creating downsampled data for volcano plots...",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    sc.pp.highly_variable_genes(adata_g, n_top_genes=2000)
adata_hvg = adata_g[:, adata_g.var["highly_variable"]]

samples = adata.obs["sample"].unique()
groups = get_groups(adata)
available_metadata = tuple(key for key in adata.obs_keys()
                           if key not in na_keys)

filtered_groups: dict[str, dict[str, anndata.AnnData]] = {}

gvol_cache: dict[str, pd.DataFrame] = {}
mvol_cache: dict[str, pd.DataFrame] = {}

group_options = dict()
for group in groups:
    group_options[group] = list(adata_g.obs[group].unique())
