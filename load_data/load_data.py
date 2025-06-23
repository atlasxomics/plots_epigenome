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

from anndata import AnnData
from functools import lru_cache
from pathlib import Path
from plotly.graph_objs.layout import Title
from plotly.subplots import make_subplots
import scipy.cluster.hierarchy as sch
from typing import Any, Dict, List, Optional

from lplots import submit_widget_state
from lplots.widgets.button import w_button
from lplots.widgets.checkbox import w_checkbox
from lplots.widgets.column import w_column
from lplots.widgets.igv import w_igv, IGVOptions
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
def adjust_pvals(
    df: pd.DataFrame,
    pval_col: str = "pvals_adj",
    threshold: float = 0.05,
    display_pval: bool = True
) -> pd.DataFrame:
    """
    Prepare adjusted p‐values for plotting or filtering.

    If `display_pval=True`, zeros and NaNs in `pval_col` are replaced by
    a small positive number so that –log10(p) or color scales don’t blow up.
    If `display_pval=False`, rows where `pval_col` is zero or NaN are dropped.

    Parameters
    ----------
    df
        Input DataFrame (will be copied).
    pval_col
        Name of the adjusted‐p‐value column.
    threshold
        Any nonzero pval above this is treated as too big; we’ll use eps instead.
    display_pval
        If True, replace 0/NaN with eps or fraction of the minimum nonzero pval;
        if False, drop 0/NaN rows entirely.

    Returns
    -------
    A new DataFrame with `pval_col` massaged as above.
    """
    df = df.copy()
    if pval_col not in df.columns:
        raise ValueError(f"Column '{pval_col}' not found in DataFrame.")

    if display_pval:
        # fill NaNs with 0
        df[pval_col] = df[pval_col].fillna(0)

        # find smallest nonzero
        nonzero = df.loc[df[pval_col] > 0, pval_col]
        min_nonzero = nonzero.min() if not nonzero.empty else None

        # choose replacement
        if min_nonzero is None or np.isnan(min_nonzero) or min_nonzero > threshold:
            replacement = np.finfo(float).eps
        else:
            replacement = min_nonzero / 10

        # replace zeros (and any negatives, just in case)
        df.loc[df[pval_col] <= 0, pval_col] = replacement

    else:
        # drop any rows where pval is zero or NaN
        df = df[df[pval_col].notna() & (df[pval_col] > 0)].copy()

    return df


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


def filter_anndata(
    adata: AnnData, group: str, subgroup: List[str], mem=False
) -> AnnData:
    if mem:
        return adata[adata.obs[group] == subgroup].to_memory()
    return adata[adata.obs[group] == subgroup]


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


def get_feature_heatmap(df, features, rank_by="scores"):
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


def get_top_n_heatmap(df, rank_by="scores", n_top=5):
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
    display_gm=True
):
    """Using sc.get.rank_genes_groups_df, make dataframe for volcano plot.
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

    if not display_gm:
        df = df[~df["names"].str.startswith("Gm")]

    return df


def plot_neighborhood_groups(
  group_adatas: Dict["str", anndata.AnnData],
  title: str,
  key: str = "cluster",
  method: str = "None",
  uns_key: Optional[str] = None,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
):
  groups = list(group_adatas.keys())
  num_groups = len(groups)

  # Dynamically determine optimal grid layout based on number of groups
  if num_groups <= 2:
    num_cols = num_groups
    num_rows = 1
    base_width = 550 * num_cols
    base_height = 506 * num_rows
  elif num_groups <= 4:
    num_cols = 2
    num_rows = (num_groups + 1) // 2
    base_width = 550 * num_cols
    base_height = 506 * num_rows
  elif num_groups <= 9:
    num_cols = 3
    num_rows = (num_groups + 2) // 3
    base_width = 400 * num_cols
    base_height = 360 * num_rows
  else:
    num_cols = 4
    num_rows = (num_groups + 3) // 4
    base_width = 336 * num_cols
    base_height = 302 * num_rows

  # Create subplots with calculated rows and columns
  combined_fig = make_subplots(
    rows=num_rows,
    cols=num_cols,
    subplot_titles=groups,
    horizontal_spacing=0.05,  # Reduced horizontal spacing
    vertical_spacing=0.05     # Reduced vertical spacing
  )

  # Create a shared colorscale range for all subplots if not provided
  if vmin is None or vmax is None:
    all_mins = []
    all_maxs = []
    for group in groups:
      if uns_key in group_adatas[group].uns:
        data = group_adatas[group].uns[uns_key][mode]
        all_mins.append(np.nanmin(data))
        all_maxs.append(np.nanmax(data))

    if vmin is None and all_mins:
      vmin = min(all_mins)
    if vmax is None and all_maxs:
      vmax = max(all_maxs)

  # Create dynamic colorscale with white at zero
  abs_max = max(abs(vmin) if vmin is not None else 0, abs(vmax) if vmax is not None else 0)

  # If all values are positive or all negative, create appropriate one-sided colorscale
  if vmin is not None and vmax is not None:
    if vmin >= 0:  # All positive values
      custom_colorscale = [
        [0, 'white'],
        [1, 'red']
      ]
    elif vmax <= 0:  # All negative values
      custom_colorscale = [
        [0, 'blue'],
        [1, 'white']
      ]
    else:  # Mixed positive and negative values
      # Calculate the midpoint (0) in the normalized scale
      midpoint = abs(vmin) / (abs(vmin) + abs(vmax))

      # Create a colorscale with white at the midpoint
      custom_colorscale = [
        [0, 'blue'],
        [midpoint, 'white'],
        [1, 'red']
      ]
  else:
    # Default to RdBu_r if bounds not determined
    custom_colorscale = "RdBu_r"

  # Loop through each sample
  for i, group in enumerate(groups):
    row = (i // num_cols) + 1
    col = (i % num_cols) + 1

    # Get data directly rather than creating a full figure
    adata = group_adatas[group]

    # Validate input
    if uns_key not in adata.uns:
      continue

    # Get the data
    data = adata.uns[uns_key][mode]
    categories = adata.obs[key].cat.categories

    # Apply clustering or sorting if requested
    if method != "None":
      # Create a copy of data for clustering
      cluster_data = np.array(data, dtype=float)
      cluster_data = np.nan_to_num(cluster_data, nan=0.0, posinf=0.0, neginf=0.0)

      try:
        # Compute linkage matrices
        row_linkage = sch.linkage(data, method=method)
        col_linkage = sch.linkage(data.T, method=method)

        # Get dendrograms
        row_dendrogram = sch.dendrogram(row_linkage, no_plot=True)
        col_dendrogram = sch.dendrogram(col_linkage, no_plot=True)

        # Reorder data according to clustering
        row_order = row_dendrogram['leaves']
        col_order = col_dendrogram['leaves']
        data = data[row_order][:, col_order]
        categories = categories[row_order]
      except Exception as e:
        print(f"Warning: Clustering failed for {group} ({str(e)}). Proceeding without clustering.")
        method = "None"

    # If method is None, sort data numerically by category
    if method == "None":
      # Sort categories and data numerically
      sorted_indices = np.argsort([float(cat) if cat.replace('.', '', 1).isdigit() else float('inf') for cat in categories])
      data = data[sorted_indices]
      categories = categories[sorted_indices]

    # Create heatmap - only the first subplot will show a colorbar
    heatmap = go.Heatmap(
      z=data,
      x=categories,
      y=categories,
      colorscale=custom_colorscale,
      showscale=(i == 0),  # Only show colorbar for the first subplot
      textfont={"size": 12},  # Smaller font for dense plots
      hoverongaps=False,
      zmin=vmin,
      zmax=vmax,
      colorbar=dict(
        title=mode,
        tickformat='.2f',
        len=0.9,
        x=1.02,  # Position colorbar slightly to the right
        yanchor="middle"
      ) if i == 0 else None
    )

    # Add trace to combined figure
    combined_fig.add_trace(heatmap, row=row, col=col)

    # Update axes for each subplot - make more compact
    combined_fig.update_xaxes(
      showgrid=False,
      title=None,
      side='bottom',
      tickfont=dict(size=10),  # Smaller font size for tick labels
      row=row,
      col=col
    )
    combined_fig.update_yaxes(
      showgrid=False,
      title=None,
      autorange='reversed',
      tickfont=dict(size=10),  # Smaller font size for tick labels
      row=row,
      col=col
    )

  # For very large numbers of samples, increase the base size
  if num_groups > 16:
    base_width = max(base_width, 1000)
    base_height = max(base_height, 900)

  # Update layout once for all subplots
  combined_fig.update_layout(
    title={
      'text': title,
      'x': 0.5,  # Center the title
      'xanchor': 'center',
      'yanchor': 'top',
      'font': {'size': 18}  # Slightly larger font size for title
    },
    plot_bgcolor='rgba(0,0,0,0)',
    autosize=False,  # Explicitly set size
    width=base_width,
    height=base_height,
    margin=dict(l=80, r=80, t=100, b=40),  # Tighter margins
    showlegend=False
  )

  # Update subplot titles with smaller font
  for i in range(len(combined_fig.layout.annotations)):
    combined_fig.layout.annotations[i].font.size = 14

  return combined_fig


def plot_umap_for_samples(
  adata,
  samples,
  color_by='cluster',
  pt_size=3,
  coords="spatial",
  flipY=True,
  color_scheme='bright',
  show_cluster=None,
  vmin=None,
  vmax=None
):
    import numpy as np
    import pandas as pd
    from plotly.subplots import make_subplots

    # Determine if color_by is discrete or continuous
    obs_values = adata.obs[color_by]
    is_discrete = pd.api.types.is_categorical_dtype(obs_values) or obs_values.dtype.name == 'category' or obs_values.nunique() < 30

    print(f"Coloring by: {color_by} (Discrete: {is_discrete})")

    if is_discrete:
        obs_groups = sorted(obs_values.unique())
        colors = generate_color_palette(len(obs_groups), color_scheme)
        group_color_map = {obs_groups[i]: colors[i] for i in range(len(obs_groups))}
    else:
        group_color_map = None

    num_rows = (len(samples) - 1) // 3 + 1
    num_cols = min(len(samples), 3)

    combined_fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        subplot_titles=samples,
        horizontal_spacing=0,
        vertical_spacing=0.05
    )

    shown_clusters = set() if is_discrete else None

    if flipY:
        flipped_y = adata.obsm['spatial'].copy()
        flipped_y[:, 1] = -flipped_y[:, 1]
        adata.obsm['spatial_flippedY'] = flipped_y
        coords = 'spatial_flippedY'

    for i, sample in enumerate(samples):
        row = i // 3 + 1
        col = i % 3 + 1

        sample_data = adata[adata.obs['sample'] == sample]

        def apply_trace_settings(trace):
            if is_discrete:
                trace['marker']['color'] = group_color_map.get(trace.name, 'black')
                if trace.name not in shown_clusters:
                    shown_clusters.add(trace.name)
                else:
                    trace.showlegend = False
            else:
                if vmin is not None:
                    trace['marker']['cmin'] = vmin
                if vmax is not None:
                    trace['marker']['cmax'] = vmax

        if show_cluster is not None:
            highlight_mask = sample_data.obs['cluster'] == show_cluster
            highlight_data = sample_data[highlight_mask]
            background_data = sample_data[~highlight_mask]

            bg_fig = snap.pl.umap(
                background_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in bg_fig.data:
                trace['marker']['color'] = 'lightgrey'
                trace['marker']['opacity'] = 0.3
                trace.showlegend = False
                combined_fig.add_trace(trace, row=row, col=col)

            fg_fig = snap.pl.umap(
                highlight_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fg_fig.data:
                apply_trace_settings(trace)
                combined_fig.add_trace(trace, row=row, col=col)

        else:
            fig = snap.pl.umap(
                sample_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fig.data:
                apply_trace_settings(trace)
                combined_fig.add_trace(trace, row=row, col=col)

    subplot_width = 500
    subplot_height = 500
    combined_fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        autosize=False,
        width=subplot_width * num_cols,
        height=subplot_height * num_rows,
        margin=dict(l=0, r=0, t=25, b=0, pad=0),
        legend=dict(
            font=dict(size=18),
            itemsizing='constant',
            itemwidth=30,
            xanchor='right',
            yanchor='top',
            x=1.1
        )
    )

    if not is_discrete:
        combined_fig.update_coloraxes(
            colorscale='Spectral_r',
            colorbar_title=color_by,
            cmin=vmin,
            cmax=vmax
        )

    for i in range(1, num_rows + 1):
        for j in range(1, num_cols + 1):
            combined_fig.update_xaxes(
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
  pval_key="pvals_adj",
  l2fc_key="logfoldchanges",
  names_key="name",
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
        x=vol_df[l2fc_key],
        y=-np.log10(vol_df[pval_key]),
        mode='markers',
        marker=dict(
            size=5,
            color=np.where(
              (vol_df[pval_key] < pvals_adj_threshold)
              & (abs(vol_df[l2fc_key]) > log2fc_threshold),
              np.where(vol_df[l2fc_key] > 0, 'red', 'blue'),
              'grey'
            )
        ),
          text=vol_df[names_key],
        hoverinfo='text',
        hovertemplate='<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(Adj P-value): %{y:.2f}<br>Adj P-value: %{customdata:.2e}<extra></extra>',
        customdata=vol_df[pval_key]
    ))

    # Add labels for top n points by p-value, highest and lowest log2 fold change
    top_pvals = vol_df.nsmallest(top_n, pval_key)
    top_logfc = vol_df.nlargest(top_n, l2fc_key)
    top_logfc_neg = vol_df.nsmallest(top_n, l2fc_key)
    lowest_logfc = vol_df.nlargest(top_n, l2fc_key)
    top_points = pd.concat([top_pvals, top_logfc, top_logfc_neg, lowest_logfc]).drop_duplicates()

    annotations = []
    for i, row in top_points.iterrows():
        annotations.append(
          dict(
              x=row[l2fc_key],
              y=-np.log10(row[pval_key]),
              text=row[names_key],
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
        xaxis_title=l2fc_key,
        yaxis_title=f"-Log10({pval_key})",
        showlegend=False,
        width=plot_width,
        height=plot_height
    )

    return fig_volcano_plot


def plot_ranked_feature_plotly(
    df: pd.DataFrame,
    y_col: str,
    x_col: Optional[str] = None,
    label_col: Optional[str] = None,
    color_col: Optional[str] = None,
    n_labels: int = 30,
    colorscale: str = "Viridis",
    marker_size: int = 10,
    text_font_size: int = 10,
    ascending: bool = False,
    title: Optional[str] = None,
    y_label: Optional[str] = None,
) -> go.Figure:
    """
    Interactive Plotly scatter of y_col vs. x_col (or auto-ranked), colored by color_col or y_col,
    with top and bottom n_labels points labeled using annotations, no legend, optional title, 
    and custom axis labels.
    """
    if label_col is None:
        raise ValueError("`label_col` must be provided to draw text labels.")
    if color_col is not None and color_col not in df.columns:
        raise ValueError(f"`color_col` {color_col!r} not in DataFrame columns.")

    df_plot = df.copy()
    # Determine x-axis or compute rank
    if x_col is None:
        df_plot["rank"] = df_plot[y_col].rank(method="first", ascending=ascending)
        x = "rank"
    else:
        if x_col not in df_plot.columns:
            raise ValueError(f"`x_col` {x_col!r} not in DataFrame columns.")
        x = x_col

    # Decide which column to color by
    ccol = color_col if color_col is not None else y_col
    cbar_title = ccol.replace("_", " ").title()
    
    df_plot[x] = df_plot[x].astype(float)
    df_plot[y_col] = df_plot[y_col].astype(float)    
    df_plot[ccol] = df_plot[ccol].astype(float)

    # Create the single scatter trace for all points
    fig = go.Figure(
        go.Scatter(
            x=df_plot[x],
            y=df_plot[y_col],
            mode='markers',
            showlegend=False,
            marker=dict(
                size=marker_size,
                color=df_plot[ccol],
                colorscale=colorscale,
                showscale=True,
                opacity=0.8,
                colorbar=dict(title=cbar_title),
            ),
            hovertemplate=(
                f"{x}: %{{x}}<br>"
                f"{y_col}: %{{y}}<br>"
                f"{label_col}: %{{customdata}}<br>"
                f"{cbar_title}: %{{marker.color}}"
                "<extra></extra>"
            ),
            customdata=df_plot[label_col].astype(str).values,
        )
    )
    
    # Sort by y value to get top/bottom points
    df_sorted = df_plot.sort_values(y_col, ascending=not ascending)
    
    # Select top and bottom points for labeling
    top_points = df_sorted.head(n_labels)
    bottom_points = df_sorted.tail(n_labels) if n_labels > 0 else pd.DataFrame()

    # Function to add annotations with alternating left-right positions
    def add_point_annotations(points_df, is_top=True):
        if points_df.empty:
            return

        # Define annotation positions that alternate left and right
        # Using different offsets for top vs bottom points for better spacing
        positions = [
            {'ax': 40, 'ay': 10},    # right
            {'ax': -40, 'ay': 10},   # left
            {'ax': 40, 'ay': -10},   # left
            {'ax': -40, 'ay': -10},    # right
        ]

        for i, (_, row) in enumerate(points_df.iterrows()):
            position = positions[i % len(positions)]

            fig.add_annotation(
                x=row[x],
                y=row[y_col],
                text=str(row[label_col]),
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1,
                arrowcolor="#636363",
                ax=position['ax'],
                ay=position['ay'],
                font=dict(size=text_font_size),
                bgcolor="rgba(255, 255, 255, 0.7)",
                bordercolor="#c7c7c7",
                borderwidth=1,
                borderpad=2,
                standoff=2,
            )

    # Add annotations for top and bottom points
    add_point_annotations(top_points, is_top=True)
    add_point_annotations(bottom_points, is_top=False)

    # Calculate padding for x-axis
    x_min, x_max = df_plot[x].min(), df_plot[x].max()
    pad = (x_max - x_min) * 0.05

    # Determine axis titles
    xaxis_title = "Rank" if x_col is None else x.replace("_", " ").title()
    yaxis_title = y_label if y_label is not None else y_col.replace("_", " ").title()

    # Final layout
    layout_kwargs = dict(
        showlegend=False,
        xaxis=dict(
            title=xaxis_title,
            range=[x_min - pad, x_max + pad],
            automargin=True,
        ),
        yaxis=dict(
            title=yaxis_title,
            automargin=True,
        ),
        template="plotly_white",
        margin=dict(l=80, r=80, t=60, b=80),  # Increased margins for labels
    )
    if title:
        layout_kwargs["title"] = dict(text=title, x=0.5)

    fig.update_layout(**layout_kwargs)
    return fig


def plotly_heatmap(
  adata: AnnData,
  key: str = "cluster",
  title: str = "",
  method: str = "None",
  colorscale: str = "RdBu_r",
  width: Optional[int] = 700,
  height: Optional[int] = 700,
  uns_key: Optional[str] = None,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
  **kwargs: Any,
) -> go.Figure:

  # Validate input
  if key not in adata.obs:
      raise ValueError(f"Key '{key}' not found in adata.obs")
  if uns_key not in adata.uns:
      raise ValueError(f"Key '{uns_key}' not found in adata.uns")
   # Ensure the key column is categorical
  if not pd.api.types.is_categorical_dtype(adata.obs[key]):
      # Convert to categorical if it's not already
      adata.obs[key] = adata.obs[key].astype('category')

  # Retrieve data from .uns if specified
  if uns_key:
      # Retrieve the data from .uns
      array = adata.uns[uns_key][mode]

      # Create a new AnnData object with the retrieved array
      ad = AnnData(
          X=array,
          obs={
              key: pd.Categorical(adata.obs[key].cat.categories)
          }
      )
  else:
      # Use original adata if no uns_key is provided
      ad = adata

  # Process data
  data = ad.X
  categories = ad.obs[key].cat.categories

  # Set color scale range
  if vmin is None:
      vmin = np.nanmin(data)
  if vmax is None:
      vmax = np.nanmax(data)

  if method != "None":

    # Create a copy of data for clustering
    cluster_data = np.array(data, dtype=float)
    cluster_data = np.nan_to_num(cluster_data, nan=0.0, posinf=0.0, neginf=0.0)

    try:
      # Compute linkage matrices
      row_linkage = sch.linkage(data, method=method)
      col_linkage = sch.linkage(data.T, method=method)

      # Get dendrograms
      row_dendrogram = sch.dendrogram(row_linkage, no_plot=True)
      col_dendrogram = sch.dendrogram(col_linkage, no_plot=True)

      # Reorder data according to clustering
      row_order = row_dendrogram['leaves']
      col_order = col_dendrogram['leaves']
      data = data[row_order][:, col_order]
      categories = categories[row_order]
    except Exception as e:
      print(f"Warning: Clustering failed ({str(e)}). Proceeding without clustering.")
      method = "None"

  # If method is None, sort data numerically by category
  if method == "None":
      # Sort categories and data numerically
      sorted_indices = np.argsort([float(cat) if cat.replace('.', '', 1).isdigit() else float('inf') for cat in categories])
      data = data[sorted_indices]
      categories = categories[sorted_indices]

  fig = go.Figure()

  # Create heatmap
  heatmap = go.Heatmap(
      z=data,
      x=categories,
      y=categories,
      colorscale=colorscale,
      showscale=True,
      textfont={"size": 12},
      hoverongaps=False,
      zmin=vmin,  # Set minimum color scale value
      zmax=vmax,  # Set maximum color scale value
      colorbar=dict(
          title=key,
          tickformat='.2f',
          len=0.9
      ),
      **kwargs
  )
  fig.add_trace(heatmap)
  # Update layout
  fig.update_layout(
      title={
          'text': title,
          'x': 0.5,
          'xanchor': 'center',
          'yanchor': 'top'
      },
      width=width,
      height=height,
      showlegend=False,
      xaxis=dict(
          showgrid=False,
          side='bottom'
      ),
      yaxis=dict(
          showgrid=False,
          autorange='reversed'
      ),
  )

  return fig


def squidpy_analysis(
  adata: anndata.AnnData, cluster_key: str = "cluster"
) -> anndata.AnnData:
  """Perform squidpy Neighbors enrichment analysis.
  """

  if not adata.obs["cluster"].dtype.name == "category":
      adata.obs["cluster"] = adata.obs["cluster"].astype("category")

  sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=4, n_rings=1)
  sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

  return adata


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

samples = adata.obs["sample"].unique()
groups = get_groups(adata)

for data in [adata_g, adata_m]:
    for group in groups:
        if adata_g.obs[group].dtype != object:  # Ensure groups are str
            adata_g.obs[group] = adata_g.obs[group].astype(str)

available_metadata = tuple(key for key in adata.obs_keys()
                           if key not in na_keys)

filtered_groups: dict[str, dict[str, anndata.AnnData]] = {}

gvol_cache: dict[str, pd.DataFrame] = {}
mvol_cache: dict[str, pd.DataFrame] = {}

group_options = dict()
for group in groups:
    group_options[group] = list(adata_g.obs[group].unique())

clusters = group_options["cluster"]

# Stuff for IGV  ------------------------------------------------------------

coverages_dict = {}
for group in groups:
    for file in data_path.value.iterdir():
        if file.path.endswith(f"{group}_coverages"):
            coverages_dict[group] = file

if len(coverages_dict) > 0:
    w_text_output(
      content=f"Found coverage folders for {' '.join(list(coverages_dict.keys()))}",
      appearance={"message_box": "success"}
    )
    submit_widget_state()
else:
    w_text_output(
        content="No coverage folders were found for project...",
        appearance={"message_box": "warning"}
    )
