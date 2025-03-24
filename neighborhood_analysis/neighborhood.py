from typing import Any, Optional, Dict
import numpy as np
import plotly.graph_objects as go
from anndata import AnnData
import scipy.cluster.hierarchy as sch
from plotly.subplots import make_subplots
import pandas as pd


# --------------------------------------------------------------------------------
def filter_anndata(
    adata: AnnData, group: str, subgroup: List[str], mem=False
) -> AnnData:
    if mem:
        return adata[adata.obs[group] == subgroup].to_memory()
    return adata[adata.obs[group] == subgroup]


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

# --------------------------------------------------------------------------------

if not adata_g:
    w_text_output(
        content="No data gene activity data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_g, anndata.AnnData):
    w_text_output(
       content="No gene activity data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

neighbor_groups = [g for g in groups if g != "cluster"]
group_dict = {g: adata.obs[g].unique() for g in neighbor_groups}

neigh_group_by = w_select(
  label="subplot groups",
  default="all",
  options=tuple(["all"] + list(group_dict.keys())),
  appearance={
    "detail": "(all, sample, condition)",
    "help_text": "Subgroup plots from neighborhood data"
  }
)

mode = w_select(
  label="displayed data",
  default="zscore",
  options=("zscore", "count"),
  appearance={
    "help_text": "Data to be plotted"
  }
)

clustering_method = w_select(
  label="hierarchical clustering method",
  default="None",
  options=(
    "None",
    "single",
    "complete",
    "average",
    "weighted",
    "centroid",
    "median",
    "ward",
  ),
  appearance={
    "help_text": "see scipy.cluster.hierarchy.linkage"
  }
)

scale_max = w_text_input(
  label="colorscale maximum",
  default=None,
  appearance={
    "help_text": "Maximum value of colorscale"
  }
)

scale_min = w_text_input(
  label="colorscale minimum",
  default=None,
  appearance={
  "help_text": "Maximum value of colorscale"
  }
)

w_row(items=[neigh_group_by, mode, clustering_method, scale_max, scale_min])

vmax = int(scale_max.value) if scale_max.value and scale_max.value.strip().isdigit() else None
try:  # Handle negative values
  vmin = int(scale_min.value) if scale_min.value.strip() != '' else None
except ValueError:
  vmin = None  # Fallback if the value can't be converted to an integer

# --------------------------------------------------------------------------------

if neigh_group_by.value == "all":
  try:
    neigh_heatmap = plotly_heatmap(
      adata_g,
      uns_key="cluster_nhood_enrichment",
      title=f"{neigh_group_by.value} cells: Neighborhood Enrichment",
      method=clustering_method.value,
      mode=mode.value,
      vmax=vmax,
      vmin=vmin
    )
  except ValueError:
    w_text_output(
      content=f"No neighborhoods found in object; computing neighborhoods...",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    adata_g = squidpy_analysis(adata_g)
    neigh_heatmap = plotly_heatmap(
      adata_g,
      uns_key="cluster_nhood_enrichment",
      title=f"{neigh_group_by.value} cells: Neighborhood Enrichment",
      method=clustering_method.value,
      mode=mode.value,
      vmax=vmax,
      vmin=vmin
    )

elif neigh_group_by.value in ["sample", "condition"]:

  group = neigh_group_by.value
  sub_groups = group_dict[group]
  if group not in filtered_groups:
    filtered_adatas: dict[str, anndata.AnnData] = {}

    for sg in sub_groups:
      w_text_output(
        content=f"Computing spatial neighbrohoods for {sg}...",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      filtered_adata = filter_anndata(adata, group, sg)
      squidpy_analysis(filtered_adata)

      filtered_adatas[sg] = filtered_adata

    filtered_groups[group] = filtered_adatas

  neigh_heatmap = plot_neighborhood_groups(
    filtered_groups[group],
    f"Neighborhoods by {group}",
    uns_key="cluster_nhood_enrichment",
    method=clustering_method.value,
    mode=mode.value,
    vmax=vmax,
    vmin=vmin
  )


else:
  raise KeyError("Group by not expected value")
