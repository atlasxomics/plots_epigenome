from typing import Any
import numpy as np
import plotly.graph_objects as go
from anndata import AnnData
import scipy.cluster.hierarchy as sch
from plotly.subplots import make_subplots
import pandas as pd


# --------------------------------------------------------------------------------
def filter_anndata(
    adata: anndata.AnnData, group: str, subgroup: List[str]
) -> anndata.AnnData:
    return adata[adata.obs[group] == subgroup]

def plotly_heatmap(
  adata: AnnData,
  key: str = "cluster",
  title: str = "",
  method: str = "None",
  colorscale: str = "RdBu_r",
  width: Optional[int] = 800,
  height: Optional[int] = 800,
  uns_key: Optional[str] = None,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
  **kwargs: Any,
) -> go.Figure:
  """
  Create an interactive heatmap using Plotly.
  
  Parameters
  ----------
  adata : AnnData
      Annotated data matrix
  key : str
      Key for observations in adata.obs
  title : str
      Title of the plot
  method : Optional[str]
      Clustering method for dendrograms; if "None" sorts
      numerically. Options:
       - single
       - complete
       - average
       - weighted
       - centroid
       - median
       - ward
  colorscale : str
      Plotly colorscale name
  width : Optional[int]
      Figure width in pixels (default: 800)
  height : Optional[int]
      Figure height in pixels (default: 800)
  uns_key : Optional[str]
      Key in adata.uns to retrieve data from
  mode : str
      Mode for data retrieval (default: 'zscore')
  vmin : Optional[float]
      Minimum value for color scale. If None, uses data minimum
  vmax : Optional[float]
      Maximum value for color scale. If None, uses data maximum
  **kwargs : Any
      Additional kwargs passed to go.Heatmap
      
  Returns
  -------
  go.Figure
      Plotly figure object
  """
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
      textfont={"size": 10},
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

from typing import Dict


def plot_neighborhood_groups(
  group_adatas: Dict["str", anndata.AnnData],
  key: str = "cluster",
  method: str = "None",
  uns_key: Optional[str] = None,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
):
  groups = list(group_adatas.keys())
  
  # Calculate number of rows needed
  num_rows = (len(groups) - 1) // 2 + 1
  num_cols = min(len(groups), 2)
  
  # Create subplots with calculated rows and columns
  combined_fig = make_subplots(
    rows=num_rows, cols=num_cols, subplot_titles=groups, horizontal_spacing=0.05, vertical_spacing=0.1
  )

  # Loop through each sample
  for i, group in enumerate(groups):
    row = i // 2 + 1
    col = i % 2 + 1
    
    # Create heatmap plot for the subset
    fig = plotly_heatmap(
      group_adatas[group],
      uns_key="cluster_nhood_enrichment",
      title=f"{group}: Neighborhood Enrichment",
      method=method,
      mode=mode,
      vmax = vmax,
      vmin = vmin
    )
    
    # Add each trace to the combined figure
    for trace in fig.data:
      combined_fig.add_trace(trace, row=row, col=col)

  subplot_width = 500
  subplot_height = 500

  combined_fig.update_layout(
      plot_bgcolor='rgba(0,0,0,0)',
      autosize=True,
      margin=dict(
        pad=10
      )
    )

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
    content="No data loaded...",
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
      vmax = vmax,
      vmin = vmin
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
      vmax = vmax,
      vmin = vmin
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
    uns_key="cluster_nhood_enrichment",
    method=clustering_method.value,
    mode=mode.value,
    vmax = vmax,
    vmin = vmin
  )


else:
  raise KeyError("Group by not expected value")
