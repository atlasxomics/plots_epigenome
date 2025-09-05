new_data_signal()

w_text_output(content="""

# Neighborhood Analysis

Explore spatial neighborhood enrichment among clusters, either for **all cells** or split by a subgroup (e.g., **sample** or **condition**).  See squipy [neighbors enrichment analysis](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_nhood_enrichment.html) for more information.


<details>
<summary><i>details</i></summary>

Each heatmap cell reflects how often cells from **cluster A** neighbor cells from **cluster B** compared with chance.  You can view values as **z-scores** (standardized enrichment; recommended) or **counts** (raw neighborhood counts).  Optionally, you can use hierarchical clustering to reorder rows/columns to reveal structure.

### Controls

1. **subplot groups** 
   - Options: **all**, **sample**, **condition**  
   - **all**: one heatmap using all cells.
   - **sample / condition**: one heatmap per subgroup (faceted).

2. **displayed data**
   - Options: **zscore**, **count**  
   - **zscore**: standardized neighborhood enrichment (best for comparisons). 
   - **count**: raw neighbor counts (scale depends on dataset size).

3. **hierarchical clustering method** 
   - Options: **None**, **single**, **complete**, **average**, **weighted**, **centroid**, **median**, **ward**  
   - Choose **None** to keep the original cluster order, or a method to cluster rows/columns and group similar patterns.

4. **colorscale maximum / minimum** 
   - Optional numeric limits for the heatmap color range (e.g., max = `5`, min = `-2`).  
   - Leave blank to auto-scale.
</details>

""")

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

neigh_button = w_button(label="Update Neighborhood Plots")

if neigh_group_by.value is not None and neigh_button.value:
  
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
  
  w_plot(source=neigh_heatmap)