new_data_signal()

w_text_output(content="""

# Neighborhood Analysis

Explore spatial neighborhood enrichment among clusters, either for **all cells** or split by a subgroup (e.g., **sample** or **condition**).  See squipy [neighbors enrichment analysis](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_nhood_enrichment.html) for more information.


<details>
<summary><i>details</i></summary>

Each heatmap cell reflects how often cells from **cluster A** neighbor cells from **cluster B** compared with chance.  You can view values as **z-scores** (standardized enrichment; recommended) or **counts** (raw neighborhood counts).  Optionally, you can facet plots by categorical observations, including custom annotations.

### Controls

1. **display type**
   - Options: **heatmap**, **radial**
   - **heatmap**: the matrix view (this cell).
   - **radial**: the "spoke" view (shown in the cell below).

2. **subplot groups** 
   - Options: **all** plus categorical observations (for example **sample**, **condition**, or custom annotations)
   - **all**: one heatmap using all cells.
   - **categorical observation**: one heatmap per subgroup (faceted).
   - Workflow groupings use precomputed results. Custom annotations are computed when first displayed and may take longer.

3. **displayed data**
   - Options: **zscore**, **count**  
   - **zscore**: standardized neighborhood enrichment (best for comparisons). 
   - **count**: raw neighbor counts (scale depends on dataset size).

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


def load_precomputed_neighborhood_groups(adata):
  """Return lightweight AnnData objects from workflow-precomputed results."""
  root = adata.uns.get("cluster_nhood_enrichment_by_group")
  if not isinstance(root, dict) or root.get("schema_version") != 1:
    return {}

  result = {}
  for group_entry in root.get("groups", {}).values():
    group_key = str(group_entry["group_key"])
    subgroups = {}
    for subgroup_entry in group_entry.get("subgroups", {}).values():
      group_value = str(subgroup_entry["group_value"])
      categories = pd.Index(
        np.asarray(subgroup_entry["cluster_categories"]).astype(str)
      )
      obs = pd.DataFrame({
        "cluster": pd.Categorical(
          categories,
          categories=categories,
          ordered=True,
        )
      })
      neighborhood_adata = anndata.AnnData(obs=obs)
      neighborhood_adata.uns["cluster_nhood_enrichment"] = {
        "zscore": np.asarray(subgroup_entry["zscore"]),
        "count": np.asarray(subgroup_entry["count"]),
      }
      subgroups[group_value] = neighborhood_adata
    result[group_key] = subgroups

  return result


def make_lightweight_neighborhood_adata(adata, group, subgroup):
  """Subset only metadata and coordinates, never the backed feature matrix."""
  if "spatial_offset" not in adata.obsm:
    raise KeyError(
      "Expected `spatial_offset` in adata.obsm (created in the Select Data step)."
    )

  mask = (adata.obs[group] == subgroup).to_numpy()
  obs_keys = ["cluster"]
  if "sample" in adata.obs.columns:
    obs_keys.append("sample")

  subset_obs = adata.obs.loc[mask, obs_keys].copy()
  for obs_key in obs_keys:
    if pd.api.types.is_categorical_dtype(subset_obs[obs_key]):
      subset_obs[obs_key] = subset_obs[obs_key].cat.remove_unused_categories()

  lightweight_adata = anndata.AnnData(obs=subset_obs)
  lightweight_adata.obsm["spatial_offset"] = np.asarray(
    adata.obsm["spatial_offset"]
  )[mask].copy()
  return lightweight_adata


precomputed_neighborhoods = load_precomputed_neighborhood_groups(adata_g)

neighbor_groups = [
  key for key in adata_g.obs_keys()
  if key != "cluster" and key not in na_keys and (
    pd.api.types.is_object_dtype(adata_g.obs[key]) or
    pd.api.types.is_categorical_dtype(adata_g.obs[key])
  )
]
group_dict = {g: adata_g.obs[g].dropna().unique() for g in neighbor_groups}

# Display toggle: choose the heatmap (this cell) or the radial 'spoke' view
# (the cell below). Both cells read `neigh_display`; only the selected one
# renders its controls and plot.
neigh_display = w_select(
  label="display type",
  key="neigh_display",
  default="heatmap",
  options=("heatmap", "radial"),
  appearance={
    "help_text": "Switch between the heatmap and radial 'spoke' views of neighborhood enrichment."
  }
)

# When 'radial' is selected, hide the heatmap controls/plot; the radial cell
# below renders instead. The selector stays visible so you can switch back.
if neigh_display.value != "heatmap":
  exit()

neigh_group_by = w_select(
  label="subplot groups",
  default="all",
  options=tuple(["all"] + list(group_dict.keys())),
  appearance={
    "detail": "(all, categorical observation)",
    "help_text": "Facet neighborhood plots by any categorical observation."
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
  "help_text": "Minimum value of colorscale"
  }
)

w_row(items=[neigh_group_by, mode, scale_max, scale_min])

if (
  neigh_group_by.value not in (None, "all")
  and neigh_group_by.value not in precomputed_neighborhoods
):
  w_text_output(
    content=(
      "This custom annotation was added after the workflow ran. Its spatial "
      "neighborhoods will be computed when first displayed and may take longer."
    ),
    appearance={"message_box": "warning"}
  )

neigh_button = w_button(label="Update Neighborhood Plots")

if neigh_group_by.value is not None and neigh_button.value:
  vmax = int(scale_max.value) if scale_max.value and scale_max.value.strip().isdigit() else None
  try:  # Handle negative values
    vmin = int(scale_min.value) if scale_min.value.strip() != '' else None
  except ValueError:
    vmin = None  # Fallback if the value can't be converted to an integer
  
  # --------------------------------------------------------------------------------
  
  if neigh_group_by.value == "all":
    sample_key = "sample" if "sample" in groups else None
    if "cluster_nhood_enrichment" not in adata_g.uns:
      w_text_output(
        content="Computing neighborhoods for all cells...",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      squidpy_analysis(adata_g, sample_key=sample_key)
    else:
      w_text_output(
        content="Using existing neighborhood enrichment for all cells...",
        appearance={"message_box": "info"}
      )
      submit_widget_state()

    neigh_heatmap, neigh_data = plotly_heatmap(
      adata_g,
      uns_key="cluster_nhood_enrichment",
      title=f"{neigh_group_by.value} cells: Neighborhood Enrichment",
      mode=mode.value,
      vmax=vmax,
      vmin=vmin
    )

    neigh_data = pd.DataFrame(neigh_data)
  
  
  elif neigh_group_by.value in group_dict:
  
    group = neigh_group_by.value
    sub_groups = group_dict[group]
    if group in precomputed_neighborhoods:
      filtered_groups[group] = precomputed_neighborhoods[group]
    elif group not in filtered_groups:
      filtered_adatas: dict[str, anndata.AnnData] = {}

      filtered_groups[group] = filtered_adatas

    filtered_adatas = filtered_groups[group]

    for sg in sub_groups:
      if sg not in filtered_adatas:
        filtered_adata = make_lightweight_neighborhood_adata(
          adata_g, group, sg
        )
        filtered_adatas[sg] = filtered_adata

      filtered_adata = filtered_adatas[sg]
      if "cluster_nhood_enrichment" not in filtered_adata.uns:
        w_text_output(
          content=f"Computing spatial neighborhoods for {sg}...",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        sample_key = "sample" if "sample" in filtered_adata.obs else None
        squidpy_analysis(filtered_adata, sample_key=sample_key)
      else:
        w_text_output(
          content=f"Using existing neighborhood enrichment for {sg}...",
          appearance={"message_box": "info"}
        )
        submit_widget_state()

    neigh_heatmap, neigh_data = plot_neighborhood_groups(
      filtered_groups[group],
      f"Neighborhoods by {group}",
      uns_key="cluster_nhood_enrichment",
      mode=mode.value,
      vmax=vmax,
      vmin=vmin
    )

  else:
    raise KeyError("Group by not expected value")
  
  w_plot(source=neigh_heatmap)
  w_table(source=neigh_data)
