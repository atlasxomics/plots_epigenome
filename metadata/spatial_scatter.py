w_text_output(content="""

# Spatial Scatter Plot (Cell Data)

<details>
<summary><i>details</i></summary>
Display all cells from the experiment arranged spatially or by UMAP coordinates; color cells by Categorical or numberic metadata values.  Available values correspond to the .obs column names in the AnnData Object.
</details>

""")

new_data_signal()

if not adata:
  w_text_output(
    content="No data data loaded...",
    appearance={"message_box": "warning"}
  )
  exit()

meta_color_by = w_select(
  label="metadata",
  default="cluster",
  options=available_metadata,
  appearance={
    "help_text": "Select the data to color points/cells by."
  }
)
meta_coords = w_select(
  label="plot coordinates",
  default="spatial",
  options=obsm_keys,
  appearance={
    "help_text": "Select how to arrange points/cells."
  }
  
)
meta_pt_size = w_select(
  label="point size",
  default="2.5",
  options=tuple(str(i) for i in np.linspace(0.5, 7, 14)),
 appearance={
    "help_text": "Select the size of the displayed points."
  }
)

meta_color = w_select(
  label="color scheme",
  default="bright",
  options=all_colors,
  appearance={
    "help_text": " "
  }
)

meta_highlight = w_select(
    label="highlight spatial cluster",
    default=None,
    options=tuple(clusters + ["None"]),
    appearance={
      "help_text": "Highlight cells belonging to a specific cluster for spatial coordinates."
    }
)

meta_max = w_text_input(
  label="color scale max",
  default="",
  appearance={
    "help_text": "Set color scale maximum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

meta_min = w_text_input(
  label="color scale min",
  default="",
  appearance={
    "help_text": "Set color scale minimum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

w_row(items=[meta_color_by, meta_coords, meta_pt_size, meta_color, meta_highlight, meta_min, meta_max])

meta_flipy = w_checkbox(
  label="flip y",
  default=True,
  appearance={
    "description": "Flip vertical axis."
  }
)

meta_button = w_button(label="Update Spatial Scatter Plot")

if meta_color_by.value is not None and meta_button.value:
  
  meta_min_val = safe_float(meta_min.value, "Cannot convert legend min to float; ignoring...")
  meta_max_val = safe_float(meta_max.value, "Cannot convert legend max to float; ignoring...")

  if meta_min_val is not None and meta_max_val is not None:
      if meta_max_val <= meta_min_val:
          w_text_output(
              content="Legend max is less than or equal to min; ignoring...",
              appearance={"message_box": "warning"}
          )
  
  if meta_coords.value == "X_umap":
  
    # Check if color_by is discrete or continuous
    obs_groups = sorted(adata.obs[meta_color_by.value].unique())
    is_discrete = len(obs_groups) < 30
  
    temp_fig = snap.pl.umap(
      adata,
      use_rep=meta_coords.value,
      show=False,
      color=meta_color_by.value,
      marker_size=float(meta_pt_size.value),
      width=1300,
      height=800
    )
  
    if is_discrete:
      meta_fig = custom_plotly(
        temp_fig,
        color_scheme=meta_color.value,
        width=1300,
        height=800
  
      )
    else:
      obs_values = adata.obs[meta_color_by.value]
      meta_min_val = meta_min_val if meta_min_val is not None else obs_values.min()
      meta_max_val = meta_max_val if meta_max_val is not None else obs_values.max()

      if meta_max_val <= meta_min_val:
          w_text_output(
              content="Legend max is less than or equal to min; ignoring...",
              appearance={"message_box": "warning"}
          )

      meta_fig = temp_fig.update_coloraxes(
          colorscale='Spectral_r',
          colorbar_title=meta_color_by.value,
          cmin=meta_min_val,
          cmax=meta_max_val,
        )
  
    temp_fig = None
  
  elif meta_coords.value == "spatial":
  
    meta_fig = plot_umap_for_samples(
      adata,
      samples,
      color_by=meta_color_by.value,
      pt_size=float(meta_pt_size.value),
      coords=meta_coords.value,
      color_scheme=meta_color.value,
      flipY=meta_flipy.value,
      show_cluster=meta_highlight.value if meta_highlight.value != "None" else None,
      vmin=meta_min_val,
      vmax=meta_max_val,
    )
  
  w_plot(source=meta_fig)