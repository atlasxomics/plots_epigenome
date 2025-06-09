w_text_output(content="""

# Spatial Scatter Plot (Cell Data)

<details>
<summary><i>details</i></summary>
Display all cells from the experiment arranged spatially or by UMAP coordinates; color cells by Categorical or numberic metadata values.  Available values correspond to the .obs column names in the AnnData Object.
</details>

""")


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
    "help_text": "Set color scale maximum for spatial plots."
  }
)

meta_min = w_text_input(
  label="color scale min",
  default="",
  appearance={
    "help_text": "Set color scale minimum for spatial plots."
  }
)

w_row(items=[meta_color_by, meta_coords, meta_pt_size, meta_color, meta_highlight, meta_max, meta_min])

meta_flipy = w_checkbox(
  label="flip y",
  default=False,
  appearance={
    "description": "Flip vertical axis."
  }
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
    meta_fig = temp_fig.update_coloraxes(
        colorscale='Spectral_r',
        colorbar_title=meta_color_by.value,
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
    vmin=float(meta_min.value) if meta_min.value != "" else None,
    vmax=float(meta_max.value) if meta_max.value != "" else None,
  )

meta_fig
