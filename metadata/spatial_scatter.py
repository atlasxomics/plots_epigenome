w_text_output(content=f"""

# Spatial Scatter Plot (Cell Data)

<details>
<summary><i>details</i></summary>
Display all cells from the experiment arranged spatially or by UMAP coordinates; color cells by Categorical or numberic metadata values.  Available values correspond to the .obs column names in the AnnData Object.
</details>

""")


if not adata:
  w_text_output(
    content="No data gene activity data loaded...",
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

w_row(items=[meta_color_by, meta_coords, meta_pt_size, meta_color])

if meta_coords.value == "X_umap":
  temp_fig = snap.pl.umap(
    adata,
    use_rep=meta_coords.value,
    show=False,
    color=meta_color_by.value,
    marker_size=float(meta_pt_size.value),
  )
  meta_fig = custom_plotly(
    temp_fig,
    color_scheme=meta_color.value,
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
  )

meta_fig
