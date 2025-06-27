w_text_output(content="""

#  Spatial Scatter (Motif Deviation)

<details>
<summary><i>details</i></summary>

In this section, you can visualize differential accessibility of motifs.  The motif deviation matrix was created via the [pychromvar](https://pychromvar.readthedocs.io/en/stable/index.html), a Python implementation of the R package [chromvar](https://www.nature.com/articles/nmeth.4401); please see the publication for more information. 

</details>

""")

if not adata_m:
    w_text_output(
       content="No motif data selected...",
       appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

motifs = w_multi_select(
    label="select motifs",
    options=available_motifs,
    appearance={
      "placeholder": "MA0143.5.SOX2",
      "help_text": "You must click 'Run' after selecting motifs to run Cell.",
      "description": "You must click 'Run' after selecting motifs to run Cell; selecting multiple motifs may take a minute or more to display."
    }
)

motif_coords = w_select(
    label="plot coordinates",
    default="spatial",
    options=obsm_keys,
    appearance={
      "help_text": "Select how to arrange points/cells."
    }
)

motif_pt_size = w_select(
    label="point size",
    default="2.5",
    options=tuple(str(i) for i in np.linspace(0.5, 7, 14)),
    appearance={
      "help_text": "Select the size of the displayed points."
    }
)

motif_highlight = w_select(
    label="highlight spatial cluster",
    default=None,
    options=tuple(clusters + ["None"]),
    appearance={
      "help_text": "Highlight cells belonging to a specific cluster for spatial coordinates."
    }
)

motif_max = w_text_input(
  label="color scale max",
  default="",
  appearance={
    "help_text": "Set color scale maximum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

motif_min = w_text_input(
  label="color scale min",
  default="",
  appearance={
    "help_text": "Set color scale minimum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

w_row(items=[motifs, motif_coords, motif_pt_size, motif_highlight, motif_min, motif_max])

motif_flipy = w_checkbox(
  label="flip y",
  default=True,
  appearance={
    "description": "Flip vertical axis."
  }
)

motifs_button = w_button(label="Update Spatial Scatter")

if motifs.value is not None and motifs_button.value:
  
  if motif_min.value != "" and motif_max.value != "":
      try:
          motif_min_val = float(motif_min.value)
          motif_max_val = float(motif_max.value)
          if motif_max_val <= motif_min_val:
              w_text_output(
                  content="Legend max is less than or equal to min; ignoring...",
                  appearance={"message_box": "warning"}
              )
      except (TypeError, ValueError):
          w_text_output(
              content="Cannot convert legend min or max to float; ignoring...",
              appearance={"message_box": "warning"}
          )
          motif_min_val = ""
          motif_max_val = ""
  
  # Check if genes has a value.
  if motifs._signal.sample().__class__.__name__ in ["Nothing", "NoneType"]:
    w_text_output(
      content="Please select motifs for plotting.",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    exit(0)
  
  w_text_output(
      content=f"Adding motif deviation score for {motifs._signal.sample()} to .obs",
      appearance={"message_box": "info"}
  )
  submit_widget_state()
  
  try:
    if len(motifs._signal.sample()) == 1:
        name = motifs._signal.sample()[0]
        convert_feature_expression(adata_m, name)
    elif len(motifs._signal.sample()) > 1:
        sc.tl.score_genes(adata_m, motifs._signal.sample(), use_raw=False)
        name = "score"
    else:
        print("No motifs detected")
        exit()
  except TypeError:
    w_text_output(
      content="Please select motif(s) to plot.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
  
  
  if motif_coords.value == "X_umap":
      fig_motifs = snap.pl.umap(
          adata_m,
          use_rep=motif_coords.value,
          show=False,
          color=name,
          marker_size=float(motif_pt_size.value),
          width=1300,
          height=800,
      )
      fig_motifs.update_coloraxes(colorscale='Spectral_r')
  elif motif_coords.value == "spatial":
      fig_motifs = plot_umap_for_samples(
          adata_m,
          samples,
          color_by=name,
          pt_size=float(motif_pt_size.value),
          coords=motif_coords.value,
          flipY=motif_flipy.value,
          show_cluster=motif_highlight.value if motif_highlight.value != "None" else None,
          vmin=float(motif_min.value) if motif_min.value != "" else None,
          vmax=float(motif_max.value) if motif_max.value != "" else None,
      )
  
  w_plot(source=fig_motifs)
