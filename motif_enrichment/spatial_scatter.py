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

w_row(items=[motifs, motif_coords, motif_pt_size])

motifs_signal_value = motifs._signal.sample()

# Check if genes has a value.
if motifs_signal_value.__class__.__name__ in ["Nothing", "NoneType"]:
  w_text_output(
    content="Please select motifs for plotting.",
    appearance={"message_box": "info"}
  )
  submit_widget_state()
  exit(0)

submit_widget_state()
w_text_output(
    content=f"Adding motif deviation score for {motifs_signal_value} to .obs",
    appearance={"message_box": "info"}
)

try:
  if len(motifs_signal_value) == 1:
      name = motifs_signal_value[0]
      convert_feature_expression(adata_m, name)
  elif len(motifs_signal_value) > 1:
      sc.tl.score_genes(adata_m, motifs_signal_value, use_raw=False)
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
        marker_size=float(motif_pt_size.value)
    )
    fig_motifs.update_coloraxes(colorscale='Spectral_r')
elif motif_coords.value == "spatial":
    fig_motifs = plot_umap_for_samples(
        adata_m,
        samples,
        color_by=name,
        pt_size=float(motif_pt_size.value),
        coords=motif_coords.value
    )

print(fig_motifs)
