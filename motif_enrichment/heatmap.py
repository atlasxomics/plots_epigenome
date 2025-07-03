w_text_output(content="""

# Heatmap (Motif Deviation)

<details>
<summary><i>details</i></summary>
Plot a heatmap for differentially regulated motifs.
</details>

""")

new_data_signal()

if not adata_m:
    w_text_output(
        content="No motif deviation data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif deviation data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

mhm_group = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select grouping for heatmap."
  }
)

mhm_button = w_button(label="Update Heatmap")

if mhm_group.value is not None and mhm_button.value:

  mhm_group = mhm_group.value

  # Handle "condition" case to support ArchR or old Snap
  if mhm_group == "condition":

      mhm_group = "conditions1"

      possible_keys = [
          "motif_per_conditions1_hm",
          "motif_per_condition_1_hm"
      ]
      for key in possible_keys:
          if key in adata_m.uns:
              motifs_heatmap_df = adata_m.uns[key]
              break
      else:
          w_text_output(
              content=f"No motifs heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
  else:
      key = f"motif_per_{mhm_group}_hm"
      if key not in adata_m.uns:
          w_text_output(
              content=f"No motifs heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
      motifs_heatmap_df = adata_m.uns[key]

  if "Unnamed: 0" in motifs_heatmap_df.columns:
    motifs_heatmap_df = motifs_heatmap_df.set_index("Unnamed: 0")
  
  title="Motif Score Heatmap"
  
  motifs_heatmap = px.imshow(
    motifs_heatmap_df,
    color_continuous_scale='Spectral_r',
    aspect='auto',
    origin='lower'
  )
  
  motifs_heatmap.update_layout(
      title=title,
      xaxis_title=mhm_group,
      yaxis_title="Motifs",
      coloraxis_colorbar=dict(title="Log2FC")
  )
  
  motifs_heatmap.update_xaxes(
    side="bottom",
    tickmode='array',
    tickvals=list(range(len(motifs_heatmap_df.columns))),
    ticktext=motifs_heatmap_df.columns,
    tickangle=45
  )
  motifs_heatmap.update_yaxes(autorange="reversed")
  
  w_plot(source=motifs_heatmap)
