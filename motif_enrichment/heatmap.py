w_text_output(content="""

# Motif Deviation Heatmap

<details>
<summary><i>details</i></summary>

Plot a heatmap for differential motifs.

<br>

The plot can display either the top n motifs, ranked by a user-defined parameter (values to plot), or a manually selected list of motifs. The selections for the values to plot parameters are set by the scanpy function make_gene_matrix; see the Returns section of the linked documentation for a description of the values.

<br>

Manually inputting motifs will overwrite the 'top n motifs' parameter.</details>

""")

if not adata_m:
  w_text_output(content="No motif data loaded...",  appearance={"message_box": "warning"})
  exit()

hm_values = ("scores", "pvals", "pvals_adj")

n_motifs = w_text_input(
  label="top n motifs",
  default="5",
  appearance={
    "help_text": "Set number of motifs to plot for each group."
  }
)
hm_value = w_select(
  label="value to plot",
  default="scores",
  options=hm_values,
  appearance={
    "help_text": "Select values to color the heatmap by."
  }
)
hm_group = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select grouping for heatmap."
  }
)

hm_motifs = w_multi_select(
  label="select motifs",
  default=None,
  options=available_motifs,
  appearance={
    "detail": "optional",
    "help_text": "Select a list of motifs to display; selecting motifs will overwrite the `top n motifs` parameter."
  }
)

w_row(items=[n_motifs, hm_value, hm_group, hm_motifs])

motifs_rank_df = sc.get.rank_genes_groups_df(
  adata_m,
  group=None,
  key=f"{hm_group.value}_motifs"
)

title="Motif Deviation Heatmap"
if hm_motifs.value == None or len(hm_motifs.value) == 0:
  motif_heatmap_df = get_top_n_heatmap(
    motifs_rank_df, hm_value.value, int(n_motifs.value)
  )
  title = f"Top {n_motifs.value} motifs by {hm_value.value} for each {hm_group.value}"
elif len(hm_motifs.value) > 0:
  motif_heatmap_df = get_feature_heatmap(
    motifs_rank_df, hm_motifs.value, hm_value.value
  )
  title = f"User-defined motifs by {hm_value.value} for each {hm_group.value}"

fig_11503 = px.imshow(
  motif_heatmap_df,
  color_continuous_scale='Spectral_r',
  aspect='auto',
  origin='lower'
)

fig_11503.update_layout(
    title=title,
    xaxis_title="Motifs",
    yaxis_title=hm_group.value,
    coloraxis_colorbar_title=hm_value.value
)

fig_11503.update_xaxes(side="bottom", tickmode='array', tickvals=list(range(len(motif_heatmap_df.columns))), ticktext=motif_heatmap_df.columns, tickangle=45)
fig_11503.update_yaxes(autorange="reversed")

fig_11503.show()
