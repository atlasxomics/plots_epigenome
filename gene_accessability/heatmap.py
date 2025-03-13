w_text_output(content="""

# Gene Activity Heatmap

<details>
<summary><i>details</i></summary>
Plot a heatmap for differentially regulated genes.

The plot can display either the top n genes, ranked by a user-defined parameter (`values to plot`), or a manually selected list of genes.  The selections for the `values to plot` parameters are set by the scanpy function [make_gene_matrix](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html); see the Returns section of the linked documentation for a description of the values. 

> Manually inputting genes will overwrite the `top n genes` parameter.
</details>

""")

if not adata_g:
  w_text_output(
    content="No data gene activity data loaded...",
    appearance={"message_box": "warning"}
  )
  exit()

hm_values = ("scores", "logfoldchanges", "pvals", "pvals_adj")

n_genes = w_text_input(
  label="top n genes",
  default="5",
  appearance={
  "help_text": "Set number of genes to plot for each group."
  }
)
ghm_value = w_select(
  label="value to plot",
  default="logfoldchanges",
  options=hm_values,
  appearance={
    "help_text": "Select values to color the heatmap by."
  }
)
ghm_group = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select grouping for heatmap."
  }
)

hm_genes = w_multi_select(
  label="select genes",
  default=None,
  options=available_genes,
  appearance={
    "detail": "optional",
    "help_text": "Select a list of genes to display; selecting genes will overwrite the `top n genes` parameter."
  }
)

w_row(items=[n_genes, ghm_value, ghm_group, hm_genes])

ghm_displayGm = w_checkbox(
  label="display pseudogenes (mouse only)",
  default=True,
  appearance={
    "description": "Whether to display pseudogenes; for mouse data only."
  }
)

rank_df = sc.get.rank_genes_groups_df(
  adata_g,
  group=None,
  key=f"{ghm_group.value}_genes"
)

title="Gene Score Heatmap"
if hm_genes.value == None or len(hm_genes.value) == 0:
  genes_heatmap_df = get_top_n_heatmap(
    rank_df, ghm_value.value, int(n_genes.value)
  )
  title = f"Top {n_genes.value} genes by {ghm_value.value} for each {ghm_group.value}"
elif len(hm_genes.value) > 0:
  genes_heatmap_df = get_feature_heatmap(
    rank_df, hm_genes.value, ghm_value.value
  )
  title = f"User-defined genes by {ghm_value.value} for each {ghm_group.value}"

if not ghm_displayGm.value:
    genes_heatmap_df = genes_heatmap_df.loc[:, ~genes_heatmap_df.columns.str.startswith("Gm")]

gene_heatmap = px.imshow(
  genes_heatmap_df,
  color_continuous_scale='Spectral_r',
  aspect='auto',
  origin='lower'
)

gene_heatmap.update_layout(
    title=title,
    xaxis_title="Genes",
    yaxis_title=ghm_group.value,
    coloraxis_colorbar_title=ghm_value.value
)

gene_heatmap.update_xaxes(
  side="bottom",
  tickmode='array',
  tickvals=list(range(len(genes_heatmap_df.columns))),
  ticktext=genes_heatmap_df.columns,
  tickangle=45
)
gene_heatmap.update_yaxes(autorange="reversed")

gene_heatmap.show()
