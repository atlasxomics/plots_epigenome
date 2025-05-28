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

ghm_group = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select grouping for heatmap."
  }
)

ghm_group = ghm_group.value
if ghm_group == "condition":
  ghm_group = "conditions1"

genes_heatmap_df = adata_g.uns[f"genes_per_{ghm_group}_hm"]
if "Unnamed: 0" in genes_heatmap_df.columns:
  genes_heatmap_df = genes_heatmap_df.set_index("Unnamed: 0")

title="Gene Score Heatmap"

gene_heatmap = px.imshow(
  genes_heatmap_df,
  color_continuous_scale='Spectral_r',
  aspect='auto',
  origin='lower'
)

gene_heatmap.update_layout(
    title=title,
    xaxis_title=ghm_group,
    yaxis_title="Genes",
    coloraxis_colorbar=dict(title="log2 fold change"),
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
