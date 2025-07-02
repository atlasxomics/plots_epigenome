w_text_output(content="""

# Heatmap (Gene Accessibility)

<details>
<summary><i>details</i></summary>
Plot a heatmap for differentially regulated genes.

The plot can display either the top n genes, ranked by a user-defined parameter (`values to plot`), or a manually selected list of genes.  The selections for the `values to plot` parameters are set by the scanpy function [make_gene_matrix](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html); see the Returns section of the linked documentation for a description of the values. 

> Manually inputting genes will overwrite the `top n genes` parameter.
</details>

""")

new_data_signal()

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

ghm_button = w_button(label="Update Heatmap")

if ghm_group.value is not None and ghm_button.value:

  ghm_group = ghm_group.value

  # Handle "condition" case to support ArchR or old Snap
  if ghm_group == "condition":

      ghm_group = "conditions1"

      possible_keys = [
          "genes_per_conditions1_hm",
          "genes_per_condition_1_hm"
      ]
      for key in possible_keys:
          if key in adata_g.uns:
              genes_heatmap_df = adata_g.uns[key]
              break
      else:
          w_text_output(
              content=f"No genes heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
  else:
      key = f"genes_per_{ghm_group}_hm"
      if key not in adata_g.uns:
          w_text_output(
              content=f"No genes heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
      genes_heatmap_df = adata_g.uns[key]

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
      coloraxis_colorbar=dict(title="Log2FC")
  )

  gene_heatmap.update_xaxes(
    side="bottom",
    tickmode='array',
    tickvals=list(range(len(genes_heatmap_df.columns))),
    ticktext=genes_heatmap_df.columns,
    tickangle=45
  )
  gene_heatmap.update_yaxes(autorange="reversed")
  
  w_plot(source=gene_heatmap)