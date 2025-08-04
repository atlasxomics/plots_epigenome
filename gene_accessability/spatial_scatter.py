w_text_output(content="""

# Spatial Scatter Plot (Gene Accessibility)

<details>
<summary><i>details</i></summary>
Display all cells from the experiment arranged spatially or by UMAP coordinates; color cells by the gene activity score. 

The gene activity score is calculated via the SnapATAC2 funcition [`make_gene_matrix`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.make_gene_matrix.html).  For multiple genes, the scanpy function [`score_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) is used. 

> WARNING: Selecting multiple genes may take a minute or more to display.  We are actively working to speed-up this step.

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

genes = w_multi_select(
  label="select genes",
  options=available_genes,
  appearance={
    "placeholder": "GAPDH",
    "help_text": "You must click 'Run' after selecting genes to run Cell.",
    "description": "You must click 'Run' after selecting genes to run Cell; selecting multiple genes may take a minute or more to display."
  }
)

genes_coords = w_select(
  label="plot coordinates",
  default="spatial",
  options=obsm_keys,
  appearance={
    "help_text": "Select how to arrange points/cells."
  }
)
genes_pt_size = w_select(
  label="point size",
  default="2.5",
  options=tuple(str(i) for i in np.linspace(0.5, 7, 14)),
  appearance={
    "help_text": "Select the size of the displayed points."
  }
)

gene_highlight = w_select(
    label="highlight spatial cluster",
    default=None,
    options=tuple(clusters + ["None"]),
    appearance={
      "help_text": "Highlight cells belonging to a specific cluster for spatial coordinates."
    }
)

gene_max = w_text_input(
  label="color scale max",
  default="",
  appearance={
    "help_text": "Set color scale maximum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

gene_min = w_text_input(
  label="color scale min",
  default="",
  appearance={
    "help_text": "Set color scale minimum for spatial plots; must set min AND max for custom thresholds to display."
  }
)

w_row(items=[genes, genes_coords, genes_pt_size, gene_highlight, gene_min, gene_max])

genes_flipy = w_checkbox(
  label="flip y",
  default=True,
  appearance={
    "description": "Flip vertical axis."
  }
)

genes_button = w_button(label="Update Spatial Scatter")

if genes.value is not None and genes_button.value:

  gene_min_val = safe_float(gene_min.value, "Cannot convert legend min to float; ignoring...")
  gene_max_val = safe_float(gene_max.value, "Cannot convert legend max to float; ignoring...")

  if gene_min_val is not None and gene_max_val is not None:
      if gene_max_val <= gene_min_val:
          w_text_output(
              content="Legend max is less than or equal to min; ignoring...",
              appearance={"message_box": "warning"}
          )
  
  # Check if genes has a value.
  if genes.value.__class__.__name__ in ["Nothing", "NoneType"]:
      w_text_output(
        content="Please select genes for plotting.",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      exit(0)

  w_text_output(
    content=f"Adding gene score for {genes.value} to .obs",
    appearance={"message_box": "info"}
  )
  submit_widget_state()
  
  if len(genes.value) == 1:
      name = genes.value[0]
      convert_feature_expression(adata_g, name)
  elif len(genes.value) > 1:
      sc.tl.score_genes(adata_g, genes.value,  use_raw=False)
      name = "score"
  else:
      print("No genes detected")
      exit()
  
  print(adata_g.obs[name])

  if genes_coords.value == "X_umap":
      fig_genes = snap.pl.umap(
        adata_g,
        use_rep=genes_coords.value,
        show=False,
        color=name,
        marker_size=float(genes_pt_size.value),
        width=1300,
        height=800,
      )

      gene_values = adata_g.obs[name]
      gene_min_val = gene_min_val if gene_min_val is not None else gene_values.min()
      gene_max_val = gene_max_val if gene_max_val is not None else gene_values.max()

      if gene_max_val <= gene_min_val:
          w_text_output(
              content="Legend max is less than or equal to min; ignoring...",
              appearance={"message_box": "warning"}
          )
      
      fig_genes.update_coloraxes(
        colorscale='Spectral_r', cmin=gene_min_val, cmax=gene_max_val
      )

  elif genes_coords.value == "spatial":
      fig_genes = plot_umap_for_samples(
        adata_g,
        samples,
        color_by=name,
        pt_size=float(genes_pt_size.value),
        coords=genes_coords.value,
        flipY=genes_flipy.value,
        show_cluster=gene_highlight.value if gene_highlight.value != "None" else None,
        vmin=gene_min_val,
        vmax=gene_max_val,
      )
  
  w_plot(source=fig_genes)