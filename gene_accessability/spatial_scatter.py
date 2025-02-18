w_text_output(content=f"""

# Spatial Scatter Plot (Gene Accessibility)

<details>
<summary><i>details</i></summary>
Display all cells from the experiment arranged spatially or by UMAP coordinates; color cells by the gene activity score. 

The gene activity score is calculated via the SnapATAC2 funcition [`make_gene_matrix`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.make_gene_matrix.html).  For multiple genes, the scanpy function [`score_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html) is used. 

> WARNING: Selecting multiple genes may take a minute or more to display.  We are actively working to speed-up this step.

</details>

""")

if not adata_g:
  w_text_output(
    content="No data gene activity data loaded...",
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
  default="1",
  options=tuple(str(i) for i in np.linspace(0.5, 7, 14)),
  appearance={
    "help_text": "Select the size of the displayed points."
  }
)


w_row(items=[genes, genes_coords, genes_pt_size])

genes_signal_value = genes._signal.sample()

# Check if genes has a value.
if genes_signal_value.__class__.__name__ in ["Nothing", "NoneType"]:
  w_text_output(
    content="Please select genes for plotting.",
    appearance={"message_box": "info"}
  )
  submit_widget_state()
  exit(0)
  
submit_widget_state()
w_text_output(
  content=f"Adding gene score for {genes_signal_value} to .obs",
  appearance={"message_box": "info"}
)


if len(genes_signal_value) == 1:
  name = genes_signal_value[0]
  convert_feature_expression(adata_g, name)
elif len(genes_signal_value) > 1:
  sc.tl.score_genes(adata_g, genes_signal_value,  use_raw=False)
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
    marker_size=float(genes_pt_size.value)
  )
  fig_genes.update_coloraxes(colorscale='Spectral_r')
elif genes_coords.value == "spatial":
  fig_genes = plot_umap_for_samples(
    adata_g,
    samples,
    color_by=name,
    pt_size=float(genes_pt_size.value),
    coords=genes_coords.value,
  )

print(fig_genes)