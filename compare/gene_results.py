new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

wf_results_signal()
if wf_results_signal.sample() == True:

  if "gene" in results_dict.keys():
    
    w_text_output(content="""## Differential Gene Accessibility""")
    
    g_group = w_select(
        label="Group",
        default=results_dict["gene"]["group_name"].unique()[0],
        options=tuple(results_dict["gene"]["group_name"].unique()),
        key="g_group"
    )
    
    g_pvals_adj_threshold = w_text_input(
      label="pval adjust threshold",
      default="0.05",
      key="g_pvals_adj_threshold"
    )
    
    g_log2fc_threshold = w_text_input(
      label="log2fc threshold",
      default="0.01",
      key="g_log2fc_threshold"
    )
    
    gcompare_rankby = w_select(
        label="Rank By",
        default="Log2FC",
        options=tuple(['Log2FC', 'FDR', 'MeanDiff']),
        key="gcompare_rankby"
    )
    
    gcompare_colorby = w_select(
        label="Color By",
        default="FDR",
        options=tuple(['Log2FC', 'FDR', 'MeanDiff']),
        key="gcompare_colorby"
    )
    
    with w_grid(
      columns=4,
      key="gene_results_controls_grid"
    ) as gene_controls_grid:
      gene_controls_grid.add(item=g_pvals_adj_threshold, col_span=1)
      gene_controls_grid.add(item=g_log2fc_threshold, col_span=1)
      gene_controls_grid.add(item=gcompare_rankby, col_span=1)
      gene_controls_grid.add(item=gcompare_colorby, col_span=1)
    
    # ----------------------------------------------------------------------
    
    if g_group.value is not None:
      
      
      group_g = g_group.value
      df_g = results_dict["gene"]
      df_g = df_g[df_g["group_name"] == group_g]
      
      gvol = plot_volcano(
        df_g,
        float(g_pvals_adj_threshold.value),
        float(g_log2fc_threshold.value),
        "GroupA",
        "GroupB",
        pval_key="FDR",
        l2fc_key="Log2FC",
        names_key="name",
        plot_width=750,
        plot_height=640,
        top_n=2
      )

      grank = plot_ranked_feature_plotly(
          df_g,
          y_col=gcompare_rankby.value,
          x_col=None,
          n_labels=4,
          label_col="name",
          color_col=gcompare_colorby.value,
          colorscale="PuBu_r",
          marker_size=6,
          title="",
          y_label=gcompare_rankby.value
      )

      # Align the two plots on a shared x-axis baseline: give both the same
      # plotting-area height and identical top/bottom margins so, side by side
      # in the row, their x-axes sit at the same vertical position. autosize +
      # width=None lets the row size their widths equally.
      gvol.update_layout(
          height=640,
          autosize=True,
          width=None,
          margin=dict(l=80, r=80, t=60, b=80),
      )
      grank.update_layout(
          height=640,
          autosize=True,
          width=None,
          margin=dict(l=80, r=80, t=60, b=80),
      )

      # Remove the alias left by older versions of this reactive cell. Plot
      # widgets require each figure to have exactly one global variable name.
      globals().pop("_fig", None)

      gvol_plot = w_plot(source=gvol, key="gvol_plot")
      grank_plot = w_plot(source=grank, key='grank_plot')

      with w_grid(
        columns=2,
        key="gene_results_plot_grid"
      ) as gene_results_grid:
        gene_results_grid.add(item=gvol_plot, col_span=1)
        gene_results_grid.add(item=grank_plot, col_span=1)

      g_table = w_table(source=df_g, key="gene_results_table")

  else:
    w_text_output(
      content="No differential gene analysis found; please check Execution logs.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()
