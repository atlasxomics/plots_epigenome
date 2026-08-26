new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

wf_results_signal()
if wf_results_signal.sample() == True:
  
  if "motif" in results_dict.keys():
  
    w_text_output(content="""## Differential Motif Enrichment""")
    
    m_group = w_select(
        label="Group",
        default=results_dict["motif"]["group_name"].unique()[0],
        options=tuple(results_dict["motif"]["group_name"].unique()),
        key="m_group"
    )
    
    m_pvals_adj_threshold = w_text_input(
      label="pval adjust threshold",
      default="0.05",
      key="m_pvals_adj_threshold"
    )
    
    m_meandiff_threshold = w_text_input(
      label="MeanDiff threshold",
      default="0.01",
      key="m_meandiff_threshold"
    )
    
    mcompare_rankby = w_select(
        label="Rank By",
        default="MeanDiff",
        options=tuple(['FDR', 'MeanDiff']),
        key="mcompare_rankby"
    )
    
    mcompare_colorby = w_select(
        label="Color By",
        default="FDR",
        options=tuple(['FDR', 'MeanDiff']),
        key="mcompare_colorby"
    )

    
    with w_grid(
      columns=4,
      key="motif_results_controls_grid"
    ) as motif_controls_grid:
      motif_controls_grid.add(item=m_pvals_adj_threshold, col_span=1)
      motif_controls_grid.add(item=m_meandiff_threshold, col_span=1)
      motif_controls_grid.add(item=mcompare_rankby, col_span=1)
      motif_controls_grid.add(item=mcompare_colorby, col_span=1)
    
    # ----------------------------------------------------------------------
    
    if m_group.value is not None:
      
      
      group_m = m_group.value
      df_m = results_dict["motif"]
      df_m = df_m[df_m["group_name"] == group_m]
      
      mvol = plot_volcano(
        df_m,
        float(m_pvals_adj_threshold.value),
        float(m_meandiff_threshold.value),
        "GroupA",
        "GroupB",
        pval_key="FDR",
        l2fc_key="MeanDiff",
        names_key="name",
        plot_width=750,
        plot_height=640,
        top_n=2
      )
      
      mrank = plot_ranked_feature_plotly(
          df_m,
          y_col=mcompare_rankby.value,
          x_col=None,
          n_labels=4,
          label_col="name",
          color_col=mcompare_colorby.value,
          colorscale="PuBu_r",
          marker_size=6,
          title="",
          y_label=mcompare_rankby.value
      )

      # Align the two plots on a shared x-axis baseline: give both the same
      # plotting-area height and identical top/bottom margins so, side by side
      # in the row, their x-axes sit at the same vertical position. autosize +
      # width=None lets the row size their widths equally.
      mvol.update_layout(
          height=640,
          autosize=True,
          width=None,
          margin=dict(l=80, r=80, t=60, b=80),
      )
      mrank.update_layout(
          height=640,
          autosize=True,
          width=None,
          margin=dict(l=80, r=80, t=60, b=80),
      )

      # Remove the alias left by older versions of this reactive cell. Plot
      # widgets require each figure to have exactly one global variable name.
      globals().pop("_fig", None)

      mvol_plot = w_plot(source=mvol, key="mvol_plot")
      mrank_plot = w_plot(source=mrank, key='mrank_plot')

      with w_grid(
        columns=2,
        key="motif_results_plot_grid"
      ) as motif_results_grid:
        motif_results_grid.add(item=mvol_plot, col_span=1)
        motif_results_grid.add(item=mrank_plot, col_span=1)

      m_table = w_table(source=df_m, key="motif_results_table")
  
  else:
    w_text_output(
      content="No differential motif analysis found; please check Execution logs.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()

else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()
