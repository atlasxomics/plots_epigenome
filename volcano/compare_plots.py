new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

compare_signal()
if compare_signal.sample() == True and choose_compare_data.value is not None:

  adata_comp = h5data_dict[choose_compare_data.value]
  
  c_condition = w_select(
      label="Condition",
      default=conditions[0],
      options=tuple(conditions),
      appearance={
        "help_text": "Select condition for comparison."
      },
      key="c_condition"
  )

  c_cluster = w_select(
      label="Cluster",
      default="All",
      options=tuple(clusters + ["All"]),
      appearance={
          "help_text": "Filter data to a specific cluster."
      },
      key="c_cluster"
  )
  c_all_opts = w_row(
    items=[c_condition, c_cluster],
    key="c_all_opts"
  )

  c_pvals_adj_threshold = w_text_input(
    label="pval adjust Threshold",
    default="0.05",
    key="c_pvals_adj_threshold"
  )
  
  c_log2fc_threshold = w_text_input(
    label="Difference Metric Threshold",
    default="0.01",
    appearance={"help_text": "Log2FC for genes, MeanDiff for motifs."},
    key="c_log2fc_threshold"
  )
  

  c_rankby = w_select(
      label="Rank By",
      default=rankby_default,
      options=tuple(rankby_opts),
      key="c_rankby"
  )
  
  c_colorby = w_select(
      label="Color By",
      default="p_val_adj",
      options=tuple(rankby_opts),
      key="c_colorby"
  )
  
  with w_grid(
    columns=4,
    key="compare_plots_controls_grid"
  ) as compare_controls_grid:
    compare_controls_grid.add(item=c_pvals_adj_threshold, col_span=1)
    compare_controls_grid.add(item=c_log2fc_threshold, col_span=1)
    compare_controls_grid.add(item=c_rankby, col_span=1)
    compare_controls_grid.add(item=c_colorby, col_span=1)
  
  # ----------------------------------------------------------------------
  
  if c_condition.value is not None:
    
    c_df = adata_comp.uns[f"volcano_1_{c_condition.value}"]
    c_df = c_df[c_df["cluster"] == c_cluster.value]
    if len(c_df) == 0:
      w_text_output(
         content=f"There is no volcano plot for cluster {c_cluster.value} because it contains more than 90% of one of the conditions. Please check Proportion plot.",
         appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit(0)

    try:
      c_df.drop(["Significance"], axis=1, inplace=True)
    except:
      print("No Significance column found")

    c_rank_metric = c_rankby.value
    if c_rank_metric in ["p_val", "p_val_adj"]:
      c_df[f"-log10{c_rank_metric}"] = -np.log10(c_df[c_rank_metric])
      c_rank_metric = f"-log10{c_rank_metric}"
    
    c_vol = plot_volcano(
      c_df,
      float(c_pvals_adj_threshold.value),
      float(c_log2fc_threshold.value),
      c_condition,
      "rest",
      pval_key="p_val",
      l2fc_key=diff_metric,
      names_key="gene",
      title="",
      plot_width=750,
      plot_height=640,
      top_n=2
    )
    c_rank = plot_ranked_feature_plotly(
        c_df,
        y_col=c_rank_metric,
        x_col=None,
        n_labels=4,
        label_col="gene",
        color_col=c_colorby.value,
        colorscale="PuBu_r",
        marker_size=6,
        title=f"Condition {c_condition.value} versus rest (cluster {c_cluster.value})",
        y_label=c_rank_metric
    )

    c_vol.update_layout(
      height=640,
      autosize=True,
      width=None,
      margin=dict(l=80, r=80, t=60, b=80),
    )
    c_rank.update_layout(
      height=640,
      autosize=True,
      width=None,
      margin=dict(l=80, r=80, t=60, b=80),
    )

    # Remove the alias left by older versions of this reactive cell. Plot
    # widgets require each figure to have exactly one global variable name.
    globals().pop("_fig", None)

    c_vol_plot = w_plot(source=c_vol, key="c_vol_plot")
    c_rank_plot = w_plot(source=c_rank, key="c_rank_plot")

    with w_grid(
      columns=2,
      key="compare_plots_plot_grid"
    ) as compare_results_grid:
      compare_results_grid.add(item=c_vol_plot, col_span=1)
      compare_results_grid.add(item=c_rank_plot, col_span=1)

    c_table = w_table(source=c_df, key="compare_plots_table")

else:
  w_text_output(
    content="  "
  )
  submit_widget_state()
