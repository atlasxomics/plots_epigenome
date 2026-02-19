new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

heatmap_signal()
if heatmap_signal.sample() is True and choose_heatmap_data.value is not None:

  adata_hm = h5data_dict[choose_heatmap_data.value]
  hm_feats = choose_heatmap_data.value

  hm_group_widget = w_select(
    label="group",
    default="cluster",
    options=tuple(groups),
    appearance={
      "detail": "(cluster, sample, condition)",
      "help_text": "Select grouping for heatmap."
    }
  )
  hm_group = hm_group_widget.value

  try:
    feature_label, stats_key, stats_df, _ = resolve_heatmap_stats_table(
      adata_hm, hm_feats, hm_group
    )
    group_col, feature_col, sig_col = get_heatmap_stats_columns(
      stats_df, stats_key
    )
  except ValueError as e:
    w_text_output(
      content=str(e),
      appearance={"message_box": "warning"}
    )
    exit()

  hm_sig_threshold = w_text_input(
    label="Significance threshold (FDR / adjusted p-value)",
    default="0.01",
    appearance={"help_text": "Marker significance cutoff."}
  )

  hm_effect_threshold = w_text_input(
    label="Log2FoldChange threshold",
    default="0.5" if hm_feats == "gene" else "0.1",
    appearance={
      "help_text": "Directional Log2FoldChange cutoff for selecting markers."
    }
  )

  hm_effect_direction = w_select(
    label="Log2FoldChange direction",
    default="positive",
    options=("positive", "negative", "absolute"),
    appearance={
      "help_text": "Direction used when selecting top markers per group; use 'positive' for upregulated, 'negative' for downregulated."
    }
  )

  hm_z_clip = w_text_input(
    label="Gene z-score clip (max/min)",
    default="2.0",
    appearance={
      "help_text": "Max/min of color scale, only applies to gene heatmaps."
    }
  )

  hm_top_n = w_text_input(
    label="Top features per group",
    default="25",
    appearance={"help_text": "Used only when feature list is empty."}
  )

  hm_feature_list = w_text_input(
    label="Feature list (optional, comma-separated)",
    default="",
    appearance={"help_text": "If provided, this list overrides Top features per group."}
  )

  controls_row2 = w_row(items=[hm_sig_threshold, hm_effect_threshold, hm_effect_direction])
  controls_row3 = w_row(items=[hm_z_clip, hm_top_n, hm_feature_list])

  value_metric = "Log2FC"
  rank_metric = "Log2FC"
  effect_direction = hm_effect_direction.value

  sig_threshold = safe_float(
    hm_sig_threshold.value,
    warn_msg="Significance threshold must be numeric; defaulting to 0.01."
  )
  if sig_threshold is None:
    sig_threshold = 0.01

  effect_threshold = safe_float(
    hm_effect_threshold.value,
    warn_msg="Effect-size threshold must be numeric; using default."
  )
  if effect_threshold is None:
    effect_threshold = 0.5 if hm_feats == "gene" else 0.1

  z_clip = safe_float(
    hm_z_clip.value,
    warn_msg="Row z-score clip must be numeric; defaulting to 2.0."
  )
  if z_clip is None or z_clip <= 0:
    z_clip = 2.0

  try:
    top_n = int(hm_top_n.value)
    if top_n < 1:
      raise ValueError
  except (TypeError, ValueError):
    w_text_output(
      content="Top features per group must be a positive integer; defaulting to 25.",
      appearance={"message_box": "warning"}
    )
    top_n = 25

  feature_input = hm_feature_list.value if hm_feature_list.value is not None else ""
  try:
    work_df = prepare_heatmap_work_df(
      stats_df=stats_df,
      group_col=group_col,
      feature_col=feature_col,
      value_metric=value_metric,
      rank_metric=rank_metric,
      sig_col=sig_col,
      sig_threshold=sig_threshold,
    )
    work_df, selected = select_archr_like_heatmap_features(
      work_df=work_df,
      group_col=group_col,
      feature_col=feature_col,
      rank_metric=rank_metric,
      top_n=top_n,
      effect_threshold=effect_threshold,
      effect_direction=effect_direction,
      feature_input=feature_input,
    )
    heatmap_df, legend_title = build_archr_like_heatmap_df(
      work_df=work_df,
      hm_feats=hm_feats,
      group_col=group_col,
      feature_col=feature_col,
      value_metric=value_metric,
      sig_col=sig_col,
      z_clip=z_clip,
    )
  except ValueError as e:
    w_text_output(
      content=str(e),
      appearance={"message_box": "warning"}
    )
    exit()

  if len(selected) > 0:
    ordered_rows = [f for f in selected if f in heatmap_df.index]
    ordered_rows += [f for f in heatmap_df.index if f not in ordered_rows]
    heatmap_df = heatmap_df.reindex(ordered_rows)

  title = f"{feature_label.capitalize()} heatmap by {hm_group}"
  if sig_col is not None:
    title += f" | {sig_col}<={sig_threshold}"
  if effect_direction == "positive":
    title += f" | {rank_metric}>={effect_threshold}"
  elif effect_direction == "negative":
    title += f" | {rank_metric}<=-{effect_threshold}"
  else:
    title += f" | |{rank_metric}|>={effect_threshold}"

  heatmap = px.imshow(
    heatmap_df,
    color_continuous_scale="Spectral_r",
    aspect="auto",
    origin="lower"
  )

  heatmap.update_layout(
    title=title,
    xaxis_title=hm_group,
    yaxis_title=f"{feature_label}s",
    coloraxis_colorbar=dict(
      title=legend_title,
      title_side="right"
    )
  )

  heatmap.update_xaxes(
    side="bottom",
    tickmode="array",
    tickvals=list(range(len(heatmap_df.columns))),
    ticktext=heatmap_df.columns,
    tickangle=45
  )
  heatmap.update_yaxes(autorange="reversed")

  w_plot(source=heatmap)
