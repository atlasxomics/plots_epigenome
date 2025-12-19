new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

heatmap_signal()
if heatmap_signal.sample() == True and choose_heatmap_data.value is not None:

  adata_hm = h5data_dict[choose_heatmap_data.value]
  hm_feats = choose_heatmap_data.value
  
  hm_group = w_select(
    label="group",
    default="cluster",
    options=tuple(groups),
    appearance={
      "detail": "(cluster, sample, condition)",
      "help_text": "Select grouping for heatmap."
    }
  )

  hm_group = hm_group.value

  if hm_feats == "gene":
    hm_legend = "Log2FC"
    feat_key = f"{hm_feats}s"
  elif hm_feats == "motif":
    hm_legend = "row-normalized –log10(adjusted p-values)[0-100]"
    feat_key = f"{hm_feats}"

  # Handle "condition" case to support ArchR or old Snap
  if hm_group == "condition":

      hm_group = "conditions1"

      if hm_feats == "gene":
        possible_keys = [
            f"{hm_feats}s_per_conditions1_hm",
            f"{hm_feats}s_per_condition_1_hm"
        ]
      elif hm_feats == "motif":
        possible_keys = [
            f"{hm_feats}_per_conditions1_hm",
            f"{hm_feats}_per_condition_1_hm"
        ]
      for key in possible_keys:
          if key in adata_hm.uns:
              heatmap_df = adata_hm.uns[key]
              break
      else:
          w_text_output(
              content=f"No features heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
  else:
      key = f"{feat_key}_per_{hm_group}_hm"
      if key not in adata_hm.uns:
          w_text_output(
              content=f"No features heatmap found for key: {key}",
              appearance={"message_box": "warning"}
          )
          exit()
      heatmap_df = adata_hm.uns[key]

  if "Unnamed: 0" in heatmap_df.columns:
    heatmap_df = heatmap_df.set_index("Unnamed: 0")
  
  title=f"{hm_feats} Score Heatmap by {hm_group}"
  
  heatmap = px.imshow(
    heatmap_df,
    color_continuous_scale='Spectral_r',
    aspect='auto',
    origin='lower'
  )
  
  heatmap.update_layout(
      title=title,
      xaxis_title=hm_group,
      yaxis_title=f"{hm_feats}s",
      coloraxis_colorbar=dict(
        title=hm_legend,
        title_side="right")
  )

  heatmap.update_xaxes(
    side="bottom",
    tickmode='array',
    tickvals=list(range(len(heatmap_df.columns))),
    ticktext=heatmap_df.columns,
    tickangle=45
  )
  heatmap.update_yaxes(autorange="reversed")
  
  w_plot(source=heatmap)
