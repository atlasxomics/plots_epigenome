new_data_signal()
if not adata:
  exit()

if h5_button.value:
    if choose_h5_data.value is not None:
      adata_h5 = h5data_dict[choose_h5_data.value]
      viewer = w_h5(ann_data=adata_h5)
    else:
      w_text_output(
        content="Please specify features (genes, motifs) to display from dropdown above",
        appearance={"message_box": "warning"}
      )