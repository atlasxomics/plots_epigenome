new_data_signal()

if not adata:
  exit()

confirm_cells = w_button(label="Confirm Cell Selections")

if confirm_cells.value:
  if len(groupA_cells) > 0 and len(groupB_cells) > 0:
    shared = set(groupA_cells) & set(groupB_cells)
    if shared:
      w_text_output(
        content=f"Selections share {len(shared)} cells; please confirm selections do not overlap.",
        appearance={"message_box": "danger"}
      )
      submit_widget_state()
      barcodes_signal(False)
    else:
      barcodes_signal(True)
      w_text_output(
        content=f"Selections are ready for the Workflow!.",
        appearance={"message_box": "success"}
      )
      submit_widget_state()
  else:
      w_text_output(
        content=f"Group A contains {len(groupA_cells)} cells, Group B contains {len(groupB_cells)}; please ensure selections for both groups are greater than 0.",
        appearance={"message_box": "danger"}
      )
      submit_widget_state()
      barcodes_signal(False)
