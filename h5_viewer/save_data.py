new_data_signal()

if not adata:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

h5_viewer_signal()
if h5_viewer_signal.sample() is True:

  w_text_output(content="""
  Click **Save H5 Data** to save your custom H5 Viewer annotations; new annotations with be available in your next session. <br>  Click **Synch H5 Data** to ensure both motif and gene objects have the same categorical annotations.
  """)
  

  save_button = w_button(label="Save H5 Data")
  save_warning = w_text_output(content="""_This operation may take a couple minutes._""")

  save_col = w_column(items=[save_button, save_warning])

  # Dropdown for selecting source object for synchronization
  sync_source_dropdown = w_select(
    label="Sync Source (object to copy from)",
    options=["", "Gene object", "Motif object"],
    default="",
    appearance={"placeholder": "Select source object..."}
  )

  synch_button = w_button(label="Synch H5 Data")

  synch_col = w_column(items=[sync_source_dropdown, synch_button])

  with w_grid(columns=4) as grid_save:
    grid_save.add(item=save_col, col_span=3)
    grid_save.add(item=synch_col, col_span=1)


  if save_button.value:
      if choose_h5_data.value == "gene":
        save_path = adata_g_path
      elif choose_h5_data.value == "motif":
        save_path = adata_m_path

      w_text_output(
        content="Writing data to disk...",
        key="writing_message",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      try:
        adata_h5.write(save_path.name())
      except:
        w_text_output(
          content="Write to disk failed...",
          key="writing_failed",
          appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

      upload_message = w_text_output(
        content="Uploading data to Latch Data...",
        key="upload_message",
        appearance={"message_box": "info"}
      )
      try:
        save_path.upload_from(Path(save_path.name()))
      except:
        w_text_output(
          content="Upload failed...",
          key="upload_failed",
          appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

      w_text_output(
        content="Upload success!",
        key="upload_success",
        appearance={"message_box": "success"}
      )
      submit_widget_state()

  if synch_button.value:
    # Check if a source object is selected
    if not sync_source_dropdown.value or sync_source_dropdown.value == "":
      w_text_output(
        content="Please select a source object from the dropdown before synchronizing.",
        key="sync_warning",
        appearance={"message_box": "warning"}
      )
    else:
      try:
        # Determine which object to use as source based on dropdown selection
        if sync_source_dropdown.value == "Gene object":
          sync_obs_metadata(adata_g, adata_m)
          w_text_output(
            content="Synch success! Gene metadata copied to Motif. Please restart the H5 Viewer.",
            key="synch_success",
            appearance={"message_box": "success"}
          )
        elif sync_source_dropdown.value == "Motif object":
          sync_obs_metadata(adata_m, adata_g)
          w_text_output(
            content="Synch success! Motif metadata copied to Gene. Please restart the H5 Viewer.",
            key="synch_success",
            appearance={"message_box": "success"}
          )
      except ValueError as e:
        w_text_output(
          content=f"Failed to synch with exception {e}",
          appearance={"message_box": "warning"}
        )

