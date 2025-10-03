new_data_signal()

if not adata:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

h5_viewer_signal()
if h5_viewer_signal.sample() is True:

  w_text_output(content="""Click below to save your custom H5 Viewer annotations.
  
  Data will be stored in Latch Data and can be loaded back into your next Plots session.""")
  

  save_button = w_button(label="Save H5 Data")

  w_text_output(content="""_This operation may take a couple minutes._""")

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