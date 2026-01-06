new_data_signal()

if not adata_g:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

h5_viewer_signal()
if h5_viewer_signal.sample() is True:

  w_text_output(content="""
  Click **Save H5 Data** to save your custom H5 Viewer annotations; new annotations with be available in your next session.
  """)
  

  save_button = w_button(label="Save H5 Data")
  save_warning = w_text_output(content="""_This operation may take a couple minutes._""")

  save_col = w_column(items=[save_button, save_warning])

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
