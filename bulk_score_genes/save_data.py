new_data_signal()

if not adata_g:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

gene_score_done_signal()
if gene_score_done_signal.sample() is True:

  w_text_output(content="""
  Click **Save H5AD Data** to save your custom H5 Viewer annotations; new annotations with be available in your next session. <br>  Click **Copy Annotations to motif data** to add new annotations to motif AnnData object.
  """)

  save_button_gs = w_button(label="Save H5AD Data")
  save_warning_gs = w_text_output(content="""_This operation may take a couple minutes._""")

  save_col_gs = w_column(items=[save_button_gs, save_warning_gs])

  if save_button_gs.value:
      save_path_gs = adata_g_path

      w_text_output(
        content="Writing data to disk...",
        key="writing_message",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      try:
        adata_g.write(save_path_gs.name())
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
        save_path_gs.upload_from(Path(save_path_gs.name()))
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
