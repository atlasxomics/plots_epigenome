new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

w_text_output(content="""## Set Inputs for Compare Workflow""")

barcodes_signal()
if barcodes_signal.sample() == True:

  wf_name = w_text_input(
    label="Output Directory Name",
    default="",
    appearance={
      "help_text": "Name of output dirctory results save in Latch Data in the directory `compare_outs`."
    }
  )

  wf_version = w_select(
    label="Workflow Version",
    default="0.9.2-331570-e34694",
    options=(
      "0.9.2-331570-e34694",
      "0.9.2.3-0cfb87-f5d981",
    ),
    appearance={
      "help_text": "Select the compare workflow version to launch."
    }
  )

  if (wf_name.value is not None and
      len(wf_name.value) > 0 and
      archrproj_dir is not None and
      remote_bcs is not None
    ):

    params = {
        "project_name": wf_name.value,
        "groupings": LatchFile(remote_bcs.path),
        "archrproject": LatchDir(archrproj_dir.path),
    }
    
    w = w_workflow(
      wf_name="wf.__init__.compare_workflow",
      version=wf_version.value,
      params=params,
      label="Launch Workflow"
    )

    wf_exe_signal(True)
    execution = w.value
    
  else:
    w_text_output(
      content="Please set Workflow inputs.",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    wf_exe_signal(False)    
  
else:
  w_text_output(
    content="Please ensure cells are selected for Group A and Group B.",
    appearance={"message_box": "neutral"}
  )
  submit_widget_state()
  wf_exe_signal(False)    
  exit()
