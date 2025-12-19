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

  wf_genome = w_select(
    label="Reference Geneome",
    default=None,
    options=tuple(["hg38", "mm10"]),
    appearance={
      "help_text": "Reference genome for experiment."
    }
  )

  w_row(items=[wf_name, wf_genome])

  if (wf_name.value is not None and
      len(wf_name.value) > 0 and
      wf_genome.value is not None and
      archrproj_dir is not None and
      remote_bcs is not None
    ):

    params = {
        "project_name": wf_name.value,
        "groupings": LatchFile(remote_bcs.path),
        "archrproject": LatchDir(archrproj_dir.path),
        "genome": genome_dict[wf_genome.value]
    }
    
    w = w_workflow(
      wf_name="wf.__init__.compare_workflow",
      version="0.7.1-8484d6-wip-4ae938",
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


