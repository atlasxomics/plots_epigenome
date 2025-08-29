new_data_signal()

if not adata:
  w_text_output(content="   ")
  exit()

w_text_output(content="""## Compare Workflow Results""")

if not adata:
  w_text_output(content="   ")
  exit()

wf_exe_signal()
if wf_exe_signal.sample() == True:

  get_results = w_button(label="Fetch Workflow Results")
  
  w_text_output(content="Click to load in Exexcution results.")
  
  
  if get_results.value:
    if execution is not None:

      list_resp = post(
          url=config.api.execution.list,
          headers = {"Authorization": get_auth_header()},
          json={"ws_account_id": f"{workspace_account_id}"},
      ).json()
      target_execution = list_resp[str(execution.id)]
      status = target_execution.get('status')
    
      if status == 'SUCCEEDED':

          w_text_output(
              content=f"Compare Clusters successfully finished running; loading results...",
                appearance={"message_box": "success"}
          )
          submit_widget_state()
          res = await execution.wait()
          try:
            local = res.output["o0"].local_path
          except:
            w_text_output(
              content="Could not find local path for results.",
              appearance={"message_box": "warning"}
            )
            submit_widget_state()
            exit()
          
          for feat in feats:
            res_dir = f"{local}/{feat}_results"
            
            if os.path.exists(res_dir):
              res_path = f"{res_dir}/all_{feat}s.csv"
              if os.path.exists(res_path):
                  results_dict[feat] = pd.read_csv(res_path)
    
            else:
              w_text_output(
                content=f"Could not find results for {feat}s; please check Execution logs.",
                appearance={"message_box": "warning"}
              )
              submit_widget_state()
  
          if len(results_dict.keys()) > 0:
            w_text_output(
              content=f"Loaded {', '.join(list(results_dict.keys()))} results",
              appearance={"message_box": "success"}
            )
            submit_widget_state()
            wf_results_signal(True)
          else:
            w_text_output(
              content="Could not find output tables for Execution; please check Execution logs.",
              appearance={"message_box": "warning"}
            )
            submit_widget_state()
            wf_results_signal(False)

          # Coverages
          remote_dir = res.output["o0"]
          coverages_dir = LatchDir(f"{res.output['o0'].remote_path}/coverages")

          files = []
          try:
            for file in coverages_dir.iterdir():
                suffix = file.path.split(".")[-1]
                if suffix == "bw":
                    files.append(file)
            if len(files) > 0:
              wf_bigwigs_signal(True)
              w_text_output(
                content="Loaded coverage bigwigs for Group A and Group B...",
                appearance={"message_box": "success"}
              )
              submit_widget_state()
            else:
              w_text_output(
                content="No coverage file found...",
                appearance={"message_box": "warning"}
              )
              submit_widget_state()
          except ValueError:
              w_text_output(
                content="Could not find remote coverages folder; please check Execution logs.",
                appearance={"message_box": "warning"}
              )
              submit_widget_state()

      elif status in ["ABORTED", "FAILED"]:
        w_text_output(
          content="Workflow Execution failed or aborted, please check Execution logs.",
          appearance={"message_box": "warning"}
        )
        submit_widget_state()
        wf_results_signal(False)
        exit()
  
      elif execution.status in ["UNDEFINED", "RUNNING"]:
        w_text_output(
          content="Workflow still running, click button again once it has completed.",
          appearance={"message_box": "warning"}
        )
        submit_widget_state()
        wf_results_signal(False)
        exit()
      else:
        w_text_output(
          content="Unknown Execution status...",
          appearance={"message_box": "neutral"}
        )
        submit_widget_state()
        
    else:
      w_text_output(
        content="""Awaiting Workflow launch...""",
        appearance={"message_box": "neutral"}
      )
      submit_widget_state()
      wf_results_signal(False)

else:
    w_text_output(
      content="""Awaiting Workflow launch...""",
      appearance={"message_box": "neutral"}
    )
    submit_widget_state()
    wf_results_signal(False)