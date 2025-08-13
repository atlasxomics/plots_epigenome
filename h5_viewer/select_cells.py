new_data_signal()

w_text_output(content="""# Compare Clusters Workflow""")

w_text_output(content="""

  ## Select Cells for Compare Workflow
  
  <details>
  <summary><i>Instructions</i></summary>
  In H5 Viewer, use the lasso select tool to select the cells for comparison Workflow.
  </details>

""")

if not adata:
  w_text_output(
    content="No data data loaded...",
    appearance={"message_box": "warning"}
  )
  exit()

select_type = w_radio_group(
    label="Choose Input Type",
    options=["Lasso-select", "Observation Annotations"],
    default="Lasso-select"
)

if select_type.value == "Lasso-select":
  viewer.sample()

  b_text = w_text_output(content="Use the lasso-select tool in the H5 Viewer to circle the cells for Group A and Group B.  Click the 'Select Cells' button to load and view the selection.") 
  
  groupA_button = w_button(label="Select Group A Cells")
  groupB_button = w_button(label="Select Group B Cells")
  
  # Create a grid with 3 columns
  with w_grid(columns=2) as g1:
      g1.add(item=groupA_button, col_span=1)
      g1.add(item=groupB_button, col_span=1)
  
  if groupA_button.value:
  
    try:
      if viewer.value["lasso_points"] is not None: 
        groupA_cells = get_barcodes_by_lasso(
          adata=adata_h5,
          obsm_key=viewer.value["lasso_points_obsm"],
          lasso_points=viewer.value["lasso_points"][0]
        )
      else:
        groupA_cells = []
  
      a_text = w_text_output(content=f"{len(groupA_cells)} cells selected for Group A.", appearance={"message_box": "success"})
      groupA_obs = adata_h5.obs.loc[groupA_cells]
      a_table = w_table(source=groupA_obs)
  
    except NameError:
      w_text_output(content="Please ensure H5 Viewer has loaded.", appearance={"message_box": "warning"})
    
  
  if groupB_button.value:
  
    try:
      if viewer.value["lasso_points"] is not None: 
        groupB_cells = get_barcodes_by_lasso(
          adata=adata_h5,
          obsm_key=viewer.value["lasso_points_obsm"],
          lasso_points=viewer.value["lasso_points"][0]
        )
      else:
        groupB_cells = []
    
      b_text = w_text_output(content=f"{len(groupB_cells)} cells selected for Group B.", appearance={"message_box": "success"}) 
      groupB_obs = adata_h5.obs.loc[groupB_cells]
      b_table = w_table(source=groupB_obs)
  
    except NameError:
      w_text_output(content="Please ensure H5 Viewer has loaded.", appearance={"message_box": "warning"})
  
  groupA_text = w_text_output(content=f"""<details><summary><i>Group A barcodes</i></summary>{','.join(list(groupA_cells))}</details>""")
  groupB_text = w_text_output(content=f"""<details><summary><i>Group B barcodes</i></summary>{','.join(list(groupB_cells))}</details>""")
  
  with w_grid(columns=2) as g2:
      # Add plots to the grid with different spans
      g2.add(item=groupA_text, col_span=1)
      g2.add(item=groupB_text, col_span=1)

elif select_type.value == "Observation Annotations":
  try:
    viewer.sample() # trigger refresh

    choose_obs = w_select(
      label="Select Annotation",
      default=None,
      options=adata_h5.obs_keys(),
      appearance={
        "help_text": "Selection with Annotation (.obs column) to define groups by."
      }
    )

    obs_type = w_radio_group(
      label="Annotation Type",
      options=["Categorical", "Continuous"],
      default="Lasso-select"
    )

    if choose_obs.value is not None:
    
      if obs_type.value == "Categorical":
        ann_keys = adata_h5.obs[choose_obs.value].unique()

        groupA_ann = w_select(
          label="Group A Value",
          default=None,
          options=ann_keys,
          appearance={
            "help_text": "Select value for Group A."
          }
        )
      
        groupB_ann = w_select(
          label="Group B Value",
          default=None,
          options=ann_keys,
          appearance={
            "help_text": "Select value for Group B."
          }
        )

        if groupA_ann.value is not None and groupB_ann.value is not None:
          if groupA_ann.value == groupB_ann.value:
            group_id_error = w_text_output(content="Please ensure Group A and Group B values are different.", appearance={"message_box": "warning"})
            submit_widget_state()

          groupA_cells = list(adata_h5[adata_h5.obs[choose_obs.value] == groupA_ann.value].obs_names)
          groupB_cells = list(adata_h5[adata_h5.obs[choose_obs.value] == groupB_ann.value].obs_names)

          if len(groupA_cells) > 0 and len(groupB_cells) > 0:
            group_success = w_text_output(content=f"Group A ({groupA_ann.value}): {len(groupA_cells)} cells; Group B ({groupB_ann.value}): {len(groupB_cells)} cells", appearance={"message_box": "success"})
    
  except NameError:
    w_text_output(content="Please ensure H5 Viewer has loaded.", appearance={"message_box": "warning"})
    submit_widget_state()   
