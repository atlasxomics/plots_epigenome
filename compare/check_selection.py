new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

choose_group_signal()

committed_selection = choose_group_signal.sample()
selection_ready = (
  isinstance(committed_selection, dict)
  and committed_selection.get("annotation") is not None
  and "adata_h5" in globals()
  and adata_h5 is not None
  and committed_selection["annotation"] in adata_h5.obs.columns
  and committed_selection.get("group_a") is not None
  and committed_selection.get("group_b") is not None
)

if selection_ready:
  compare_obs_val = committed_selection["annotation"]
  groupA_val = committed_selection["group_a"]
  groupB_val = committed_selection["group_b"]

  if groupA_val == groupB_val:
    w_text_output(
        content="Please ensure different values are selected for Group A and Group B.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    barcodes_signal(False)
    exit()

  # Validate the committed selection and build its workflow configuration in
  # one cell, without an intermediate validation signal and cell rerun.
  groupA_cells = list(
    adata_h5.obs_names[adata_h5.obs[compare_obs_val] == groupA_val]
  )
  groupB_cells = list(
    adata_h5.obs_names[adata_h5.obs[compare_obs_val] == groupB_val]
  )

  if len(groupA_cells) == 0 or len(groupB_cells) == 0:
    w_text_output(
      content="Both selected groups must contain at least one cell.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    barcodes_signal(False)
    exit()

  shared = set(groupA_cells) & set(groupB_cells)
  if shared:
    w_text_output(
      content="The selected groups overlap; please choose non-overlapping groups.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    barcodes_signal(False)
    exit()

  w_text_output(
    content=(
      f"Group A ({groupA_val}): {len(groupA_cells)} cells; "
      f"Group B ({groupB_val}): {len(groupB_cells)} cells"
    ),
    appearance={"message_box": "success"}
  )

  try:
    wf_name._signal(None)
  except NameError:
    pass

  remote_bcs = None
  try:
    cells_payload = {"groupA": groupA_cells, "groupB": groupB_cells}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
      json.dump(cells_payload, f)
      local_cfg = f.name
    dataset_name = Path(data_path.value.path.rstrip("/")).name or "dataset"
    dataset_slug = re.sub(r"[^A-Za-z0-9._-]+", "_", dataset_name).strip("._-")
    if len(dataset_slug) == 0:
      dataset_slug = "dataset"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    remote_bcs = LPath(
      f"{data_path.value.path}/compare_config_{dataset_slug}_A{len(groupA_cells)}_B{len(groupB_cells)}_{timestamp}.json"
    )
    remote_bcs.upload_from(Path(local_cfg))
  except Exception:
    w_text_output(
      content="Failed to upload barcodes to remote.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    barcodes_signal(False)
    exit()

  barcodes_signal(committed_selection)
  w_text_output(
    content=(
      "Selections are ready for the Workflow! "
      f"Config saved to `{remote_bcs.path}`."
    ),
    appearance={"message_box": "success"}
  )
  submit_widget_state()

else:
  barcodes_signal(False)
  w_text_output(content="   ")
  submit_widget_state()
  exit()
