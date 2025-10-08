new_data_signal()
if not adata:
  w_text_output(content="   ")
  exit()

if "h5_button" in globals():  # For some reason, this cell runs before h5_button defined
  if h5_button.value:
  
      if choose_h5_data.value is not None:
        adata_h5 = h5data_dict[choose_h5_data.value]
  
        if sample_layout_button.value:
        
            proceed = True
        
            # Ensure fields are present and non-empty after stripping
            cols_val = (h5_cols.value or "").strip() if h5_cols is not None else ""
            rows_val = (h5_rows.value or "").strip() if h5_rows is not None else ""
            spacing_val = (h5_spacing.value or "").strip() if h5_spacing is not None else ""

            if h5_flipy is not None and isinstance(h5_flipy.value, bool):
                flipy_val = h5_flipy.value
            else:
                flipy_val = False
            
            valid_sort_modes = {"original", "sample", "condition"}
            if h5_sortby is not None and (h5_sortby.value in valid_sort_modes):
                sort_val = h5_sortby.value
            else:
                sort_val = "original"
        
            if cols_val and rows_val and spacing_val:
        
                try:
                    n_cols = int(cols_val)
                except (TypeError, ValueError):
                    proceed = False
                    w_text_output(
                        content="Cannot convert 'Number of Columns' input into an integer; ignoring...",
                        appearance={"message_box": "warning"}
                    )
        
                try:
                    n_rows = int(rows_val)
                except (TypeError, ValueError):
                    proceed = False
                    w_text_output(
                        content="Cannot convert 'Number of Rows' input into an integer; ignoring...",
                        appearance={"message_box": "warning"}
                    )
        
                try:
                    spacing = float(spacing_val)
                except (TypeError, ValueError):
                    proceed = False
                    w_text_output(
                        content="Cannot convert 'Spacing Between Samples' input into a float; ignoring...",
                        appearance={"message_box": "warning"}
                    )
        
                # Optional: basic range checks
                if proceed:
                    if n_cols < 1 or n_rows < 1 or spacing <= 0:
                        proceed = False
                        w_text_output(
                            content="Rows/Columns must be ≥ 1 and Spacing must be > 0; ignoring...",
                            appearance={"message_box": "warning"}
                        )
  
                if proceed:
                  total_positions = n_rows * n_cols
                  if len(samples) > total_positions:
                    proceed = False
                    w_text_output(
                      content=f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {len(samples)} samples",
                      appearance={"message_box": "warning"}
                    )
        
                if proceed:
                    new_obsm = f"spatial_offset_{n_rows}x{n_cols}-{spacing}-{('FlipY' if flipy_val else 'noFlipY')}-{sort_val}"
                    process_matrix_layout(
                        adata_h5,
                        n_rows=n_rows,
                        n_cols=n_cols,
                        tile_spacing=spacing,
                        flipy=flipy_val,
                        sample_order_mode=sort_val,
                        new_obsm_key=new_obsm
                    )
        
            else:
                w_text_output(
                    content="Please complete all fields in 'Change H5 Viewer spatial arrangement' to specify layout.",
                    appearance={"message_box": "warning"}
                )
              
        
        viewer = w_h5(ann_data=adata_h5)
        h5_viewer_signal(True)
  
      else:
        w_text_output(
          content="Please specify features (genes, motifs) to display from drop-down above",
          appearance={"message_box": "warning"}
        )

else:
  w_text_output(content="   ")
  exit()
