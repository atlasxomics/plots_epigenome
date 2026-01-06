w_text_output(content="## H5 Viewer")

w_text_output(content="""
<details>
<summary><i>H5 Viewer Instructions</i></summary>

<details>
<summary><i>H5 Viewer Instructions</i></summary>

Once initialized, the viewer will display your project cells plotted in UMAP space.

## Navigating the Viewer
- **Change coordinates:**
  - `X_umap`: UMAP coordinates
  - `X_dataset`: samples arranged in spatial layout
  - `spatial`: *not recommended* (overlays all samples at once)
- **Coloring cells:**
  - In the left panel, click the **paint bucket** icon next to any annotation to color by that variable.
  - In the *Genes of Interest* section (bottom of left panel), click **+**, type a gene name, select it, and click **Add**. Then click the paint bucket icon next to the feature to color cells by that gene or motif.

---

## Selecting Cells for the Compare Clusters Workflow

The H5 Viewer can be used to create custom annotations in the AnnData object. These annotations can then be used as input to the **Compare Clusters Workflow**.

- **Add a new annotation:**  
  In the left panel, click the **+** next to either *Continuous* or *Categorical* observations.
- **Select cells with the lasso tool:**
  - Hover over the AnnData figure, click the **lasso icon** in the Plotly toolbar (top-right).
  - Use your cursor to manually outline cells of interest.
  - To select multiple regions, **hold Shift**.
- Once cells are lasso-selected, a new **+** will appear next to the observation fields. Click it to add your selected cells to the annotation.

---

## Using Filters to Define New Observations

- **Filter cells:**  
  Click the **ellipsis (…)** next to an observation → select **Filter**, then enter your criteria.  
- **Capture filtered cells:**  
  Use the lasso tool to select the filtered cells and add them to a custom annotation as described above.  
- **Clear filters:**  
  Click the **filter icon** in the bottom-right of the display panel to remove active filters.  
- Repeat this process to create additional filter-defined groups in your AnnData object.

</details>

</details>
""")

new_data_signal()
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    exit()

if "choose_h5_data" not in globals():
  adata_h5 = adata_g

if "h5_button" in globals():
  if h5_button.value:
  
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
