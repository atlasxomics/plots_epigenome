new_data_signal()

w_text_output(content="""
# H5 Viewer: Interactive Plotting for AnnData Objects

This tab contains the **H5 Viewer**, an interactive plotting module for AnnData objects.  
To initialize the viewer, select the features (genes or motifs) you’d like to display from the drop-down menu below, then click **Start H5 Viewer**.


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
""")

# Abort if no data loaded
if not adata:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    exit()

# Choose whether to display gene or motif data
choose_h5_data = w_select(
    label="Select Data for H5 Viewer",
    default=None,
    options=h5data_dict.keys(),
    appearance={
        "help_text": "Select which features to display in the H5 Viewer."
    }
)

# Checkbox to show/hide layout controls
sample_layout_button = w_checkbox(
    label="Change H5 Viewer spatial arrangement",
    default=False,
    appearance={"description": "Toggle to specify rows, columns, and spacing."}
)

if sample_layout_button.value:

    h5_cols = w_text_input(
        label="Number of Columns",
        key="cols",
        default=None,
        appearance={"help_text": "Specify the number of columns in the layout."}
    )
    h5_rows = w_text_input(
        label="Number of Rows",
        key="rows",
        default=None,
        appearance={"help_text": "Specify the number of rows in the layout."}
    )
    h5_spacing = w_text_input(
        label="Spacing Between Samples",
        key="x",
        default="100.0",
        appearance={"help_text": "Specify the spacing between samples."}
    )

    h5_groupby = w_checkbox(
        label="Group spatial by condition",
        key="h5_groupby",
        default=False,
        appearance={"help_text": "Specify the spacing between samples."}
    )

    if sample_layout_button.value:
      h5_condition_order = w_text_input(
          label="Order of conditions",
          key="h5_condition_order",
          default=None,
          appearance={"help_text": "Enter comma-seperated list of conditions to specify order."}
      )

    layout_col1 = w_column(items=[h5_cols, h5_rows, h5_spacing])
    layout_col2 = w_column(items=[h5_groupby, h5_condition_order])

    layout_row = w_row(items=[layout_col1, layout_col2])

# Button to start the H5 viewer
h5_button = w_button(label="Start H5 Viewer")

reset_tab = w_button(label="Reset Tab")

if reset_tab.value:
    # Reset core signals
    choose_group_signal(False)
    groupselect_signal(False)
    barcodes_signal(False)
    wf_exe_signal(False)
    wf_results_signal(False)
    wf_bigwigs_signal(False)

    # Reset metadata widgets to their defaults
    choose_h5_data._signal(None)
    sample_layout_button._signal(False)

    if "adata_h5" in globals():
      del adata_h5

    if "h5_cols" in globals():
      h5_cols._signal(None)
    if "h5_rows" in globals():
      h5_rows._signal(None)
    if "h5_spacing" in globals():
      h5_spacing._signal("100.0")    

    if "choose_obs" in globals():
      choose_obs._signal(None)

    if "groupA_ann" in globals():
      groupA_ann._signal(None)
    if "groupB_ann" in globals():
      groupB_ann._signal(None)

    if "wf_name" in globals():
      wf_name._signal("")
    if "wf_genome" in globals():
      wf_genome._signal(None)

    # Reset the update plot button
    if "h5_button" in globals():
      h5_button._signal(False)

    if "groupA_cells" in globals():
      groupA_cells = []
    if "groupB_cells" in globals():
      groupB_cells = []

    # Reset reset tab to avoid loop
    reset_tab._signal(False)

    # Ensure all cells initalize
    new_data_signal(True)

    if "load_compare_box" in globals():
      load_compare_box._signal(False)
    
    if "load_compare_button" in globals():
      load_compare_button._signal(False)

    if "compare_genome" in globals():
      compare_genome._signal(None)
      
    if "compare_path" in globals():
      compare_path._signal(None)

# Use a two-column grid: left placeholder, right button
with w_grid(columns=4) as grid_top:
    grid_top.add(item=h5_button, col_span=3)
    grid_top.add(item=reset_tab, col_span=1)