w_text_output(content="## H5 Viewer Settings")

new_data_signal()

# Abort if no data loaded
if not adata_g:
    w_text_output(
        content=" ",
    )
    exit()

# Choose whether to display gene or motif data
choose_h5_data = w_select(
    label="Select Data for H5 Viewer",
    default="gene",
    options=["gene", "motif"],
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

    h5_sortby_opts = ["original", "sample"]
    if "condition" in groups:
      h5_sortby_opts.append("condition")

    h5_cols = w_text_input(
        label="Number of Columns",
        key="h5_cols",
        default=None,
        appearance={"help_text": "Specify the number of columns in the layout."}
    )
    h5_rows = w_text_input(
        label="Number of Rows",
        key="h5_rows",
        default=None,
        appearance={"help_text": "Specify the number of rows in the layout."}
    )
    h5_spacing = w_text_input(
        label="Spacing Between Samples",
        key="h5_spacing",
        default="100.0",
        appearance={"help_text": "Specify the spacing between samples."}
    )

    h5_flipy = w_checkbox(
        label="Flip Y Axis",
        key="h5_flipy",
        default=False,
        appearance={"description": "Rotate samples around the y axis."}
    )

    h5_sortby = w_select(
        label="Sort Samples By",
        key="h5_sortby",
        default="original",
        options=tuple(h5_sortby_opts),
        appearance={
          "help_text": "Sort samples alphabetically or by condition.",
          "description": "'original' maintains original order."
        }
    )

    layout_col1 = w_column(items=[h5_cols, h5_rows])
    layout_col2 = w_column(items=[h5_spacing, h5_sortby, h5_flipy])

    with w_grid(columns=5) as layout_grid:
        layout_grid.add(item=layout_col1, col_span=1)
        layout_grid.add(item=layout_col2, col_span=4)
  
# Button to start the H5 viewer
h5_button = w_button(label="Refresh H5 Viewer")

