new_data_signal()

w_text_output(content="""

  # H5 Viewer

  This tab contains an interactive plotting module for AnnData objects: the H5 Viewer.  To initalize the viewer,
  select the features you'd like to display from the dropdown menu below and click the 'Start H5 Viewer button'.
  
  <details>
  <summary><i>H5 Viewer Instructions</i></summary>
  Once the Start H5 Viewer button is clicked, the window initialize below display the project cells plotted in UMAP space.
  
  - To change the coordinates the cells are plotted on, choose from the dropdown menu in the upper left-hand corner of the main pannel. 
  
    - The dropdown wil intially be set to 'X_umap', which displays the cells according to UMAP coordinates.
    - To display all samples spatially, select '**X_dataset**'.
    - Do not use 'spatial'; this will overlay the samples on the same coordinates.
  - The left-hand pannel contains cell annotation that the plot can be colored by.  To color cells, click the paint brush icon next to the desired annotation.
  - The 'Genes of Interest' section on the right-hand pannel enable cells to be colored by feaeture scores.
   - To add a feature, click the + icon and start typing the gene name.  Select the feature and click 'Add' to add it as a cell annotation.
   - To color the cells by the feature value, click the paintbrush icon. 
   - To color cells by the average of multiple features, **hold `shift` and click the paintbrush icon for each feature**.
  </details>

""")

if not adata:
  w_text_output(
    content="No data data loaded...",
    appearance={"message_box": "warning"}
  )
  exit()

choose_h5_data = w_select(
  label="Select Data for H5 Viewer",
  default=None,
  options=["gene", "motif"],
  appearance={
    "help_text": "Select which features to display in the h5 Viewer; either selection displays metadata."
  }
)



h5_button = w_button(label="Start H5 Viewer")
