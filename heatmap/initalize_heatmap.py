new_data_signal()

w_text_output(content="""

# Heatmaps

Display heatmaps of pre-computed differential statistics across project groupings (e.g., **Clusters**, **Samples**, or **Conditions**) for either genes or motifs.

The plot displays the results of [ArchR::plotMarkerHeatmap](https://www.archrproject.com/reference/plotMarkerHeatmap.html) for genes and [ArchR::plotEnrichHeatmap](https://www.archrproject.com/reference/plotEnrichHeatmap.html) for motifs.

To initialize the plot, select one or more features (genes or motifs) from the dropdown menu below.

""")

# Abort if no data loaded
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

# Choose whether to display gene or motif data
choose_heatmap_data = w_select(
    label="Select Data for Heatmap Plots",
    default="gene",
    options=["gene", "motif"],
    appearance={
        "help_text": "Select which features to display in the heatmap."
    }
)

if choose_heatmap_data.value is not None:
  heatmap_signal(True)

