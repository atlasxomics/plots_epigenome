w_text_output(content="""

# Violin Plot

<details>
<summary><i>details</i></summary>
For a provided grouping (Clusters, Samples, Conditions), display a violin plot showing the distribution for numeric cell data, gene accessibility, or motif deviation data.
</details>

""")

new_data_signal()

if not adata:
    w_text_output(
        content="No data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata, anndata.AnnData):
    w_text_output(
       content="No AnnData loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

numeric_metadata = [data for data in available_metadata if data not in groups]

violin_data = w_select(
  label="data",
  default="tsse",
  options=tuple(numeric_metadata + available_motifs + available_genes),
  appearance={
    "help_text": "Select values to plot."
  }
)
violin_group_by = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select group to display on x-axis."
  }
)

data_type = None
if violin_data.value in numeric_metadata:
  data_type = "obs"
elif violin_data.value in available_genes:
  data_type = "gene"
elif  violin_data.value in available_motifs:
  data_type = "motif"

if not data_type:
  w_text_output(
    content="Selected data not found in AnnData object...",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()


w_row(items=[violin_data, violin_group_by])

v_button = w_button(label="Update Violin Plot")

if v_button.value and data_type is not None:

  v_adata = adata_m if data_type == "motif" else adata_g
  
  violin_df = create_violin_data(
    v_adata, violin_group_by.value, violin_data.value, data_type=data_type
  )
  
  violin_fig = px.violin(
      violin_df,
      x='group',
      y='value',
      box=True,
      points=False,
      color='group',
      color_discrete_sequence=px.colors.qualitative.Alphabet,
  )
  violin_fig.update_layout(
      title=f"Distribution of {violin_data.value} by {violin_group_by.value}",
      xaxis_title=violin_group_by.value,
      yaxis_title=violin_data.value,
      plot_bgcolor='rgba(0,0,0,0)',
      showlegend=False
  )
  violin_fig.update_xaxes(showgrid=False)
  violin_fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

  w_plot(source=violin_fig)

