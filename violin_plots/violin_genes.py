w_text_output(content="""

# Violin Plot (Gene/Cell Data)

<details>
<summary><i>details</i></summary>
For a provided grouping (Clusters, Samples, Conditions), display a violin plot showing the distribution for numeric cell data or gene accessibility data.
</details>

""")

if not adata_g:
    w_text_output(
        content="No data gene activity data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_g, anndata.AnnData):
    w_text_output(
       content="No gene activity data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

numeric_metadata = [data for data in available_metadata if data not in groups]

violin_data = w_select(
  label="data",
  default="tsse",
  options=tuple(numeric_metadata + available_genes),
  appearance={
    "help_text": "Select gene activity or metadata/QC data to plot."
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

data_type = "obs" if violin_data.value in numeric_metadata else "gene"

w_row(items=[violin_data, violin_group_by])

gv_button = w_button(label="Update Violin Plot")

if gv_button.value:
  
  violin_df = create_violin_data(
    adata_g, violin_group_by.value, violin_data.value, data_type=data_type
  )
  print(f"This cell has run. {violin_group_by.value, violin_data.value}")
  
  fig_11419 = px.violin(
      violin_df,
      x='group',
      y='value',
      box=True,
      points=False,
      color='group',
      color_discrete_sequence=px.colors.qualitative.Alphabet,
  )
  fig_11419.update_layout(
      title=f"Distribution of {violin_data.value} by {violin_group_by.value}",
      xaxis_title=violin_group_by.value,
      yaxis_title=violin_data.value,
      plot_bgcolor='rgba(0,0,0,0)',
      showlegend=False
  )
  fig_11419.update_xaxes(showgrid=False)
  fig_11419.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')
  
  w_plot(source=fig_11419)
