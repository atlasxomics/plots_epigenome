w_text_output(content="""

# Violin Plot (Motif/Cell Data)

<details>
<summary><i>details</i></summary>
For a provided grouping (Clusters, Samples, Conditions), display a violin plot showing the distribution for numeric cell data or motif deviation data.
</details>

""")

if not adata_m:
    w_text_output(
        content="No motif data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

numeric_metadata = [data for data in available_metadata if data not in groups]

m_violin_data = w_select(
  label="data",
  default="tsse",
  options=tuple(numeric_metadata + available_motifs),
  appearance={
    "help_text": "Select motif deviation or metadata/QC data to plot."
  }
)
m_violin_group_by = w_select(
  label="group",
  default="cluster",
  options=tuple(groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select group to display on x-axis."
  }
)

data_type = "obs" if m_violin_data.value in numeric_metadata else "gene"

w_row(items=[m_violin_data, m_violin_group_by])

m_violin_df = create_violin_data(
  adata_m, m_violin_group_by.value, m_violin_data.value, data_type=data_type
)
print(f"This cell has run. {m_violin_group_by.value, m_violin_data.value}")

m_violin_fig = px.violin(
    m_violin_df,
    x='group',
    y='value',
    box=True,
    points=False,
    color='group',
    color_discrete_sequence=px.colors.qualitative.Alphabet,
)
m_violin_fig.update_layout(
    title=f"Distribution of {m_violin_data.value} by {m_violin_group_by.value}",
    xaxis_title=m_violin_group_by.value,
    yaxis_title=m_violin_data.value,
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=False
)
m_violin_fig.update_xaxes(showgrid=False)
m_violin_fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

m_violin_fig.show()