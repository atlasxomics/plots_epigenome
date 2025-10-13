w_text_output(content="""

# Violin Plot

Generate violin plots to visualize the distribution of values across a selected grouping (e.g., **Clusters**, **Samples**, or **Conditions**).  
You can display numeric cell-level metrics, gene accessibility values, or motif deviation scores.

> If both gene and motif `AnnData` objects are available, the data dropdown will include names from both.  
> Motif names are distinguished by a numeric suffix (e.g., `CTCF-177`).

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

violin_groups = [
  key for key in adata.obs_keys() if
  pd.api.types.is_object_dtype(adata.obs[key]) or pd.api.types.is_categorical_dtype(adata.obs[key])
]

available_metadata = tuple(key for key in adata.obs_keys()
                           if key not in na_keys)
numeric_metadata = [data for data in available_metadata if data not in groups]

violin_data = w_select(
  label="data",
  default="tsse",
  options=tuple(numeric_metadata + available_motifs + available_genes),
  appearance={
    "help_text": "Select values to plot.",
    "description": "Motifs are denoted with a numberic suffix ie. CTCT-177."
  }
)

violin_group_by = w_select(
  label="group",
  default="cluster",
  options=tuple(violin_groups),
  appearance={
    "detail": "(cluster, sample, condition)",
    "help_text": "Select group to display on x-axis."
  }
)

violin_type = w_select(
  label="plot type",
  default="box",
  options=tuple(["box", "violin"]),
  appearance={
    "detail": "(box, violin)",
    "help_text": "Use box to speed up large datasets."
  }
)
violin_row = w_row(items=[violin_data, violin_group_by, violin_type])

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



if data_type is not None:

  v_adata = adata_m if data_type == "motif" else adata_g
  
  violin_df = create_violin_data(
    v_adata, violin_group_by.value, violin_data.value, data_type=data_type
  )

  if violin_type.value == "box":
    violin_fig = px.box(
        violin_df,
        x='group',
        y='value',
        points=False,
        color='group',
        color_discrete_sequence=px.colors.qualitative.Alphabet,
    )

  elif violin_type.value == "violin":
    violin_fig = px.violin(
        violin_df,
        x='group',
        y='value',
        box=True,
        points=False,
        color='group',
        color_discrete_sequence=px.colors.qualitative.Alphabet,
      )

  else:
    w_text_output(
      content="Plot type not recognized",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

  violin_fig.update_layout(
      title=f"Distribution of {violin_data.value} by {violin_group_by.value}",
      xaxis_title=violin_group_by.value,
      yaxis_title=violin_data.value,
      plot_bgcolor='rgba(0,0,0,0)',
      showlegend=False
  )
  violin_fig.update_xaxes(showgrid=False)
  violin_fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

  # Determine ordering: numeric sort if possible, else alphabetical
  v_cats = violin_df['group'].unique().tolist()
  v_category_order = sort_group_categories(v_cats)
  
  violin_fig.update_xaxes(categoryorder='array', categoryarray=v_category_order)

  w_plot(source=violin_fig)

else:
  w_output_text(content="  ")
  sumbit_widget_state()