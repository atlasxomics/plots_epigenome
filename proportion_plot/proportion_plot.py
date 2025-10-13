new_data_signal()

w_text_output(content="""

# Proportion Plot

Generate stacked bar plots where the **x-axis** represents a primary grouping (e.g., conditions) and the **stacked segments** represent sub-groups (e.g., clusters).  
The **y-axis** can be displayed as either raw counts or proportions of cells.

""")

# Abort if no data loaded
if not adata:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

prop_groups = [
  key for key in adata.obs_keys() if
  pd.api.types.is_object_dtype(adata.obs[key]) or pd.api.types.is_categorical_dtype(adata.obs[key])
]

group_by = w_select(
    label="group by",
    default="sample",
    options=tuple(prop_groups),
    appearance={
        "detail": "(cluster, sample, condition)",
        "help_text": "Select group to display on x-axis."
    }
)

stack_by = w_select(
    label="stack by",
    default="cluster",
    options=tuple(prop_groups),
    appearance={
        "detail": "(cluster, sample, condition)",
        "help_text": "Select group to stack the bars by."
    }
)
return_type = w_select(
    label="return type",
    default="proportion",
    options=("proportion", "counts"),
    appearance={
        "help_text": "Display cell counts or proportions per grouping."
    }
)

prop_row = w_row(items=[group_by, stack_by, return_type])

stacked_df = create_proportion_dataframe(
    adata, group_by.value, stack_by.value, return_type=return_type.value
)

# Create proportion plot from the stacked dataframe created above
proportion_plot = px.bar(
    stacked_df,
    x="group_by",
    y="value",
    color="stack_by",
    barmode="stack",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    title=f"Distribution of {stack_by.value} by {group_by.value}"
)

# Update layout
proportion_plot.update_layout(
    xaxis_title=group_by.value,
    yaxis_title="Proportion" if return_type.value == "proportion" else "Count",
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=True,
    legend_title=stack_by.value
)

# Update axes
proportion_plot.update_xaxes(showgrid=False)
proportion_plot.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

# Determine ordering: numeric sort if possible, else alphabetical
p_cats = stacked_df['group_by'].unique().tolist()
p_category_order = sort_group_categories(p_cats)

proportion_plot.update_xaxes(categoryorder='array', categoryarray=p_category_order)

w_plot(source=proportion_plot)
