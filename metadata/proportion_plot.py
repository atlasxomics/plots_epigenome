w_text_output(content="""

# Proportion Plot

<details>
<summary><i>details</i></summary>
Display stacked bar plots displaying groupings (x axis) and sub-groupings (stacks) of cells grouped by categorical data (i.e., Conditions with Clusters).  The y-axis can display either counts or proportions.
</details>

""")

if not adata:
    w_text_output(
      content="No data loaded...",
      appearance={"message_box": "warning"}
    )
    exit()

group_by = w_select(
    label="group by",
    default="sample",
    options=tuple(groups),
    appearance={
        "detail": "(cluster, sample, condition)",
        "help_text": "Select group to display on x-axis."
    }
)
stack_by = w_select(
    label="stack by",
    default="cluster",
    options=tuple(groups),
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

stacked_df = create_proportion_dataframe(
    adata, group_by.value, stack_by.value, return_type=return_type.value
)

w_row(
    items=[group_by, stack_by, return_type]
)

# Create proportion plot from the stacked dataframe created above
fig = px.bar(
    stacked_df,
    x="group_by",
    y="value",
    color="stack_by",
    barmode="stack",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    title=f"Distribution of {stack_by.value} by {group_by.value}"
)

# Update layout
fig.update_layout(
    xaxis_title=group_by.value,
    yaxis_title="Proportion" if return_type.value == "proportion" else "Count",
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=True,
    legend_title=stack_by.value
)

# Update axes
fig.update_xaxes(showgrid=False)
fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

fig.show()
