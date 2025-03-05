if not adata_m:
    w_text_output(
       content="No motif data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()

try:
    mvol_grouping
except NameError:
    w_text_output(
        content="Please select 'grouping' in the Cell above first.",
        appearance={"message_box": "warning"}
    )
    exit()

mvol_options = group_options[mvol_grouping.value]

mvol_group_a = w_select(
    label="group A",
    options=tuple(mvol_options),
    default=None,
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "First group for volcano plot selection."
    }
)

mvol_group_b = w_select(
  label="group B",
  default=None,
  options=tuple(mvol_options + ["All"]),
  appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "Second group for volcano plot selection; if 'All', the selected group will be compared to all other groups."
  }
)

m_pvals_adj_threshold = w_text_input(
  label="pval adjust threshold",
  default="0.01",
)

m_log2fc_threshold = w_text_input(
  label="log2fc threshold",
  default="0.01",
)

mvol_width = w_text_input(
  label="plot width",
  default="1000",
)

mvol_height = w_text_input(
  label="plot height",
  default="425",
)

w_row(items=[
    mvol_group_a,
    mvol_group_b,
    m_pvals_adj_threshold,
    m_log2fc_threshold,
    mvol_width,
    mvol_height
])

mvol_display0 = w_checkbox(
  label="display 0 p-vals",
  default=True,
  appearance={
    "description": "Whether to display features for with the p-value evaluates to 0."
  }
)

w_text_output(content="""
> Some p-values are so small that they round to 0 due to numerical limits.
To display these features on the volcano plot, we replace them with a very small value.
As a result, these points may appear as a horizontal line at the top of the plot.  They
can be toggeled with the 'display 0 p-vals' button.
""")

# Unsubscribe computation from widgets
mvol_group_a_value = mvol_group_a._signal.sample()
mvol_group_b_value = mvol_group_b._signal.sample()

# Check if groups have a value.
for value in [mvol_group_a_value, mvol_group_b_value]:
    print(value, value.__class__.__name__ )
    if value.__class__.__name__ in ["Nothing", "NoneType", "None"]:
        w_text_output(
          content="Please select groups for plotting.",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        exit(0)

if mvol_group_a_value == mvol_group_b_value:
    w_text_output(
      content="Groups to compare must be different, please select different \
        groups.",
      appearance={"message_box": "warning"}
    )
    exit()

mvol_key = f"{mvol_group_a_value}_{mvol_group_b_value}_filter-{mvol_display0.value}_motifs"

if mvol_key in mvol_cache.keys():
    mvol_df = mvol_cache[mvol_key]
else:
    mvol_df = make_volcano_df(
        adata_m,
        mvol_grouping.value,
        mvol_group_a_value,
        mvol_group_b_value,
        "motifs",
        float(m_pvals_adj_threshold.value),
        mvol_display0.value
    )
    mvol_cache[mvol_key] = mvol_df

fig_volcano_plot_m = plot_volcano(
  mvol_df,
  float(m_pvals_adj_threshold.value),
  float(m_log2fc_threshold.value),
  mvol_group_a_value,
  mvol_group_b_value,
  int(mvol_width.value),
  int(mvol_height.value)
)

# Show the plot
fig_volcano_plot_m.show()
