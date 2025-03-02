w_text_output(content="""

# Compare Conditions (Motif Volcano Plot)

<details>
<summary><i>details</i></summary>

Visualize differential activation of motifs between groups.

<br>

First, select the grouping you are interested in (cluster, sample, condition), then the group (i.e., Cluster 1). The plot will compare the selected group to a union of all other groups in the grouping (i.e., Cluster 1 versus all other clusters).

<br>

The volcano plot displayes the log2 fold change on the x-axis and the negative log10 of the adjusted p-value on the y-axis. The p-value is adjusted to account for the False Discovery Rate; see scanpy documentation for more details.

<br>

_We are working to add the ability to compare specific groups (i.e., Cluster 1 versus Cluster 2) and filter groups by other groups by other metadata (i.e., Cluster 1-health)._

""")

if not adata_m:
    w_text_output(content="No motif data loaded...",  appearance={"message_box": "warning"})
    exit()

w_text_output(
    content="Select the grouping (cluster, sample, condition) you are interested in comparing."
)

mvol_grouping = w_select(
    label="grouping",
    default="cluster",
    options=tuple(groups),
    appearance={
        "help_text": "Select categorical grouping for comparison."
    }
)

w_text_output(
    content="Group selected for comparison; navigate to Cell below to create a Volcano Plot.",
    appearance={"message_box": "info"}
)

###############################################################################

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

mvol_key = f"{mvol_group_a_value}_{mvol_group_b_value}_motifs"

if mvol_key in mvol_cache.keys():
    mvol_df = mvol_cache[mvol_key]
else:
    mvol_df = make_volcano_df(
        adata_m,
        mvol_grouping.value,
        mvol_group_a_value,
        mvol_group_b_value,
        "motifs"
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
