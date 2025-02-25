w_text_output(content="""

# Compare Conditions (Gene Volcano Plot)

<details>
<summary><i>details</i></summary>

Visualize differential gene accessibility between groups.

<br>

First, select the grouping you are interested in (cluster, sample, condition), then the group (i.e., Cluster 1). The plot will compare the selected group to a union of all other groups in the grouping (i.e., Cluster 1 versus all other clusters).

<br>

The volcano plot displayes the log2 fold change on the x-axis and the negative log10 of the adjusted p-value on the y-axis. The p-value is adjusted to account for the False Discovery Rate; see scanpy documentation for more details.

<br>

_We are working to add the ability to compare specific groups (i.e., Cluster 1 versus Cluster 2) and filter groups by other groups by other metadata (i.e., Cluster 1-health)._

</details>

""")

if not adata_g:
    w_text_output(
        content="No data gene activity data loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

w_text_output(content="Select the grouping (cluster, sample, condition) you are interested in comparing.")

gvol_grouping = w_select(
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

if not adata_g:
    w_text_output(
        content="No data gene activity data loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

try:
    gvol_grouping
except NameError:
    w_text_output(
        content="Please select 'grouping' in the Cell above first.",
        appearance={"message_box": "warning"}
    )
    exit()

gvol_options = group_options[gvol_grouping.value]

gvol_group_a = w_select(
    label="group A",
    default=None,
    options=tuple(gvol_options),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "First group for volcano plot selection."
    }
)

gvol_group_b = w_select(
    label="group B",
    default="All",
    options=tuple(gvol_options + ["All"]),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "Second group for volcano plot selection; if 'All', the selected group will be compared to all other groups."
    }
)

pvals_adj_threshold = w_text_input(
  label="pval adjust threshold",
  default="0.01",
)

log2fc_threshold = w_text_input(
  label="log2fc threshold",
  default="0.01",
)

gvol_width = w_text_input(
  label="plot width",
  default="1000",
)

gvol_height = w_text_input(
  label="plot height",
  default="425",
)

w_row(items=[
    gvol_group_a,
    gvol_group_b,
    pvals_adj_threshold,
    log2fc_threshold,
    gvol_width,
    gvol_height
])

# Unsubscribe computation from widgets
gvol_group_a_value = gvol_group_a._signal.sample()
gvol_group_b_value = gvol_group_b._signal.sample()

# Check if groups have a value.
for value in [gvol_group_a_value, gvol_group_b_value]:
    if value.__class__.__name__ in ["Nothing", "NoneType", "None"]:
        w_text_output(
          content="Please select groups for plotting.",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        exit(0)

if gvol_group_a_value == gvol_group_b_value:
    w_text_output(
      content="Groups to compare must be different, please select different \
        groups.",
      appearance={"message_box": "warning"}
    )
    exit()

gvol_key = f"{gvol_group_a_value}_{gvol_group_b_value}_genes"

if gvol_key in gvol_cache.keys():
    gvol_df = gvol_cache[gvol_key]
else:
    gvol_df = make_volcano_df(
        adata_g,
        gvol_grouping.value,
        gvol_group_a_value,
        gvol_group_b_value,
        "genes"
    )
    gvol_cache[gvol_key] = gvol_df


fig_volcano_plot = plot_volcano(
  gvol_df,
  float(pvals_adj_threshold.value),
  float(log2fc_threshold.value),
  gvol_group_a_value,
  gvol_group_b_value,
  int(gvol_width.value),
  int(gvol_height.value)
)

# Show the plot
fig_volcano_plot.show()
