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

group_options = dict()
for group in groups:
  group_options[group] = list(adata_g.obs[group].unique())

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
    w_text_output(content="Please select 'grouping' in the Cell above first.", appearance={"message_box": "warning"})
    exit()

gvol_group = w_select(
  label="subgroup",
  options=tuple(group_options[gvol_grouping.value]),
  appearance={
    "help_text": "Select group for volcano plot; by default, the selected group will be compared to all other groups."
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

w_row(items=[gvol_group, pvals_adj_threshold, log2fc_threshold, gvol_width, gvol_height])

# Volcano plot
gvol_df = sc.get.rank_genes_groups_df(
  adata_g,
  group=gvol_group.value,
  key=f"{gvol_grouping.value}_genes"
)

# Remove rows with 'pvals_adj' == 0 from the dataframe
gvol_df_filtered = gvol_df[gvol_df['pvals_adj'] != 0]

fig_volcano_plot = plot_volcano(
  gvol_df_filtered,
  float(pvals_adj_threshold.value),
  float(log2fc_threshold.value),
  gvol_group.value,
  gvol_grouping.value,
  int(gvol_width.value),
  int(gvol_height.value)
)

# Show the plot
fig_volcano_plot.show()
