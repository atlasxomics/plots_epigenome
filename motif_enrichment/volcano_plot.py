w_text_output(content=f"""

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

group_options = dict()
for group in groups:
  group_options[group] = list(adata_g.obs[group].unique())

w_text_output(content="Select the grouping (cluster, sample, condition) you are interested in comparing.")

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
  w_text_output(content="No motif data loaded...",  appearance={"message_box": "warning"})
  exit()

mvol_group = w_select(
  label="subgroup",
  options=tuple(group_options[mvol_grouping.value]),
  appearance={
    "help_text": "Select group for volcano plot; by default, the selected group will be compared to all other groups."
  }
)

m_pvals_adj_threshold = w_text_input(
  label="pval adjust threshold",
  default="0.01",
  appearance={
    "help_text": " \n"
  }
)

m_log2fc_threshold = w_text_input(
  label="log2fc threshold",
  default="0.01",
  appearance={
    "help_text": " "
  }
)

mvol_width = w_text_input(
  label="plot width",
  default="800",
  appearance={
    "help_text": " "
  }
)

mvol_height = w_text_input(
  label="plot height",
  default="600",
  appearance={
    "help_text": " "
  }
)

w_row(items=[mvol_group, m_pvals_adj_threshold, m_log2fc_threshold, mvol_width, mvol_height])

mvol_df = sc.get.rank_genes_groups_df(
  adata_m,
  group=mvol_group.value,
  key=f"{mvol_grouping.value}_motifs"
)

# Remove rows with 'pvals_adj' == 0 from the dataframe
mvol_df_filtered = mvol_df[mvol_df['pvals_adj'] != 0]

fig_volcano_plot_m = plot_volcano(
  mvol_df_filtered,
  float(m_pvals_adj_threshold.value),
  float(m_log2fc_threshold.value),
  mvol_group.value,
  mvol_grouping.value,
  int(mvol_width.value),
  int(mvol_height.value)
)


# Show the plot
fig_volcano_plot_m.show()