w_text_output(content="""

# Compare Conditions (Gene Activity Volcano Plot)

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

try:
  conditions = group_options["condition"]
except KeyError:
  w_text_output(
    content="No conditions found in experiment...",
    appearance={"message_box": "warning"}
  )
  exit()

gvol_condition = w_select(
    label="condition",
    default=conditions[0],
    options=tuple(conditions),
    appearance={
        "help_text": "Select condition for comparison."
    }
)

gvol_cluster = w_select(
    label="cluster",
    default="All",
    options=tuple(clusters + ["All"]),
    appearance={
        "help_text": "Filter data to a specific cluster."
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
    gvol_condition,
    gvol_cluster,
    pvals_adj_threshold,
    log2fc_threshold,
])

w_row(items=[gvol_width, gvol_height])

gvol_df = adata.uns[f"volcano_1_{gvol_condition.value}"]
gvol_df = gvol_df[gvol_df["cluster"] == gvol_cluster.value]
if len(gvol_df) == 0:
    w_text_output(
       content=f"There is no volcano plot for cluster {gvol_cluster.value} because it contains more than 90% of one of the conditions. Please check Proportion plot.",
       appearance={"message_box": "warning"}
    )
    exit(0)


fig_volcano_plot = plot_volcano(
    gvol_df,
    float(pvals_adj_threshold.value),
    float(log2fc_threshold.value),
    gvol_condition.value,
    "rest",
    pval_key="p_val",
    l2fc_key="avg_log2FC",
    names_key="gene",
    plot_width=int(gvol_width.value),
    plot_height=int(gvol_height.value),
    top_n=2
)

# # # Show the plot
fig_volcano_plot.show()
