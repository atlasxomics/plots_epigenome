w_text_output(content="""

# Volcano Plot (Motif Deviation)

<details>
<summary><i>details</i></summary>
</details>

""")


if not adata_m:
    w_text_output(
        content="No motif data activity data selected...",
        appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif data loaded...",
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

clusters = group_options["cluster"]

mvol_condition = w_select(
    label="condition",
    default=conditions[0],
    options=tuple(conditions),
    appearance={
        "help_text": "Select condition for comparison."
    }
)

mvol_cluster = w_select(
    label="cluster",
    default="All",
    options=tuple(clusters + ["All"]),
    appearance={
        "help_text": "Filter data to a specific cluster."
    }
)

mvol_pvals_adj_threshold = w_text_input(
  label="pval adjust threshold",
  default="0.01",
)

mvol_log2fc_threshold = w_text_input(
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
    mvol_condition,
    mvol_cluster,
    mvol_pvals_adj_threshold,
    mvol_log2fc_threshold,
])

w_row(items=[mvol_width, mvol_height])

mvol_button = w_button(label="Update Volcano Plot")

if mvol_condition.value is not None and mvol_button.value:

  mvol_df = adata_m.uns[f"volcano_1_{mvol_condition.value}"]
  mvol_df = mvol_df[mvol_df["cluster"] == mvol_cluster.value]
  try:
    mvol_df.drop(["Significance"], axis=1, inplace=True)
  except:
    print("No Significance column found")
  if len(mvol_df) == 0:
      w_text_output(
         content=f"There is no volcano plot for cluster {gvol_cluster.value} because it contains more than 90% of one of the conditions. Please check Proportion plot.",
         appearance={"message_box": "warning"}
      )
      exit(0)


  mvol_volcano_plot = plot_volcano(
      mvol_df,
      float(mvol_pvals_adj_threshold.value),
      float(mvol_log2fc_threshold.value),
      mvol_condition.value,
      "rest",
      pval_key="p_val",
      l2fc_key="MeanDiff",
      names_key="gene",
      plot_width=int(mvol_width.value),
      plot_height=int(mvol_height.value),
      top_n=2
  )
  
  w_plot(source=mvol_volcano_plot)
  w_table(source=mvol_df)
