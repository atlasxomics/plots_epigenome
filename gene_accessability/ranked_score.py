w_text_output(content="""

# Ranked Plot (Gene Accessibility)

<details>
<summary><i>details</i></summary>

For the selected condition(s), create a ranked plot of gene scores between groups.

</details>

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
  gcompare_options = group_options["condition"]
except KeyError:
  w_text_output(
    content="No conditions found in experiment...",
    appearance={"message_box": "warning"}
  )
  exit()

gcompare_values = ["p_val", "p_val_adj", "Log2FC"]

gcompare_group_a = w_select(
    label="Condition A",
    default=None,
    options=tuple(gcompare_options),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "First condition for ranked feature plot."
    }
)

gcompare_group_b = w_select(
    label="Condition B",
    default=None,
    options=tuple(gcompare_options),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "Second group for ranked feature plot."
    }
)

gcompare_cluster = w_select(
    label="cluster",
    default="All",
    options=tuple(clusters + ["All"]),
    appearance={
        "help_text": "Filter data to a specific cluster."
    }
)

gcompare_rankby = w_select(
    label="Rank By",
    default="Log2FC",
    options=tuple(gcompare_values),
    appearance={
        "help_text": "Metric to rank plot by."
    }
)

gcompare_colorby = w_select(
    label="Color By",
    default="p_val",
    options=tuple(gcompare_values),
    appearance={
        "help_text": "Metric to color plot by."
    }
)


gcompare_n_genes = w_text_input(
  label="label top n genes",
  default="4",
  appearance={
    "help_text": "Set number of top/bottom genes to plot for each group."
  }
)

w_row(items=[
    gcompare_group_a,
    gcompare_group_b,
    gcompare_rankby,
    gcompare_colorby,
    gcompare_n_genes,
])

gcompare_button = w_button(label="Update Rank Plot")
if gcompare_group_a.value is not None and gcompare_group_b.value is not None and gcompare_button.value:
  
  # Check if groups have a value.
  for value in [gcompare_group_a.value, gcompare_group_b.value]:
      if value.__class__.__name__ in ["Nothing", "NoneType", "None", "Nothing.x"]:
          w_text_output(
            content="Please select groups for plotting.",
            appearance={"message_box": "info"}
          )
          submit_widget_state()
          exit(0)
  
  if gcompare_group_a.value == gcompare_group_b.value:
      w_text_output(
        content="Groups to compare must be different, please select different \
          groups.",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()      
      exit()
  
  gcompare_df = adata_g.uns[f"volcano_1_{gcompare_group_a.value}"]
  gcompare_df = gcompare_df[gcompare_df["cluster"] == gcompare_cluster.value]
  if len(gcompare_df) == 0:
      w_text_output(
         content=f"There is no volcano plot for cluster {gcompare_cluster.value} because it contains more than 90% of one of the conditions. Please check Proportion plot.",
         appearance={"message_box": "warning"}
      )
      exit(0)
  
  gcompare_rankby = gcompare_rankby.value
  if gcompare_rankby in ["p_val", "p_val_adj"]:
      gcompare_df[f"-log10{gcompare_rankby}"] = -np.log10(gcompare_df[gcompare_rankby])
      gcompare_rankby = f"-log10{gcompare_rankby}"
  
  fig_rank_plot_g = plot_ranked_feature_plotly(
      gcompare_df,
      y_col=gcompare_rankby,
      x_col=None,
      n_labels=int(gcompare_n_genes.value),
      label_col="gene",
      color_col=gcompare_colorby.value,
      colorscale="PuBu_r",
      marker_size=6,
      title=f"Differential genes: {gcompare_group_a.value} v. {gcompare_group_b.value}",
      y_label=gcompare_rankby,
  )
  
  w_plot(source=fig_rank_plot_g)
