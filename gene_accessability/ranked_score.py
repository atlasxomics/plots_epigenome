w_text_output(content="""

# Compare Conditions (Gene Ranked Score Plot)

<details>
<summary><i>details</i></summary>

For the selected condition(s), create a ranked plot of gene scores between groups. The plot shows the raw z-scores on the y-axis, with the highest and lowest scoring genes labeled in red.  The color of each dot represents the adjusted p-value.

<br>

As with the volcano plot above, if "All" is selected for "Condition B", all genes will be compared to the union of all other conditions.  Otherwise, only the top 2,000 most variable genes are included to speed up computation.

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
if "condition" not in adata_g.obs.keys():
    w_text_output(
       content="No conditions found in experiment...",
       appearance={"message_box": "warning"}
    )
    exit()

gcompare_options = group_options["condition"]

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
    options=tuple(gcompare_options + ["All"]),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "Second group for ranked feature plot; if 'All', the selected group will be compared to all other groups."
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
    gcompare_n_genes,
])

w_text_output(content="""
> If **“All”** is selected for **“group B”**, all features are shown. Otherwise, only the
**top 2,000 most variable features** are included to speed up computation.
""")

# Unsubscribe computation from widgets
gcompare_group_a_value = gcompare_group_a._signal.sample()
gcompare_group_b_value = gcompare_group_b._signal.sample()

# Check if groups have a value.
for value in [gcompare_group_a_value, gcompare_group_b_value]:
    if value.__class__.__name__ in ["Nothing", "NoneType", "None", "Nothing.x"]:
        w_text_output(
          content="Please select groups for plotting.",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        exit(0)

if gcompare_group_a_value == gcompare_group_b_value:
    w_text_output(
      content="Groups to compare must be different, please select different \
        groups.",
      appearance={"message_box": "warning"}
    )
    exit()

gcompare_key = f"{gcompare_group_a_value}_{gcompare_group_b_value}_True_genes"
gcompare_wf_key = f"pairwise_{gcompare_group_a_value}_vs_{gcompare_group_b_value}"

if gcompare_key in gvol_cache.keys():
    print("using cache")
    gcompare_df = gvol_cache[gcompare_key]
elif gcompare_wf_key in adata_g.uns.keys():  # wf precomputed
    print("using pre-computed")
    try:
        gcompare_df = sc.get.rank_genes_groups_df(adata_g, group=None, key=gcompare_wf_key)
    except KeyError as e:
        print(f"{e}: missing data in precomputed, computing...")
        gcompare_df = make_volcano_df(
            adata_hvg,
            "condition",
            gcompare_group_a_value,
            gcompare_group_b_value,
            "genes",
            0.0,  # no pval_adj filter
            True,
        )
        gvol_cache[gcompare_key] = gcompare_df
else:
    print("computing ")
    gcompare_df = make_volcano_df(
        adata_hvg,
        "condition",
        gcompare_group_a_value,
        gcompare_group_b_value,
        "genes",
        0.0,  # no pval_adj filter
        True,
    )
    gvol_cache[gcompare_key] = gcompare_df

fig_rank_plot_g = plot_ranked_feature_plotly(
    gcompare_df,
    y_col="scores",
    x_col=None,
    n_labels=int(gcompare_n_genes.value),
    label_col="names",
    color_col="pvals_adj",
    colorscale="PuBu_r",
    marker_size=6,
    title=f"Differential genes: {gcompare_group_a_value} v. {gcompare_group_b_value}",
    y_label="Z-score"
)

fig_rank_plot_g.show()
