w_text_output(content="""

# Compare Conditions (Motif Ranked Score Plot)

<details>
<summary><i>details</i></summary>

For the selected condition(s), create a ranked plot of motif scores between groups. The plot shows the raw z-scores on the y-axis, with the highest and lowest scoring motifs labeled in red.  The color of each dot represents the adjusted p-value.

<br>

As with the volcano plot above, if "All" is selected for "Condition B", all motifs will be compared to the union of all other conditions.  Otherwise, only the top 2,000 most variable motifs are included to speed up computation.

</details>

""")

if not adata_m:
    w_text_output(
       content="No motif data selected...",
       appearance={"message_box": "warning"}
    )
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif data loaded...",
       appearance={"message_box": "warning"}
    )
    exit()
if "condition" not in adata_m.obs.keys():
    w_text_output(
       content="No conditions found in experiment...",
       appearance={"message_box": "warning"}
    )
    exit()

mcompare_options = group_options["condition"]

mcompare_group_a = w_select(
    label="Condition A",
    default=None,
    options=tuple(mcompare_options),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "First condition for ranked motif plot."
    }
)

mcompare_group_b = w_select(
    label="Condition B",
    default=None,
    options=tuple(mcompare_options + ["All"]),
    appearance={
        "help_text": "You must click 'Run' after selecting both groups to run Cell.",
        "description": "Second group for ranked motif plot; if 'All', the selected group will be compared to all other groups."
    }
)

mcompare_n_motifs = w_text_input(
  label="label top n motifs",
  default="4",
  appearance={
    "help_text": "Set number of top/bottom motifs to plot for each group."
  }
)

w_row(items=[
    mcompare_group_a,
    mcompare_group_b,
    mcompare_n_motifs,
])

w_text_output(content="""
> If **"All"** is selected for **"group B"**, all motifs are shown. Otherwise, only the
**top 2,000 most variable motifs** are included to speed up computation.
""")

# Unsubscribe computation from widgets
mcompare_group_a_value = mcompare_group_a._signal.sample()
mcompare_group_b_value = mcompare_group_b._signal.sample()

# Check if groups have a value.
for value in [mcompare_group_a_value, mcompare_group_b_value]:
    if value.__class__.__name__ in ["Nothing", "NoneType", "None", "Nothing.x"]:
        w_text_output(
          content="Please select groups for plotting.",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        exit(0)

if mcompare_group_a_value == mcompare_group_b_value:
    w_text_output(
      content="Groups to compare must be different, please select different \
        groups.",
      appearance={"message_box": "warning"}
    )
    exit()

mcompare_key = f"{mcompare_group_a_value}_{mcompare_group_b_value}_True_motifs"
mcompare_wf_key = f"pairwise_{mcompare_group_a_value}_vs_{mcompare_group_b_value}"

if mcompare_key in mvol_cache.keys():
    mcompare_df = mvol_cache[mcompare_key]
elif mcompare_wf_key in adata_m.uns.keys():  # wf precomputed
    mcompare_df = sc.get.rank_genes_groups_df(adata_m, group=None, key=mcompare_wf_key)
    mvol_cache[mcompare_key] = mcompare_df
else:
    mcompare_df = make_volcano_df(
        adata_m,
        "condition",
        mcompare_group_a_value,
        mcompare_group_b_value,
        "motifs",
        0.0,  # no pval_adj filter
        True,
    )
    mvol_cache[mcompare_key] = mcompare_df

fig_rank_plot_m = plot_ranked_feature_plotly(
    mcompare_df,
    y_col="scores",
    x_col=None,
    n_labels=int(mcompare_n_motifs.value),
    label_col="names",
    color_col="pvals_adj",
    colorscale="PuBu_r",
    marker_size=6,
    title=f"Differential motifs: {mcompare_group_a_value} v. {mcompare_group_b_value}",
    y_label="Z-score"
)

fig_rank_plot_m.show()
