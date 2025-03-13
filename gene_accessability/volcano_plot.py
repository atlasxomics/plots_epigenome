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
    default=None,
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

gvol_display0 = w_checkbox(
  label="display 0 p-vals",
  default=True,
  appearance={
    "description": "Whether to display features for with the p-value evaluates to 0."
  }
)

gvol_displayGm = w_checkbox(
  label="display pseudogenes (mouse only)",
  default=True,
  appearance={
    "description": "Whether to display pseudogenes; for mouse data only."
  }
)

w_row(items=[gvol_display0, gvol_displayGm])

w_text_output(content="""
> If **“All”** is selected for **“group B”**, all features are shown. Otherwise, only the
**top 2,000 most variable features** are included to speed up computation.

> Some p-values are so small that they round to 0 due to numerical limits.
To display these features on the volcano plot, we replace them with a very small value.
As a result, these points may appear as a horizontal line at the top of the plot.  They
can be toggeled with the 'display 0 p-vals' button.
""")

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

gvol_key = f"{gvol_group_a_value}_{gvol_group_b_value}_filter-{gvol_display0.value}_{gvol_displayGm.value}_genes"

if gvol_key in gvol_cache.keys():
    gvol_df = gvol_cache[gvol_key]
else:
    gvol_df = make_volcano_df(
        adata_hvg,
        gvol_grouping.value,
        gvol_group_a_value,
        gvol_group_b_value,
        "genes",
        float(pvals_adj_threshold.value),
        gvol_display0.value,
        gvol_displayGm.value
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
