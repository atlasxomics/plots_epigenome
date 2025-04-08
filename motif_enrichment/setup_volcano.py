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
    w_text_output(content="No motif data loaded...",
    appearance={"message_box": "warning"})
    exit()
if not isinstance(adata_m, anndata.AnnData):
    w_text_output(
       content="No motif data loaded...",
       appearance={"message_box": "warning"}
    )
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

# Reset gvol group values between grouping selections
mvol_group_a_value = None
mvol_group_b_value = None
