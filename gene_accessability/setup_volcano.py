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
