w_text_output(content="## Track Browser")

w_text_output(content="""
<details>
<summary><i>Track Browser Instructions</i></summary>

This Cell visuallizes ATAC-seq coverage for a specified grouping (cluster, sample, condition) with the IGV-Web App.
To initialize the app, select a grouping and genome below.
Please see the [IGV-Web App User Guide](https://igv.org/doc/webapp/#UserGuide/) for more information.

</details>
""")

new_data_signal()

if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    exit()

coverages_genome = w_select(
    label="Genome",
    options=("hg38", "mm10"),
    key="coverages_genome",
    default="hg38",
    appearance={
      "help_text": "Select reference genome."
    }
)

coverages_group = w_select(
    label="Coverage Group",
    options=tuple(coverages_dict.keys()),
    key="coverages_group",
    default="sample",
    appearance={
      "help_text": "Select grouping for coverage tracks."
    }
    
)

w_row(items=[coverages_genome, coverages_group])
