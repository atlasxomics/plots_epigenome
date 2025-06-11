coverages_group = w_select(
    label="Coverage group",
    options=tuple(coverages_dict.keys())
)

coverages_gene = w_select(
    label="Choose a gene to view",
    options=available_genes,
)

coverages_button = w_button(label="Update IGV Viewer")

input_column = w_column(items=[coverages_group, coverages_gene, coverages_button])

if coverages_gene.value is not None and coverages_group.value is not None and coverages_button.value:

    coverages_dir = coverages_dict[coverages_group.value]

    bedgraphs = []
    bigwigs = []
    for file in coverages_dir.iterdir():
        suffix = file.path.split(".")[-1]
        if suffix == "gz":
            bedgraphs.append(file)
        elif suffix == "bw":
            bigwigs.append(file)

    files = []
    if len(bedgraphs) + len(bigwigs) == 0:
        w_text_output(
            content="No bedgraph (.gz) or bigwig (.bw) files found in coverages directory...",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()
    elif len(bedgraphs) > 0 and len(bigwigs) == 0:
        w_text_output(
            content="Using bigwigs instead of bedgraph could make IGV much faster...",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        format = "bedgraph"
        files = bedgraphs
    elif len(bedgraphs) == 0 and len(bigwigs) > 0:
        format = "bigwig"
        files = bigwigs
    elif len(bedgraphs) > 0 and len(bigwigs) > 0:
        w_text_output(
            content="Both bedgraph (.gz) andd bigwig (.bw) files found in coverages directory; please ensure only one type is present...",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

    tracks = []
    for f in files:
        tracks.append({
            "name": f.path.split("/")[-1],
            "type": "wig",
            "format": format,
            "url": f.path,
            "autoscale": True,
            "visibilityWindow": 100000,
        })

    selected_gene = coverages_gene.value

    opts: IGVOptions = {
        "genome": "hg38",
        "locus": selected_gene,
        "tracks": tracks
    }

    igviewer = w_igv(options=opts)
