new_data_signal()

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

coverages_group = w_select(
    label="Coverage group",
    options=tuple(coverages_dict.keys()),
    appearance={
      "help_text": "Select grouping for coverage tracks."
    }
)

coverages_gene = w_select(
    label="Choose a gene to view",
    options=available_genes,
    appearance={
      "help_text": "Specify gene for coverage tracks."
    }
)

coverages_min = w_text_input(
    label="Track min (optional)",
    default="",
    appearance={
      "help_text": "Fix y-axis min for tracks"
    }
)

coverages_max = w_text_input(
    label="Track max (optional)",
    default="",
    appearance={
      "help_text": "Fix y-axis max for tracks"
    }
)

coverages_button = w_button(label="Update IGV Viewer")

input_row = w_row(items=[coverages_gene, coverages_min, coverages_max])
input_column = w_column(items=[coverages_group, input_row, coverages_button])

if coverages_min.value != "":
  try:
    float(coverages_min.value)
  except ValueError:
    w_text_output(
      content="Could not convert 'Track min' to float; please input a valid number.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()

if coverages_max.value != "":
  try:
    float(coverages_max.value)
  except ValueError:
    w_text_output(
      content="Could not convert 'Track max' to float; please input a valid number.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()

if coverages_gene.value is not None and coverages_group.value is not None and coverages_button.value:
    print(coverages_min.value, coverages_max.value)

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
        submit_widget_state()
        format = "bigwig"
        files = bigwigs

    tracks = []
    count = 0
    for f in files:
        tracks.append({
            "name": f.path.split("/")[-1],
            "type": "wig",
            "format": format,
            "url": f.path,
            "autoscale": False,
            "visibilityWindow": 100000,
        })
        if coverages_min.value != "":
          try:
            tracks[count]["min"] = float(coverages_min.value)
          except:
            print("Can't convert min")
            continue
        if coverages_max.value != "":
          try:
            tracks[count]["max"] = float(coverages_max.value)
          except:
            print("Can't convert max")
            continue

        count += 1
    print(tracks)

    selected_gene = coverages_gene.value

    opts: IGVOptions = {
        "genome": "hg38",
        "locus": selected_gene,
        "tracks": tracks,
    }

    igviewer = w_igv(options=opts)
