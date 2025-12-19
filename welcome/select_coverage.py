new_data_signal()

if not adata_g:
    w_text_output(
        content=" ",
    )
    exit()

tracks_signal()
coverages_group = w_select(
    label="Coverage Group",
    options=tuple(coverages_dict.keys()),
    default="sample",
    appearance={
      "help_text": "Select grouping for coverage tracks."
    }
    
)