w_text_output(content="## Coverage Track Settings")

new_data_signal()

if not adata_g:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

coverages_genome = w_select(
    label="Genome",
    options=("hg38", "mm10"),
    default="hg38",
    appearance={
      "help_text": "Select reference genome."
    }
)

coverages_group = w_select(
    label="Coverage Group",
    options=tuple(coverages_dict.keys()),
    default="sample",
    appearance={
      "help_text": "Select grouping for coverage tracks."
    }
    
)

tracks_signal(True)