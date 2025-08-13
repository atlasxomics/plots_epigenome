new_data_signal()

if not adata:
  w_text_output(
    content="No data data loaded...",
    appearance={"message_box": "warning"}
  )
  exit()

wf_results_signal()
if wf_results_signal.sample() == True:
  
  if "motif" in results_dict.keys():
  
    w_text_output(content="""## Differential Motif Enrichment""")
    
    m_group = w_select(
        label="Group",
        default="GroupA",
        options=tuple(["GroupA", "GroupB"]),
    )
    
    m_pvals_adj_threshold = w_text_input(
      label="pval adjust threshold",
      default="0.05",
    )
    
    m_meandiff_threshold = w_text_input(
      label="MeanDiff threshold",
      default="0.01",
    )
    
    mcompare_rankby = w_select(
        label="Rank By",
        default="MeanDiff",
        options=tuple(['FDR', 'MeanDiff']),
