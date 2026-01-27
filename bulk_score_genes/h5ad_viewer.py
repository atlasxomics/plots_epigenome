new_data_signal()
gene_score_done_signal()

# Ensure gene activity AnnData is loaded
if not adata_g or not isinstance(adata_g, AnnData):
    w_text_output(
        content="No gene activity data loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

w_h5(ann_data=adata_g)

adata_g_obs_df = adata_g.obs
w_table(label="Cell Metadata", source=adata_g_obs_df)

