new_data_signal()
sig_bansky()

if not adata:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()
viewer = w_h5(ann_data=adata_g)
adata_g_obs_df = adata_g.obs
w_table(label="Cell Metadata", source=adata_g_obs_df)

