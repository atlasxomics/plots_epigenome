new_data_signal()

try:
    import leidenalg
except ImportError:
    # Try pip first
    print("Installing leidenalg...")
    exit_code = os.system("pip install leidenalg")
    
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

sig_updated_features()
if not chk_box_bansky.value:
    if method_group.sample() == "Banksy":
            
            if "params" not in locals():
                exit()
            
            
            w = w_workflow(
              wf_name = "wf.__init__.domain_detection_wf", 
              version ="0.0.0-3dc176-wip-bbc794",
              params = params, 
              label = "Run Banksy"
            )
            
            execution = w.value
            
            if execution is not None:
                res = await execution.wait()
            
                if res.status == 'SUCCEEDED':
                    local_path = "/root/adata_final_labels.csv"
                    print(output_dir.value)
                    adata_final_labels_latchpath = LPath(f"{output_dir.value.path}/{run_name.value}/adata_final_labels.csv")
                    adata_final_labels_latchpath.download(local_path)
                    
                    df_adata_labels = pd.read_csv(local_path)
                    df_adata_labels = df_adata_labels.set_index("barcode") 
                    
                    col = df_adata_labels.columns.tolist()[0]
                    print(result_label_input.value)
                    print(adata_g.obs.columns)
                    
                    if result_label_input.value in adata_g.obs.columns:
                        print("Here")
                        adata.obs.loc[df_adata_labels.index, f"{result_label_input.value}"] = df_adata_labels[col] 
                    else:
                        df_adata_labels = df_adata_labels.rename(columns = {col:result_label_input.value})
                        adata_g.obs = adata_g.obs.join(df_adata_labels, how = "left")
                    adata_g.obs[result_label_input.value] =  adata_g.obs[result_label_input.value].astype(str)
                    adata_g.obs[result_label_input.value] =  adata_g.obs[result_label_input.value].fillna("Null")
                    
                    print(adata_g.obs)  
                    sig_bansky(True)
    else:
        btn_domains = w_button(label = "Compute spatial domains")
    
        if btn_domains.value:
            import squidpy as sq
            adata_sub = adata_g[adata_g.obs[group_key.value] == group_val.value]
    
            print(adata_sub)
            sq.gr.spatial_neighbors(adata_sub, n_neighs = int(neighbors.value), coord_type="generic") 
            w_text_output(content = "[1/2] Computed Neighbors", appearance = {"message_box":"info"})
        
            sc.tl.leiden(adata_sub, 
                         resolution=float(threshold.value), 
                         neighbors_key="spatial_neighbors", 
                         key_added=f"spatial_domains_res_{threshold.value}")
            w_text_output(content = "[2/2] Computed spatial domains", appearance = {"message_box":"info"})
        
            df = adata_sub.obs[[f"spatial_domains_res_{threshold.value}"]]
            if result_label_input.value in adata_g.obs.columns:
                print("Here")
                col = f"spatial_domains_res_{threshold.value}"
                adata.obs.loc[df.index, f"{result_label_input.value}"] = df[col]
                
            else:
                col = f"spatial_domains_res_{threshold.value}"
                df = df.rename(columns = {col:result_label_input.value})
                adata_g.obs = adata_g.obs.join(df, how = "left")
            w_text_output(content = "Succesfully computed spatial domains using squidpy", appearance = {"message_box":"success"})
    
            sig_bansky(True)
            

