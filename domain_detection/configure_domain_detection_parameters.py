
w_text_output(content = """# Spatial Domain Detection  

Spatial domain detection is the process of dividing a tissue into distinct regions (domains) by combining **gene expression profiles** with **spatial context**.  
Unlike standard clustering in scRNA-seq, which groups cells solely by transcriptomic similarity, spatial domain detection also considers **where cells are located** in the tissue.  

- Domains often align with histological structures or functional zones (e.g., periportal vs. pericentral liver regions).  
- Methods range from fast, graph-based clustering (e.g., Squidpy + Leiden) to more sophisticated probabilistic models (e.g., Banksy, BayesSpace).  
- Useful for exploring **tissue architecture**, **zonation patterns**, and **disease-specific spatial organization**.  
""")

new_data_signal()
sig_bansky = Signal(False)
sig_updated_features = Signal(False)

if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

chk_box_bansky = w_checkbox(label = "I have an existing bansky segmentation I want to use")

if chk_box_bansky.value:
    banksy_outfile = output_dir = w_ldata_picker(
        label="Select csv file containing bansky segmentation results",
        key = "bansky csv file",
    )

    result_label_input = w_text_input(
        label="Result Label",
        default="",
        key="result_label",
        appearance={"help_text": "Enter a name for the result annotation"}
    )

    btn_add_annotation = w_button(label = "Add bansky annotations")

    if btn_add_annotation.value:
        local_path = "/root/adata_final_labels.csv"
        print(output_dir.value)
        adata_final_labels_latchpath = banksy_outfile.value
        adata_final_labels_latchpath.download(local_path)
        
        df_adata_labels = pd.read_csv(local_path)
        df_adata_labels = df_adata_labels.set_index("barcode") 
        
        col = df_adata_labels.columns.tolist()[0]
        print(result_label_input.value)
        print(adata_g.obs.columns)
        
        # Always handle as string to avoid Categorical dtype issues
        if result_label_input.value in adata_g.obs.columns:
            print("Here")
            adata_g.obs[result_label_input.value] = adata_g.obs[result_label_input.value].astype(str)
            adata_g.obs.loc[df_adata_labels.index, result_label_input.value] = df_adata_labels[col].astype(str)
        else:
            df_adata_labels = df_adata_labels.rename(columns={col: result_label_input.value})
            df_adata_labels[result_label_input.value] = df_adata_labels[result_label_input.value].astype(str)
            adata_g.obs = adata_g.obs.join(df_adata_labels, how="left")
    
        # Fill missing entries with "Null"
        adata_g.obs[result_label_input.value] = adata_g.obs[result_label_input.value].astype(str).fillna("Null")
        
        print(adata_g.obs)
        sig_bansky(True)

        
else:   
    
    # Filters Section
    w_text_output(content="""## Filter Cells  
    
Subset the `AnnData` object to retain only the cell populations relevant to the analysis.  
For example, you may restrict the dataset to **hepatocyte-like cells (HSEs)** when investigating liver-specific zonation patterns.  
    
    """)
    obs = adata_g.obs
    cat_bool_cols = obs.select_dtypes(include=["category", "bool", "object"]).columns.tolist()
    
    group_key = w_select(
        label="Filter Group",
        options=cat_bool_cols,
        default=cat_bool_cols[0]
    )
    
    filter_vals = adata_g.obs[group_key.value]
    keys = [ob for ob in filter_vals.unique() if not pd.isna(ob)]
            
    group_val = w_select(
        label="Filter Value",
        options=keys,
        default=keys[0]
    )
    w_row(items = [group_key, group_val])
    
    # Parameters Section
    w_text_output(content="## Parameters")
    threshold = w_text_input(
        label="Clustering Threshold",
         appearance = {"help_text" : "Resolution for performing clustering."},
        default="0.7"
    )
    neighbors = w_text_input(
        label="Number of Neighbors",
        appearance = {"help_text" : "Number of neighbors to construct spatial neighbor graph"},
        default="15"
    )
    w_row(items = [threshold, neighbors])
    
    # Output columns name
    w_text_output(content="""## Annotation Result
    """)
    result_label_input = w_text_input(
        label="Result Label",
        default="",
        key="result_label",
        appearance={"help_text": "Enter a name for the result annotation"}
    )
    
    w_text_output(content="""## Methods for Domain Detection
- **Banksy**: A slower but more accurate method for identifying spatial domains.  
- **Squidpy**: A faster method that constructs a nearest-neighbor graph from spatial coordinates and performs clustering. This approach is well-suited for quick sanity checks.  
     """)
    
    method_group = w_radio_group(
        label="Domain Detection Methods:",
        options=["Banksy", "Squidpy"],
        default="Squidpy"
    )
    
    if method_group.value == "Banksy":
        w_text_output(content="## Output")
        run_name = w_text_input(
            label="Run Name",
            default="Run-1",
            key = "bansky",
        )
        output_dir = w_ldata_picker(
            label="Select Output Directory",
            key = "bansky out",
        )
        
        w_row(items = [run_name, output_dir])
        save_data_checkbox = w_button(label = "Save Data")
        if save_data_checkbox.value:
            adata_sub = adata_g[adata_g.obs[group_key.value] == group_val.value]
            
            if output_dir.value is None:
                w_text_output(content = "Please select an output directory to proceed further", 
                    appearance = {"message_box":"warning"})
                exit()
        
            w_text_output(content = f"Uploading to subsetted data to {output_dir.value.path}. This could take a few seconds",
                    appearance = {"message_box":"info"})
            submit_widget_state()
            
            local_path = "/root/local_file.h5ad"
            print(output_dir.value.path)
            adata_sub.write(local_path, compression = "gzip")
            latch_path = LPath(f"{output_dir.value.path}/{run_name.value}.h5ad")
            latch_path.upload_from(Path(local_path))
            
            print(latch_path.path)
            
            w_text_output(content = f"Subsetted data and uploaded to {output_dir.value.path}/{run_name.value}.h5ad",
                    appearance = {"message_box":"success"})
            params = {
                "run_name": run_name.value,
                "input_file": LatchFile(latch_path.path),
                "output_dir": LatchDir(output_dir.value.path),
                "num_neighbours": int(neighbors.value),
                "resolutions": float(threshold.value),
            }
                    
    sig_updated_features(True)
    
