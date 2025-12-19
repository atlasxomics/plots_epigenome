new_data_signal()
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    exit()

if "choose_h5_data" in globals():
  if choose_h5_data.value is not None:
    adata_h5 = h5data_dict[choose_h5_data.value]
  else:
    adata_h5 = adata_g
else:
  adata_h5 = adata_g

viewer = w_h5(ann_data=adata_h5)
h5_viewer_signal(True)

w_text_output(content="""
<details>
<summary><i>H5 Viewer Instructions</i></summary>

**Loading data**  
- Click the **Select File** icon and choose a directory containing AnnData objects from the Latch Data module.  
- The directory should contain at least one of the following files:  
  - `adata_ge.h5ad`: a SnapATAC2 `AnnData` object with `.X` as a gene accessibility matrix.  
  - `adata_motifs.h5ad`: a SnapATAC2 `AnnData` object with `.X` as a motif deviation matrix.  
- BigWig files for cluster-, sample-, and condition-level groups should be saved in the output directory under subfolders named `[group]_coverages`.
- Loading large datasets into memory may take several minutes.  
- By default, compatible files are located in `latch:///snap_outs/[project_name]/`.
- If the notebook becomes frozen, try refreshing the browser tab or clicking the Run All In Tab b button in the Run All dropdown menu.

</details>
""")
