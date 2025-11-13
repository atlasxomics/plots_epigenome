if "gene_score_done_signal" not in globals():
    gene_score_done_signal = Signal(False)

new_data_signal()
# Ensure gene activity AnnData is loaded
if not adata_g or not isinstance(adata_g, AnnData):
    w_text_output(
        content="No gene activity data loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

w_text_output(content="# Gene Scoring and Cell Type Assignment")

# --- Step 0: Filter adata_g to create adata_subset based on a selected obs value ---
w_text_output(content="""
## Select a Subset of Cells

Gene scores can vary significantly across different conditions.  
To ensure meaningful results, it's best to subset your dataset to a single condition before computing scores.

Start by selecting the metadata column that defines your conditions,  
then choose the specific condition you'd like to focus on.
""")

# Widget to choose which obs column to filter by
filter_col = w_select(
    label="Filter by metadata column",
    default="condition" if "condition" in adata_g.obs else "",
    options=tuple(adata_g.obs.columns),
    appearance={"help_text": "Select an AnnData.obs column to filter cells."}
)

# Widget to choose the value within that column
if filter_col.value is not None:
    vals = adata_g.obs[filter_col.value].dropna().unique().tolist()
    filter_val = w_select(
        label=f"Filter value for {filter_col.value}",
        default=None,
        options=tuple(vals),
        appearance={"help_text": f"Select a value in {filter_col.value} to subset cells."}
    )
else:
    filter_val = None

# Require a filter value if a column is chosen
if filter_col.value and filter_val.value is None:
    w_text_output(
        content="Please select a filter value to subset cells.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

# Create the subset AnnData
if filter_col.value and filter_val.value is not None:
    adata_subset = adata_g[adata_g.obs[filter_col.value] == filter_val.value].copy()
else:
    adata_subset = adata_g.copy()

# --- Step 1: Text input for number of cell types ---
w_text_output(content="""
## Select Number of Cell Types
Specify how many cell types you want to predict, along with the marker gene sets for each. For each cell, the tool calculates a gene set score for every cell type and assigns the cell type with the highest score as the _predicted cell type_.
""")
n_types_input = w_text_input(
    label="Number of cell types",
    default="",
    key="no_cell_types_key",
    appearance={"help_text": "Enter an integer number of cell types"}
)

# --- Step 2: Parse the input into an integer ---
try:
    n_types = int(n_types_input.value)
except (ValueError, TypeError):
    n_types = 0

# Prepare lists to store widgets
label_inputs = []
feature_selects = []

# --- Step 3: Dynamically create widgets for each cell type ---
for idx in range(n_types):
    label_input = w_text_input(
        label=f"Cell type {idx+1} label",
        default="",
        appearance={"help_text": "Enter a name for this cell type"},
        key=f"Cell type {idx+1} label"
    )
    feature_select = w_multi_select(
        label=f"Select features for cell type {idx+1}",
        options=tuple(adata_g.var_names),
        appearance={"help_text": "Choose one or more genes for this cell type"},
        key=f"Select features for cell type {idx+1}"
    )
    label_inputs.append(label_input)
    feature_selects.append(feature_select)
    w_row(items=[label_input, feature_select])

# --- NEW: Text input to store overall result label ---
w_text_output(content="""
## Annotation Result
""")

print("Columns", [col for col in adata_g.obs.columns if pd.api.types.is_categorical_dtype(adata_g.obs[col])])

result_label_select = w_select(
    label="Annotation Label",
    key="result_label_dropdown",
    options=["Create New Label"] + adata_g.obs.select_dtypes(include=["category", "object"]).columns.tolist(),
    appearance={"help_text": "Select an existing label to store the annotations, or create a new label"}
)

if result_label_select.value == "Create New Label":
  result_label_input = w_text_input(
      label="Result Label",
      default="",
      key="result_label",
      appearance={"help_text": "Enter a name for the result annotation"}
  )

