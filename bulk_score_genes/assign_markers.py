w_text_output(content="# Gene Scoring and Cell Type Assignment")

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

# --- Step 0: Filter adata_g to create adata_subset based on a selected obs value ---
w_text_output(content="""
## Select a Subset of Cells

Gene scores can vary across conditions.
Subseting your dataset to a single condition before computing scores could improve the data clarity.

Start by selecting the metadata column that defines your groups,
then choose the specific group you'd like to focus on.
""")

# Widget to choose which obs column to filter by

filter_groups = [
  key for key in adata.obs_keys() if
  (pd.api.types.is_object_dtype(adata.obs[key]) or pd.api.types.is_categorical_dtype(adata.obs[key]))
  and key != "cluster"
]
filter_groups = filter_groups + ["None"]

filter_col = w_select(
    label="Filter by metadata column",
    default=None,
    options=tuple(filter_groups),
    appearance={"help_text": "Select an categorical observation to filter cells."}
)

use_filter = False if filter_col.value == None or filter_col.value == "None" else True


# Widget to choose the value within that column
if use_filter:
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
if use_filter and filter_val.value is None:
    w_text_output(
        content="Please select a filter value to subset cells.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

# Create the subset AnnData
if use_filter and filter_val.value is not None:
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

# Text input to store overall result label
w_text_output(content="""
## Annotation Result
""")

result_labels = [
  col for col in adata_g.obs.select_dtypes(include=["category", "object"]).columns
  if col not in groups  # Don't want to overwrite cannonical groups.
]

result_label_select = w_select(
    label="Annotation Label",
    key="result_label_dropdown",
    default="Create New Label",
    options=tuple(["Create New Label"] + result_labels),
    appearance={"help_text": "Select an existing label to store the annotations, or create a new label"}
)

if result_label_select.value == "Create New Label":
  result_label_input = w_text_input(
      label="Result Label",
      default="",
      key="result_label",
      appearance={"help_text": "Enter a name for the result annotation"}
  )

# --- Validation and gating for downstream steps ---
validation_errors = []

if n_types <= 0:
    validation_errors.append("Number of cell types must be greater than 0.")

for idx in range(n_types):
    label_value = (label_inputs[idx].value or "").strip()
    if not label_value:
        validation_errors.append(f"Cell type {idx + 1} label cannot be empty.")

    features = feature_selects[idx].value
    if not features:
        validation_errors.append(f"Select at least one feature for cell type {idx + 1}.")

if result_label_select.value == "Create New Label":
    new_label = (result_label_input.value or "").strip() if result_label_input else ""
    if not new_label:
        validation_errors.append("Please enter a name for the new Annotation Label.")
else:
    if not result_label_select.value:
        validation_errors.append("Please select an Annotation Label.")

if validation_errors:
    w_text_output(
        content="Please resolve the following before continuing:\n" + "\n".join(
            f"- {msg}" for msg in validation_errors
        ),
        appearance={"message_box": "warning"}
    )
else:
    choose_group_signal(True)

