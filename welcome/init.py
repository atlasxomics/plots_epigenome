import math
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc

from anndata import AnnData

from lplots import palettes, submit_widget_state
from lplots.reactive import Signal
from lplots.widgets.button import w_button
from lplots.widgets.checkbox import w_checkbox
from lplots.widgets.column import w_column
from lplots.widgets.grid import w_grid
from lplots.widgets.h5 import w_h5
from lplots.widgets.igv import IGVOptions, w_igv
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.plot import w_plot
from lplots.widgets.row import w_row
from lplots.widgets.select import w_select
from lplots.widgets.table import w_table
from lplots.widgets.text import w_text_input, w_text_output

w_text_output(content="# **ATX Spatial Epigenomics Report**")
w_text_output(content="""

This notebook provides interactive tools for **exploratory data analysis** and **figure generation** from spatial epigenomic DBiT-seq experiments. Plotting modules are organized into tabs at the **top of this window**--move between tabs to explore results.

""")

DEFAULT_H5_CATEGORICAL_PALETTE = [
    "#C33530", "#282E66", "#43884A", "#7E2F8A", "#E48341",
    "#FAE64D", "#8E9ECD", "#B570A8", "#E0C3DA", "#9FD3E2",
    "#96C56C", "#E38180", "#9584B9", "#C25434", "#63B9A8",
    "#694D99", "#33707A", "#731F1C", "#D0A970", "#3D3D3D",
]
DEFAULT_CATEGORICAL_PALETTE_NAME = "Default H5 Viewer Palette"

if "new_data_signal" not in globals():
    new_data_signal = Signal(False)
if "refresh_h5_signal" not in globals():
    refresh_h5_signal = Signal(False)

na_keys = ["barcode", "on_off", "row", "col", "xcor", "ycor", "score"]

adata_g = None
adata_m = None
adata_g_path = None
adata_m_path = None
available_genes = []
available_motifs = []
samples = []
groups = []
coverages_dict = {}
h5data_dict = {}
adata_h5 = None
loaded_h5_data_key = None


def create_proportion_dataframe(
    adata: AnnData,
    group_by: str,
    stack_by: str,
    return_type: str = "proportion",
) -> pd.DataFrame:
    """Build stacked-bar input data from AnnData obs columns."""
    count_df = pd.crosstab(adata.obs[group_by], adata.obs[stack_by])

    if return_type == "proportion":
        result_df = count_df.div(count_df.sum(axis=1), axis=0)
    elif return_type == "counts":
        result_df = count_df
    else:
        raise ValueError("Invalid return_type. Use 'proportion' or 'counts'.")

    long_df = result_df.reset_index().melt(
        id_vars=group_by,
        value_name="value",
        var_name=stack_by,
    )
    long_df.columns = ["group_by", "stack_by", "value"]
    return long_df


def create_violin_data(
    adata: AnnData,
    group_by: str,
    plot_data: str,
    data_type: str = "obs",
) -> pd.DataFrame:
    """Create a long-form dataframe for box/violin plotting."""
    if data_type == "obs":
        values = adata.obs[plot_data]
    elif data_type in {"gene", "motif"}:
        values = adata[:, plot_data].X.flatten()
    else:
        raise ValueError("data_type must be either 'obs', 'gene', or 'motif'.")

    return pd.DataFrame({
        "group": adata.obs[group_by],
        "value": values,
    })


def drop_obs_column(adata_objects, col_to_drop: str = "orig.ident") -> None:
    """Remove a column from obs for a list of AnnData objects."""
    for adata_obj in adata_objects:
        if adata_obj is not None and col_to_drop in adata_obj.obs.columns:
            adata_obj.obs = adata_obj.obs.drop(columns=[col_to_drop])


async def get_notebook_palettes():
    try:
        palette_data = await palettes.get()
    except Exception as e:
        print(f"Unable to load notebook palettes: {e}")
        return {"categorical": [], "continuous": []}

    if not isinstance(palette_data, dict):
        return {"categorical": [], "continuous": []}

    cleaned_palettes = {}
    for kind in ("categorical", "continuous"):
        cleaned_palettes[kind] = [
            palette
            for palette in palette_data.get(kind, [])
            if palette.get("display_name") and palette.get("colors")
        ]

    return cleaned_palettes


def get_palette_selector_options(
    palette_data,
    kind: str = "categorical",
    fallback_name: str = DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    palette_names = []
    seen = {fallback_name}

    for palette in palette_data.get(kind, []):
        display_name = palette["display_name"]
        if display_name not in seen:
            palette_names.append(display_name)
            seen.add(display_name)

    return tuple([fallback_name] + palette_names)


def get_selected_palette_colors(
    palette_data,
    selected_name,
    kind: str = "categorical",
    fallback_colors=None,
    fallback_name: str = DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    if fallback_colors is None:
        fallback_colors = DEFAULT_H5_CATEGORICAL_PALETTE

    if not selected_name or selected_name == fallback_name:
        return fallback_colors

    for palette in palette_data.get(kind, []):
        if palette["display_name"] == selected_name:
            return palette["colors"]

    return fallback_colors


def build_discrete_color_map(categories, colors):
    if not colors:
        colors = DEFAULT_H5_CATEGORICAL_PALETTE

    return {category: colors[i % len(colors)] for i, category in enumerate(categories)}


def get_groups(adata: AnnData) -> List[str]:
    """Return categorical grouping columns relevant to the remaining tabs."""
    kept_groups = []
    for group in ["cluster", "sample", "condition"]:
        try:
            if len(adata.obs[group].unique()) > 1:
                kept_groups.append(group)
        except KeyError as e:
            print(f"{e}")

    return kept_groups


def process_matrix_layout(
    adata_all: AnnData,
    n_rows: int,
    n_cols: int,
    sample_key: str = "sample",
    spatial_key: str = "spatial",
    new_obsm_key: str = "X_dataset",
    tile_spacing: float = 100.0,
    flipy: bool = False,
    sample_order_mode: str = "original",
    condition_key: str = "condition",
) -> None:
    """Add a tiled spatial layout to `.obsm` for H5 viewing."""
    sample_series = adata_all.obs[sample_key].astype(str)

    if sample_order_mode == "original":
        ordered_samples = list(pd.unique(sample_series))
    elif sample_order_mode == "sample":
        ordered_samples = sorted(sample_series.unique().tolist())
    elif sample_order_mode == "condition":
        obs_tmp = adata_all.obs[[sample_key, condition_key]].copy()
        obs_tmp[sample_key] = obs_tmp[sample_key].astype(str)
        obs_tmp[condition_key] = obs_tmp[condition_key].astype(str)
        cond_per_sample = (
            obs_tmp
            .assign(_i=np.arange(len(obs_tmp)))
            .sort_values("_i")
            .groupby(sample_key, sort=False)[condition_key]
            .first()
        )
        ordered_samples = (
            cond_per_sample.reset_index()
            .sort_values([condition_key, sample_key], kind="stable")[sample_key]
            .tolist()
        )
    else:
        raise ValueError("sample_order_mode must be one of {'original','sample','condition'}")

    n_samples = len(ordered_samples)
    total_positions = n_rows * n_cols
    if n_samples > total_positions:
        raise ValueError(
            f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {n_samples} samples"
        )

    sample_masks = {
        sample_name: (sample_series == sample_name).to_numpy()
        for sample_name in ordered_samples
    }

    x_spatial = adata_all.obsm[spatial_key]
    x_new = np.empty_like(x_spatial, dtype=float)
    grid_bounds = {}
    sample_positions = {}

    for idx, sample_name in enumerate(ordered_samples):
        row = idx // n_cols
        col = idx % n_cols
        sample_positions[sample_name] = (row, col)

        xspa = x_spatial[sample_masks[sample_name]]
        l_max = xspa.max(axis=0)
        l_min = xspa.min(axis=0)

        grid_bounds[(row, col)] = {
            "width": float(l_max[0] - l_min[0]),
            "height": float(l_max[1] - l_min[1]),
            "min_x": float(l_min[0]),
            "min_y": float(l_min[1]),
            "max_x": float(l_max[0]),
            "max_y": float(l_max[1]),
        }

    row_heights = [
        max(
            (grid_bounds[(r, c)]["height"] for c in range(n_cols) if (r, c) in grid_bounds),
            default=0.0,
        )
        for r in range(n_rows)
    ]
    col_widths = [
        max(
            (grid_bounds[(r, c)]["width"] for r in range(n_rows) if (r, c) in grid_bounds),
            default=0.0,
        )
        for c in range(n_cols)
    ]

    row_y_offsets = [0.0]
    for i in range(n_rows - 1):
        row_y_offsets.append(row_y_offsets[-1] - row_heights[i] - tile_spacing)

    col_x_offsets = [0.0]
    for i in range(n_cols - 1):
        col_x_offsets.append(col_x_offsets[-1] + col_widths[i] + tile_spacing)

    for sample_name in ordered_samples:
        row, col = sample_positions[sample_name]
        mask = sample_masks[sample_name]
        xspa = x_spatial[mask].copy().astype(float)

        bounds = grid_bounds[(row, col)]
        target_x = col_x_offsets[col]
        target_y = row_y_offsets[row]

        if flipy:
            center_y = (bounds["min_y"] + bounds["max_y"]) / 2.0
            xspa[:, 1] = 2.0 * center_y - xspa[:, 1]

        dif_x = target_x - bounds["min_x"]
        dif_y = target_y - (bounds["max_y"] if not flipy else bounds["min_y"])

        xspa[:, 0] += dif_x
        xspa[:, 1] += dif_y
        x_new[mask] = xspa

    adata_all.obsm[new_obsm_key] = x_new


def rename_obs_keys(adata: AnnData) -> AnnData:
    """Add standardized obs columns expected by the remaining tabs."""
    key_map = {
        "Sample": "sample",
        "nFrags": "n_fragment",
        "Condition": "condition",
        "Clusters": "cluster",
    }

    keys = adata.obs_keys()
    for src, dest in key_map.items():
        if src in keys and dest not in keys:
            adata.obs[dest] = adata.obs[src]

    return adata


def reorder_obs_columns(adata: AnnData, first_col: str = "cluster") -> None:
    """Move a preferred obs column to the front when present."""
    if first_col in adata.obs.columns:
        new_order = [first_col] + [c for c in adata.obs.columns if c != first_col]
        adata.obs = adata.obs[new_order]


def sort_group_categories(values):
    """Sort group labels numerically if possible, else alphabetically."""
    num_vals = []
    for v in values:
        try:
            num_vals.append(float(v))
        except (ValueError, TypeError):
            return sorted(map(str, values))

    sorted_pairs = sorted(zip(values, num_vals), key=lambda x: x[1])
    return [v for v, _ in sorted_pairs]


def sync_obs_metadata(
    adata1: AnnData,
    adata2: AnnData,
    *,
    reconcile_shared: bool = True,
    reconcile_direction: str = "adata1_to_adata2",
) -> None:
    """Synchronize non-numeric obs metadata between two AnnData objects."""
    idx1 = pd.Index(adata1.obs_names)
    idx2 = pd.Index(adata2.obs_names)

    if not idx1.equals(idx2) and set(idx1) != set(idx2):
        raise ValueError("AnnData objects must contain exactly the same cells (obs_names).")

    if reconcile_direction not in {"adata1_to_adata2", "adata2_to_adata1"}:
        raise ValueError("reconcile_direction must be 'adata1_to_adata2' or 'adata2_to_adata1'.")

    obs1 = adata1.obs
    obs2 = adata2.obs

    nonnum1 = [c for c in obs1.columns if not pd.api.types.is_numeric_dtype(obs1[c])]
    nonnum2 = [c for c in obs2.columns if not pd.api.types.is_numeric_dtype(obs2[c])]

    for col in [c for c in nonnum1 if c not in obs2.columns]:
        adata2.obs[col] = obs1[col].reindex(adata2.obs_names)

    for col in [c for c in nonnum2 if c not in obs1.columns]:
        adata1.obs[col] = obs2[col].reindex(adata1.obs_names)

    if not reconcile_shared:
        return

    common_idx = pd.Index(adata1.obs_names)
    shared_nonnum = [c for c in nonnum1 if c in nonnum2]
    for col in shared_nonnum:
        s1 = obs1[col].reindex(common_idx)
        s2 = obs2[col].reindex(common_idx)

        try:
            differs = not s1.equals(s2)
        except Exception:
            differs = True

        if not differs:
            continue

        if reconcile_direction == "adata1_to_adata2":
            adata2.obs[col] = obs1[col].reindex(adata2.obs_names)
        else:
            adata1.obs[col] = obs2[col].reindex(adata1.obs_names)
