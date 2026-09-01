new_data_signal()

# Only render this radial "spoke" view when the display toggle defined in the
# Neighborhood Analysis cell above is set to "radial"; otherwise render nothing
# (the heatmap cell handles the "heatmap" case). Reading `.value` subscribes
# this cell so it reruns when the toggle changes.
_neigh_display_val = neigh_display.value if "neigh_display" in globals() else "heatmap"
if _neigh_display_val != "radial":
    exit()

w_text_output(content="""
## Radial Neighborhood Enrichment

An alternative "spoke" view of the same spatial neighborhood enrichment shown in the heatmap above. Each subplot is one **focal cluster**; spokes point to every other cluster arranged around the circle. Every drawn spoke extends to the circle's **perimeter**, **color** encodes the sign of the z-score (**red** = positive/co-localized, **blue** = negative/avoided), and **thickness** scales with the magnitude of the z-score (thicker = stronger). Hover a spoke to read its z-score. Use the **z-score threshold** to control which connections are drawn.
""")

if not adata_g:
    w_text_output(
        content="No gene activity data selected...",
        appearance={"message_box": "warning"}
    )
    exit()

if "filtered_groups" not in globals():
    filtered_groups = {}


def load_precomputed_neighborhood_groups_radial(adata):
    """Rebuild lightweight AnnData objects from workflow-precomputed results."""
    root = adata.uns.get("cluster_nhood_enrichment_by_group")
    if not isinstance(root, dict) or root.get("schema_version") != 1:
        return {}
    result = {}
    for group_entry in root.get("groups", {}).values():
        group_key = str(group_entry["group_key"])
        subgroups = {}
        for subgroup_entry in group_entry.get("subgroups", {}).values():
            group_value = str(subgroup_entry["group_value"])
            categories = pd.Index(
                np.asarray(subgroup_entry["cluster_categories"]).astype(str)
            )
            obs = pd.DataFrame({
                "cluster": pd.Categorical(
                    categories, categories=categories, ordered=True
                )
            })
            nadata = anndata.AnnData(obs=obs)
            nadata.uns["cluster_nhood_enrichment"] = {
                "zscore": np.asarray(subgroup_entry["zscore"]),
                "count": np.asarray(subgroup_entry["count"]),
            }
            subgroups[group_value] = nadata
        result[group_key] = subgroups
    return result


def make_lightweight_neighborhood_adata_radial(adata, group, subgroup):
    """Subset only metadata and coordinates, never the backed feature matrix."""
    if "spatial_offset" not in adata.obsm:
        raise KeyError(
            "Expected `spatial_offset` in adata.obsm (created in Select Data)."
        )
    mask = (adata.obs[group] == subgroup).to_numpy()
    obs_keys = ["cluster"]
    if "sample" in adata.obs.columns:
        obs_keys.append("sample")
    subset_obs = adata.obs.loc[mask, obs_keys].copy()
    for obs_key in obs_keys:
        if pd.api.types.is_categorical_dtype(subset_obs[obs_key]):
            subset_obs[obs_key] = subset_obs[obs_key].cat.remove_unused_categories()
    lightweight = anndata.AnnData(obs=subset_obs)
    lightweight.obsm["spatial_offset"] = np.asarray(
        adata.obsm["spatial_offset"]
    )[mask].copy()
    return lightweight


def ordered_enrichment_matrix(adata_src, uns_key, mode):
    """Return (zscore_matrix, mode_matrix, numerically-ordered labels)."""
    cat_series = adata_src.obs["cluster"]
    if not pd.api.types.is_categorical_dtype(cat_series):
        cat_series = cat_series.astype("category")
    categories = list(cat_series.cat.categories)
    ordered = sort_group_categories([str(c) for c in categories])
    str_cats = [str(c) for c in categories]
    order_idx = [str_cats.index(c) for c in ordered]

    z = np.asarray(adata_src.uns[uns_key]["zscore"], dtype=float)
    z = z[np.ix_(order_idx, order_idx)]
    m = np.asarray(adata_src.uns[uns_key][mode], dtype=float)
    m = m[np.ix_(order_idx, order_idx)]
    return z, m, ordered


def build_radial_neighborhood_fig(
    z_mat, mode_mat, labels, title, thr=0.0, show_self=False, rmax=None,
    zmax=None,
):
    """Grid of polar 'spoke' plots, one per focal cluster.

    Every drawn spoke reaches the perimeter (fixed length); the z-score
    threshold controls which connections are shown, color encodes sign, and
    spoke thickness scales with the magnitude of the z-score.
    """
    n = len(labels)

    # Fixed radius so every spoke reaches the circle's edge.
    if rmax is None or not np.isfinite(rmax) or rmax <= 0:
        rmax = 1.0

    # Magnitude scale for spoke thickness: normalize |z| against the largest
    # off-diagonal |z| (shared across subgroups when zmax is passed in).
    if zmax is None:
        off_abs = np.abs(z_mat).copy()
        np.fill_diagonal(off_abs, 0.0)
        finite_abs = off_abs[np.isfinite(off_abs)]
        zmax = float(np.nanmax(finite_abs)) if finite_abs.size else 1.0
    if not np.isfinite(zmax) or zmax <= 0:
        zmax = 1.0

    # Angular width range (category units; a full sector is ~1.0).
    w_min, w_max = 0.04, 0.55

    ncols = min(4, n)
    nrows = math.ceil(n / ncols)
    specs = [[{"type": "polar"} for _ in range(ncols)] for _ in range(nrows)]
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        specs=specs,
        subplot_titles=[f"C{lab}" for lab in labels],
        horizontal_spacing=0.04,
        vertical_spacing=0.11,
    )

    for i, focal in enumerate(labels):
        r_row = i // ncols + 1
        c_col = i % ncols + 1
        rs, thetas, colors, widths, zvals = [], [], [], [], []
        for j, other in enumerate(labels):
            if (not show_self) and (other == focal):
                continue
            z_val = z_mat[i, j]
            # Hide connections below the threshold; draw the rest full-length.
            if abs(z_val) < thr:
                continue
            rs.append(rmax)
            thetas.append(str(other))
            colors.append("#C33530" if z_val >= 0 else "#282E66")
            # Thickness encodes |z|: linearly map [thr, zmax] -> [w_min, w_max].
            denom = (zmax - thr) if zmax > thr else 1.0
            frac = (abs(z_val) - thr) / denom
            frac = max(0.0, min(1.0, frac))
            widths.append(w_min + frac * (w_max - w_min))
            zvals.append(float(z_val))
        fig.add_trace(
            go.Barpolar(
                theta=thetas,
                r=rs,
                marker=dict(color=colors, line=dict(width=0)),
                width=widths,
                customdata=zvals,
                showlegend=False,
                hovertemplate=(
                    f"focal C{focal}"
                    "<br>neighbor C%{theta}"
                    "<br>z-score %{customdata:.2f}<extra></extra>"
                ),
            ),
            row=r_row,
            col=c_col,
        )

    fig.update_polars(
        radialaxis=dict(
            range=[0, rmax], showticklabels=False, ticks="", showline=False
        ),
        angularaxis=dict(
            direction="clockwise",
            rotation=90,
            categoryorder="array",
            categoryarray=[str(lab) for lab in labels],
            tickfont=dict(size=7),
        ),
        bgcolor="white",
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=16)),
        width=250 * ncols,
        height=270 * nrows,
        margin=dict(l=20, r=20, t=80, b=20),
    )
    # Lift each subplot's cluster label a few pixels off the plot top for space.
    for ann in fig.layout.annotations:
        ann.font.size = 12
        ann.yshift = 10
    return fig


radial_precomputed = load_precomputed_neighborhood_groups_radial(adata_g)

radial_group_dict = {
    key: adata_g.obs[key].dropna().unique()
    for key in list(adata_g.obs.columns)
    if key != "cluster" and key not in na_keys and (
        pd.api.types.is_object_dtype(adata_g.obs[key])
        or pd.api.types.is_categorical_dtype(adata_g.obs[key])
    )
}

radial_group_by = w_select(
    label="subplot groups",
    key="radial_group_by",
    default="all",
    options=tuple(["all"] + list(radial_group_dict.keys())),
    appearance={
        "detail": "(all, categorical observation)",
        "help_text": "Facet radial plots by any categorical observation."
    }
)

radial_mode = w_select(
    label="value metric",
    key="radial_mode",
    default="zscore",
    options=("zscore", "count"),
    appearance={
        "help_text": "Metric reported in the values table; spoke color/thickness always follow the z-score."
    }
)

radial_thr = w_text_input(
    label="z-score threshold",
    key="radial_thr",
    default="0",
    appearance={
        "help_text": "Only connections with absolute z-score at or above this value are drawn."
    }
)

radial_show_self = w_checkbox(
    label="Show self-enrichment spoke",
    key="radial_show_self",
    default=False,
    appearance={"description": "Include the focal cluster's own diagonal value."}
)

w_row(items=[radial_group_by, radial_mode, radial_thr, radial_show_self])

if (
    radial_group_by.value not in (None, "all")
    and radial_group_by.value not in radial_precomputed
):
    w_text_output(
        content=(
            "This custom annotation was added after the workflow ran. Its spatial "
            "neighborhoods will be computed when first displayed and may take longer."
        ),
        appearance={"message_box": "warning"}
    )

radial_button = w_button(label="Update Radial Plots", key="radial_button")

if radial_group_by.value is not None and radial_button.value:
    try:
        thr_val = float(radial_thr.value) if radial_thr.value not in (None, "") else 0.0
    except (TypeError, ValueError):
        w_text_output(
            content="z-score threshold must be numeric; defaulting to 0.",
            appearance={"message_box": "warning"}
        )
        thr_val = 0.0
    if thr_val < 0:
        thr_val = 0.0
    show_self_val = bool(radial_show_self.value)
    mode_val = radial_mode.value

    radial_figs = []
    radial_table_rows = []

    if radial_group_by.value == "all":
        sample_key = "sample" if "sample" in groups else None
        if "cluster_nhood_enrichment" not in adata_g.uns:
            w_text_output(
                content="Computing neighborhoods for all cells...",
                appearance={"message_box": "info"}
            )
            submit_widget_state()
            squidpy_analysis(adata_g, sample_key=sample_key)
        else:
            w_text_output(
                content="Using existing neighborhood enrichment for all cells...",
                appearance={"message_box": "info"}
            )
            submit_widget_state()

        z_mat, mode_mat, labels = ordered_enrichment_matrix(
            adata_g, "cluster_nhood_enrichment", mode_val
        )
        radial_figs.append(
            build_radial_neighborhood_fig(
                z_mat, mode_mat, labels,
                "All cells: Radial Neighborhood Enrichment",
                thr=thr_val, show_self=show_self_val,
            )
        )
        for i, focal in enumerate(labels):
            for j, other in enumerate(labels):
                radial_table_rows.append({
                    "subgroup": "all",
                    "focal_cluster": focal,
                    "neighbor_cluster": other,
                    "zscore": z_mat[i, j],
                    mode_val: mode_mat[i, j],
                })

    else:
        group = radial_group_by.value
        sub_groups = radial_group_dict[group]
        if group in radial_precomputed:
            filtered_groups[group] = radial_precomputed[group]
        elif group not in filtered_groups:
            filtered_groups[group] = {}
        group_adatas = filtered_groups[group]

        for sg in sub_groups:
            if sg not in group_adatas:
                group_adatas[sg] = make_lightweight_neighborhood_adata_radial(
                    adata_g, group, sg
                )
            if "cluster_nhood_enrichment" not in group_adatas[sg].uns:
                w_text_output(
                    content=f"Computing spatial neighborhoods for {sg}...",
                    appearance={"message_box": "info"}
                )
                submit_widget_state()
                sk = "sample" if "sample" in group_adatas[sg].obs else None
                squidpy_analysis(group_adatas[sg], sample_key=sk)

        # Shared thickness scale across all subgroups for fair comparison.
        matrices = {}
        shared_zmax = 0.0
        for sg in sub_groups:
            z_mat, mode_mat, labels = ordered_enrichment_matrix(
                group_adatas[sg], "cluster_nhood_enrichment", mode_val
            )
            matrices[str(sg)] = (z_mat, mode_mat, labels)
            off_abs = np.abs(z_mat).copy()
            np.fill_diagonal(off_abs, 0.0)
            finite_abs = off_abs[np.isfinite(off_abs)]
            if finite_abs.size:
                shared_zmax = max(shared_zmax, float(np.nanmax(finite_abs)))
        if shared_zmax <= 0:
            shared_zmax = 1.0

        for sg in sort_group_categories([str(s) for s in sub_groups]):
            z_mat, mode_mat, labels = matrices[sg]
            radial_figs.append(
                build_radial_neighborhood_fig(
                    z_mat, mode_mat, labels,
                    f"{group} = {sg}: Radial Neighborhood Enrichment",
                    thr=thr_val, show_self=show_self_val, zmax=shared_zmax,
                )
            )
            for i, focal in enumerate(labels):
                for j, other in enumerate(labels):
                    radial_table_rows.append({
                        "subgroup": sg,
                        "focal_cluster": focal,
                        "neighbor_cluster": other,
                        "zscore": z_mat[i, j],
                        mode_val: mode_mat[i, j],
                    })

    # Render figures by list index so no single global variable is shared
    # between multiple w_plot widgets (which would trip the duplicate-source guard).
    for _idx in range(len(radial_figs)):
        w_plot(source=radial_figs[_idx], key=f"radial_plot_{_idx}")

    radial_enrichment_df = pd.DataFrame(radial_table_rows)
    w_table(
        label="Radial neighborhood enrichment values",
        source=radial_enrichment_df,
        key="radial_enrichment_table",
    )
