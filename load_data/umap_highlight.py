def plot_umap_for_samples(
  adata,
  samples,
  color_by='cluster',
  pt_size=3,
  coords="spatial",
  flipY=True,
  color_scheme='bright',
  show_cluster=None
):
    import numpy as np
    import pandas as pd

    # Determine if color_by is discrete or continuous
    obs_values = adata.obs[color_by]
    is_discrete = pd.api.types.is_categorical_dtype(obs_values) or obs_values.dtype.name == 'category' or obs_values.nunique() < 30

    print(f"Coloring by: {color_by} (Discrete: {is_discrete})")

    # For discrete, build color map
    if is_discrete:
        obs_groups = sorted(obs_values.unique())
        colors = generate_color_palette(len(obs_groups), color_scheme)
        group_color_map = {obs_groups[i]: colors[i] for i in range(len(obs_groups))}
    else:
        group_color_map = None

    # Plot grid layout
    num_rows = (len(samples) - 1) // 3 + 1
    num_cols = min(len(samples), 3)

    combined_fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        subplot_titles=samples,
        horizontal_spacing=0,
        vertical_spacing=0.05
    )

    shown_clusters = set() if is_discrete else None

    # Flip Y if needed
    if flipY:
        flipped_y = adata.obsm['spatial'].copy()
        flipped_y[:, 1] = -flipped_y[:, 1]
        adata.obsm['spatial_flippedY'] = flipped_y
        coords = 'spatial_flippedY'

    for i, sample in enumerate(samples):
        row = i // 3 + 1
        col = i % 3 + 1

        sample_data = adata[adata.obs['sample'] == sample]

        if show_cluster is not None:
            # Split data into highlighted (in cluster) and masked (other clusters)
            highlight_mask = sample_data.obs['cluster'] == show_cluster
            highlight_data = sample_data[highlight_mask]
            background_data = sample_data[~highlight_mask]

            # Background layer (greyed out)
            bg_fig = snap.pl.umap(
                background_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in bg_fig.data:
                trace['marker']['color'] = 'lightgrey'
                trace['marker']['opacity'] = 0.3
                trace.showlegend = False
                combined_fig.add_trace(trace, row=row, col=col)

            # Foreground layer (highlighted cluster)
            fg_fig = snap.pl.umap(
                highlight_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fg_fig.data:
                if is_discrete:
                    trace['marker']['color'] = group_color_map.get(trace.name, 'black')
                    if trace.name not in shown_clusters:
                        shown_clusters.add(trace.name)
                    else:
                        trace.showlegend = False
                combined_fig.add_trace(trace, row=row, col=col)

        else:
            # No masking — full sample
            fig = snap.pl.umap(
                sample_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fig.data:
                if is_discrete:
                    trace['marker']['color'] = group_color_map.get(trace.name, 'black')
                    if trace.name not in shown_clusters:
                        shown_clusters.add(trace.name)
                    else:
                        trace.showlegend = False
                combined_fig.add_trace(trace, row=row, col=col)

    # Layout
    subplot_width = 500
    subplot_height = 500
    combined_fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        autosize=False,
        width=subplot_width * num_cols,
        height=subplot_height * num_rows,
        margin=dict(l=0, r=0, t=25, b=0, pad=0),
        legend=dict(
            font=dict(size=18),
            itemsizing='constant',
            itemwidth=30,
            xanchor='right',
            yanchor='top',
            x=1.1
        )
    )

    combined_fig.update_coloraxes(colorscale='Spectral_r', colorbar_title=color_by)

    for i in range(1, num_rows + 1):
        for j in range(1, num_cols + 1):
            combined_fig.update_xaxes(
                scaleratio=1,
                showticklabels=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                row=i,
                col=j
            )
            combined_fig.update_yaxes(
                showticklabels=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                row=i,
                col=j
            )

    return combined_fig
