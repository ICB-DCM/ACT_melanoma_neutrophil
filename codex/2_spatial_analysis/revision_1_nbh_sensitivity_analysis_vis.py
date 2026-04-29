#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :cd8_mel_codex_analysis -> revision_1_nbh_sensitivity_analysis_vis.py
# @Author : Gemma van der Voort
# @Time   : 27.03.26 16:44
# @Desc   : uses the nbh smoothing diagnostics files generated in revision_0_nbh_based_regions_senstivity_analysis.py
# analyses amount of iterations needed and the stability of the pixels.
# idea: even if there are still changes, if these changes are small cycling reassignments, we can justify stopping prior
# to full convergence.
# '''=================================================
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import os
import ast
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import (processed_data_dir, fig_dir, modifier, modifier_base, parent_dir,
                              processed_spatial_data_dir, data_repo_path)
from utils_codex import map_cell_type, spatial_metacluster_filtering
# import stylesheet
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style_old.mplstyle"))
fig_dir_full = f'{fig_dir}/manuscipt_plots/revision'
img_modifier = 'v4_corrected'  # for the data corrections done at the phenotyping level.
version_key = 'v4_8'  # version of the metacluster assignment
metacluster_key_orig = f'Metacluster {version_key}'
modifier = ''
max_iterations = [5, 10, 11, 12, 13, 14, 15, 20, 25, 30, 50]
#%% load data
dict_diagnostics = {}
for max_iteration in max_iterations:
    metacluster_col_name = f"{metacluster_key_orig}_filtered_0.5_max_it_{max_iteration}_{version_key}"
    df = pd.read_csv(f"{processed_spatial_data_dir}/revision/metacluster_filtering_diagnostics_{metacluster_col_name}.csv")
    # fix list to string conversion issue that occurred when saving:
    cols = ['changed_cell_ids', 'old_labels', 'new_labels', 'transition_counts']
    for col in cols:
        df[col] = df[col].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else x
        )
    dict_diagnostics[max_iteration] = df
#%% load adata for total cell amounts per image
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{version_key}.h5ad")
df_10 = dict_diagnostics[10]
# df_5 = dict_diagnostics[5]
# df_11 = dict_diagnostics[11]
# df_12 = dict_diagnostics[12]
# df_13 = dict_diagnostics[13]
# df_14 = dict_diagnostics[14]
# df_15 = dict_diagnostics[15]
# df_20 = dict_diagnostics[20]
# df_30 = dict_diagnostics[30]
df_50 = dict_diagnostics[50]
#%% check for cycling reassignments in the final iterations.

def ensure_deserialized_diagnostics(
    diagnostics_df,
    cols=('changed_cell_ids', 'old_labels', 'new_labels', 'transition_counts')
):
    diagnostics_df = diagnostics_df.copy()

    for col in cols:
        if col not in diagnostics_df.columns:
            continue

        diagnostics_df[col] = diagnostics_df[col].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else x
        )

    return diagnostics_df


def explode_spatial_filter_diagnostics(diagnostics_df):
    records = []

    for _, row in diagnostics_df.iterrows():
        cell_ids = row['changed_cell_ids']
        old_labels = row['old_labels']
        new_labels = row['new_labels']

        if not isinstance(cell_ids, list) or len(cell_ids) == 0:
            continue

        if not (len(cell_ids) == len(old_labels) == len(new_labels)):
            raise ValueError(
                f"Mismatched diagnostics lengths for image={row['image']}, "
                f"iteration={row['iteration']}: "
                f"{len(cell_ids)=}, {len(old_labels)=}, {len(new_labels)=}"
            )

        for cell_id, old_label, new_label in zip(cell_ids, old_labels, new_labels):
            records.append({
                'image': row['image'],
                'iteration': row['iteration'],
                'cell_id': cell_id,
                'old_label': old_label,
                'new_label': new_label,
            })

    return pd.DataFrame(records)


def compute_bounce_events_by_iteration(diagnostics_df):
    """
    Return one row per bouncing cell event at iteration t.

    A bounce at iteration t is:
        t-1: A -> B
        t:   B -> A
    """
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)
    diag_long = explode_spatial_filter_diagnostics(diagnostics_df)

    if diag_long.empty:
        return pd.DataFrame(columns=[
            'image',
            'iteration',
            'cell_id',
            'old_label_prev',
            'new_label_prev',
            'old_label',
            'new_label',
            'bounce_pair'
        ])

    bounce_records = []

    for image, df_img in diag_long.groupby('image'):
        for cell_id, df_cell in df_img.groupby('cell_id'):
            df_cell = df_cell.sort_values('iteration').copy()

            if len(df_cell) < 2:
                continue

            prev_row = None
            for _, row in df_cell.iterrows():
                if prev_row is not None:
                    consecutive = row['iteration'] == prev_row['iteration'] + 1
                    is_reverse = (
                        prev_row['old_label'] == row['new_label']
                        and prev_row['new_label'] == row['old_label']
                        and prev_row['old_label'] != prev_row['new_label']
                    )

                    if consecutive and is_reverse:
                        label_a = prev_row['old_label']
                        label_b = prev_row['new_label']
                        bounce_pair = tuple(sorted([str(label_a), str(label_b)]))

                        bounce_records.append({
                            'image': image,
                            'iteration': row['iteration'],
                            'cell_id': cell_id,
                            'old_label_prev': prev_row['old_label'],
                            'new_label_prev': prev_row['new_label'],
                            'old_label': row['old_label'],
                            'new_label': row['new_label'],
                            'bounce_pair': bounce_pair
                        })

                prev_row = row

    return pd.DataFrame(bounce_records)

def compute_bouncing_by_iteration(diagnostics_df):
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)
    bounce_events_df = compute_bounce_events_by_iteration(diagnostics_df)

    per_iter = diagnostics_df[['image', 'iteration', 'n_changes']].copy()
    per_iter['n_bouncing_cells'] = 0

    if not bounce_events_df.empty:
        bounce_counts = (
            bounce_events_df.groupby(['image', 'iteration'])['cell_id']
            .nunique()
            .reset_index(name='n_bouncing_cells')
        )

        per_iter = per_iter.drop(columns=['n_bouncing_cells']).merge(
            bounce_counts,
            on=['image', 'iteration'],
            how='left'
        )
        per_iter['n_bouncing_cells'] = per_iter['n_bouncing_cells'].fillna(0).astype(int)

    per_iter['pct_bouncing_of_total'] = np.where(
        per_iter['n_changes'] > 0,
        100 * per_iter['n_bouncing_cells'] / per_iter['n_changes'],
        0.0
    )

    return per_iter.sort_values(['image', 'iteration']).reset_index(drop=True)


def get_first_all_bouncing_iteration_per_image(bouncing_df):
    """
    First iteration per image where all remaining changes are bouncing:
        n_changes > 0 and pct_bouncing_of_total == 100
    """
    df = bouncing_df.copy()

    first_hits = (
        df[(df['n_changes'] > 0) & (df['pct_bouncing_of_total'] == 100)]
        .sort_values(['image', 'iteration'])
        .groupby('image', as_index=False)
        .first()[['image', 'iteration']]
        .rename(columns={'iteration': 'first_all_bouncing_iteration'})
    )

    return first_hits


def get_last_image_all_bouncing_iteration(bouncing_df):
    """
    Iteration where the last image first reaches the all-bouncing regime.
    """
    first_hits = get_first_all_bouncing_iteration_per_image(bouncing_df)

    if first_hits.empty:
        return None, first_hits

    last_iteration = first_hits['first_all_bouncing_iteration'].max()

    return last_iteration, first_hits

def summarize_metric_by_iteration(df, value_col):
    summary = (
        df.groupby('iteration')[value_col]
        .agg(['mean', 'std', 'count'])
        .reset_index()
        .rename(columns={'count': 'n'})
    )

    summary['std'] = summary['std'].fillna(0.0)
    summary['sem'] = np.where(summary['n'] > 1, summary['std'] / np.sqrt(summary['n']), 0.0)
    summary['ci_low'] = summary['mean'] - 1.96 * summary['sem']
    summary['ci_high'] = summary['mean'] + 1.96 * summary['sem']

    return summary

def plot_pct_bouncing_per_iteration(
    diagnostics_df,
    plot_mode='individual_and_mean',
    figsize=(8, 6),
    annotate_last_all_bouncing=True
):
    valid_modes = ['individual_and_mean', 'mean_and_ci']
    if plot_mode not in valid_modes:
        raise ValueError(f"plot_mode must be one of {valid_modes}, got {plot_mode}")

    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)
    bouncing_df = compute_bouncing_by_iteration(diagnostics_df)
    summary = summarize_metric_by_iteration(bouncing_df, value_col='pct_bouncing_of_total')

    fig, ax = plt.subplots(figsize=figsize)

    if plot_mode == 'individual_and_mean':
        for image, df_img in bouncing_df.groupby('image'):
            df_img = df_img.sort_values('iteration')
            ax.plot(
                df_img['iteration'],
                df_img['pct_bouncing_of_total'],
                alpha=0.25,
                linewidth=1.0
            )

        ax.plot(
            summary['iteration'],
            summary['mean'],
            linewidth=2.5,
            label='mean'
        )
        ax.legend(frameon=False)

    elif plot_mode == 'mean_and_ci':
        ax.plot(
            summary['iteration'],
            summary['mean'],
            linewidth=2.5,
            label='mean'
        )
        ax.fill_between(
            summary['iteration'],
            summary['ci_low'],
            summary['ci_high'],
            alpha=0.25,
            label='95% confidence interval'
        )
        ax.legend(frameon=False)

    if annotate_last_all_bouncing:
        last_iter, first_hits = get_last_image_all_bouncing_iteration(bouncing_df)

        if last_iter is not None:
            ax.axvline(
                x=last_iter,
                linestyle='--',
                linewidth=1.5
            )

            ymax = ax.get_ylim()[1]
            x_offset = 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0])
            y_offset = 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])

            ax.text(
                last_iter + x_offset,
                0.95 * ymax - y_offset,
                f"{last_iter}",
                va='top',
                ha='left'
            )
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Oscillating changes (% of total changes)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    return fig
#%%
def plot_bounce_transition_heatmap_at_stability(
    diagnostics_df,
    figsize=(8, 7),
    normalize=False
):
    """
    Plot a transition matrix of bounce pairs at the first all-bouncing iteration
    for each image.

    Parameters
    ----------
    diagnostics_df : pd.DataFrame
        Diagnostics dataframe from spatial_metacluster_filtering().
    figsize : tuple
        Figure size.
    normalize : bool
        If True, divide counts by total number of bounce events included.

    Returns
    -------
    fig : matplotlib.figure.Figure
    transition_matrix : pd.DataFrame
    selected_bounce_events_df : pd.DataFrame
    """
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)
    bouncing_df = compute_bouncing_by_iteration(diagnostics_df)
    bounce_events_df = compute_bounce_events_by_iteration(diagnostics_df)
    first_hits = get_first_all_bouncing_iteration_per_image(bouncing_df)

    label_map = {
        'T-cell zone': 'T cell zone',
        'B-cell follicle': 'B cell follicle',
        'Medulla-Interfollicular zone-SCS': 'Medulla/SCS\nInterfollicular zone',
    }

    if first_hits.empty or bounce_events_df.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'no all-bouncing regime detected', ha='center', va='center')
        ax.axis('off')
        return fig, pd.DataFrame(), pd.DataFrame()

    selected_bounce_events_df = bounce_events_df.merge(
        first_hits,
        left_on=['image', 'iteration'],
        right_on=['image', 'first_all_bouncing_iteration'],
        how='inner'
    )

    if selected_bounce_events_df.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'no bounce events found at all-bouncing iteration', ha='center', va='center')
        ax.axis('off')
        return fig, pd.DataFrame(), pd.DataFrame()

    # use unordered pairs because A<->B and B<->A are the same bounce relationship
    selected_bounce_events_df[['label_1', 'label_2']] = pd.DataFrame(
        selected_bounce_events_df['bounce_pair'].tolist(),
        index=selected_bounce_events_df.index
    )

    transition_matrix = pd.crosstab(
        selected_bounce_events_df['label_1'],
        selected_bounce_events_df['label_2']
    )

    # rename labels
    transition_matrix.index = transition_matrix.index.map(label_map)
    transition_matrix.columns = transition_matrix.columns.map(label_map)

    all_labels = sorted(set(transition_matrix.index).union(set(transition_matrix.columns)))
    transition_matrix = transition_matrix.reindex(
        index=all_labels,
        columns=all_labels,
        fill_value=0
    )

    # make symmetric for display
    transition_matrix = transition_matrix + transition_matrix.T
    np.fill_diagonal(transition_matrix.values, 0)

    if normalize:
        total = transition_matrix.values.sum()
        if total > 0:
            transition_matrix = (transition_matrix / total) * 100

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(transition_matrix.values, aspect='auto')

    ax.set_xticks(np.arange(len(transition_matrix.columns)))
    ax.set_xticklabels(transition_matrix.columns, rotation=90)
    ax.set_yticks(np.arange(len(transition_matrix.index)))
    ax.set_yticklabels(transition_matrix.index)

    ax.set_xlabel('')
    ax.set_ylabel('')

    cbar = fig.colorbar(im, ax=ax)
    if normalize:
        cbar.set_label('Oscillating changes\n[% of total changes]')
    else:
        cbar.set_label('Oscillating changes')

    plt.tight_layout()

    return fig, transition_matrix, selected_bounce_events_df


def _plot_metric_with_mode(
    df,
    value_col,
    ylabel,
    plot_mode='individual_and_mean',
    alpha_individual=0.25,
    linewidth_individual=1.0,
    linewidth_mean=2.5,
    figsize=(8, 6)
):
    valid_modes = ['individual_and_mean', 'mean_and_ci']
    if plot_mode not in valid_modes:
        raise ValueError(f"plot_mode must be one of {valid_modes}, got {plot_mode}")

    fig, ax = plt.subplots(figsize=figsize)

    summary = summarize_metric_by_iteration(df, value_col=value_col)

    if plot_mode == 'individual_and_mean':
        for image, df_img in df.groupby('image'):
            df_img = df_img.sort_values('iteration')
            ax.plot(
                df_img['iteration'],
                df_img[value_col],
                alpha=alpha_individual,
                linewidth=linewidth_individual
            )

        ax.plot(
            summary['iteration'],
            summary['mean'],
            linewidth=linewidth_mean,
            label='mean'
        )
        ax.legend(frameon=False)

    elif plot_mode == 'mean_and_ci':
        ax.plot(
            summary['iteration'],
            summary['mean'],
            linewidth=linewidth_mean,
            label='mean'
        )
        ax.fill_between(
            summary['iteration'],
            summary['ci_low'],
            summary['ci_high'],
            alpha=0.25,
            label='95% confidence interval'
        )
        ax.legend(frameon=False)

    ax.set_xlabel('iteration')
    ax.set_ylabel(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # sns.despine(ax=ax)
    plt.tight_layout()

    return fig

def plot_total_changes_per_iteration(
    diagnostics_df,
    plot_mode='individual_and_mean',
    figsize=(8, 6),
    relative=False,
    n_cells_per_image=None
):
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)

    df_plot = diagnostics_df[['image', 'iteration', 'n_changes']].copy()

    if relative:
        if n_cells_per_image is None:
            raise ValueError("n_cells_per_image must be provided when relative=True")

        df_plot['n_cells'] = df_plot['image'].map(n_cells_per_image)

        if df_plot['n_cells'].isna().any():
            missing = df_plot.loc[df_plot['n_cells'].isna(), 'image'].unique()
            raise ValueError(f"Missing cell counts for images: {missing}")

        df_plot['value'] = 100 * df_plot['n_changes'] / df_plot['n_cells']
        ylabel = 'changes (% of cells)'

    else:
        df_plot['value'] = df_plot['n_changes']
        ylabel = 'number of changes'

    fig = _plot_metric_with_mode(
        df=df_plot.rename(columns={'value': 'metric'}),
        value_col='metric',
        ylabel=ylabel,
        plot_mode=plot_mode,
        figsize=figsize
    )

    return fig



# bounce_summary_df5, bouncing_events_df5 = quantify_bouncing_in_final_two_iterations(df_5)
# bounce_summary_df10, bouncing_events_df10 = quantify_bouncing_in_final_two_iterations(df_10)
# bounce_summary_df11, bouncing_events_df11 = quantify_bouncing_in_final_two_iterations(df_11)
# bounce_summary_df12, bouncing_events_df12 = quantify_bouncing_in_final_two_iterations(df_12)
# bounce_summary_df13, bouncing_events_df13 = quantify_bouncing_in_final_two_iterations(df_13)
# bounce_summary_df14, bouncing_events_df14 = quantify_bouncing_in_final_two_iterations(df_14)
# bounce_summary_df15, bouncing_events_df15 = quantify_bouncing_in_final_two_iterations(df_15)
# bounce_summary_df20, bouncing_events_df20 = quantify_bouncing_in_final_two_iterations(df_20)
# none in 5, none in 10.
#%%
# bounce_summary_df['pct_bouncing_events_of_total_changes'] = (
#     100 * (2 * bounce_summary_df['n_bouncing_cells'])
#     / bounce_summary_df['total_changes_final_two_iterations']
# )
# explore the df
diag_long = explode_spatial_filter_diagnostics(df_10)

# how often each cell changes
cell_change_counts = (
    diag_long.groupby(['image', 'cell_id'])
    .size()
    .reset_index(name='n_changes')
)

# focus on unstable cells
unstable_cells = cell_change_counts[cell_change_counts['n_changes'] >= 3]

print(unstable_cells.sort_values('n_changes', ascending=False).head())
#%% plot
fig1 = plot_total_changes_per_iteration(
    diagnostics_df=df_50,
    figsize=(9, 7),
    plot_mode='individual_and_mean')
fig1.savefig(f"{fig_dir_full}/spatial_filtering_total_changes_per_iteration_individual_and_mean_{modifier}.pdf", dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2 = plot_pct_bouncing_per_iteration(
    diagnostics_df=df_50,
    plot_mode='individual_and_mean',
    figsize=(9, 7),
    annotate_last_all_bouncing=True)
fig2.savefig(f"{fig_dir_full}/spatial_filtering_pct_bouncing_per_iteration_individual_and_mean_{modifier}.pdf", dpi=300, bbox_inches='tight')
plt.close(fig2)

fig3 = plot_total_changes_per_iteration(
    diagnostics_df=df_50,
    figsize=(9, 7),
    plot_mode='mean_and_ci')
fig3.savefig(f"{fig_dir_full}/spatial_filtering_total_changes_per_iteration_mean_and_ci_{modifier}.pdf", dpi=300, bbox_inches='tight')
plt.close(fig3)

fig4 = plot_pct_bouncing_per_iteration(
    diagnostics_df=df_50,
    plot_mode='mean_and_ci',
    figsize=(9, 7),
    annotate_last_all_bouncing=True)
fig4.savefig(f"{fig_dir_full}/spatial_filtering_pct_bouncing_per_iteration_mean_and_ci_{modifier}.pdf", dpi=300, bbox_inches='tight')
plt.close(fig4)

#%% which transitions heatmap
fig_heatmap, transition_matrix, selected_bounce_events_df = plot_bounce_transition_heatmap_at_stability(
    diagnostics_df=df_50,
    normalize=False
)
fig_heatmap.savefig(f"{fig_dir_full}/spatial_filtering_bounce_transition_heatmap_at_stability_{modifier}.pdf",
                    dpi=300, bbox_inches='tight')
plt.close(fig_heatmap)
fig_heatmap2, transition_matrix, selected_bounce_events_df = plot_bounce_transition_heatmap_at_stability(
    diagnostics_df=df_50,
    normalize=True
)
fig_heatmap2.savefig(f"{fig_dir_full}/spatial_filtering_bounce_transition_heatmap_at_stability_normalised_{modifier}.pdf",
                    dpi=300, bbox_inches='tight')
#%% relvative changes per image
def get_n_cells_per_image(adata, library_key='dataset_name'):
    return adata.obs[library_key].value_counts().to_dict()

n_cells_per_image = get_n_cells_per_image(adata)

fig = plot_total_changes_per_iteration(
    diagnostics_df=df_50,
    plot_mode='individual_and_mean',
    relative=True,
    n_cells_per_image=n_cells_per_image
)
fig.savefig(f"{fig_dir_full}/spatial_filtering_total_changes_per_iteration_individual_and_mean_relative_{modifier}.pdf",
            dpi=300, bbox_inches='tight')
plt.close(fig)
#%% percentage of cells still makign a useful change at iteration 10
def pct_non_oscillating_changes_at_iteration(
    diagnostics_df,
    iteration_of_interest=10
):
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)

    # long format
    diag_long = explode_spatial_filter_diagnostics(diagnostics_df)

    if diag_long.empty:
        return 0.0

    # --- 1. all cells that ever changed
    all_changed_cells = set(diag_long['cell_id'].unique())

    # --- 2. cells that are part of ANY bouncing event (future oscillators)
    bounce_events_df = compute_bounce_events_by_iteration(diagnostics_df)

    if bounce_events_df.empty:
        bouncing_cells = set()
    else:
        bouncing_cells = set(bounce_events_df['cell_id'].unique())

    # --- 3. cells changing at the iteration of interest
    cells_at_iter = set(
        diag_long.loc[
            diag_long['iteration'] == iteration_of_interest,
            'cell_id'
        ]
    )

    # --- 4. remove oscillators
    useful_cells = cells_at_iter - bouncing_cells

    # --- 5. compute percentage
    if len(all_changed_cells) == 0:
        pct = 0.0
    else:
        pct = 100 * len(useful_cells) / len(all_changed_cells)

    # optional: return diagnostics too
    return {
        'iteration': iteration_of_interest,
        'n_all_changed_cells': len(all_changed_cells),
        'n_cells_at_iteration': len(cells_at_iter),
        'n_bouncing_cells': len(bouncing_cells),
        'n_useful_cells': len(useful_cells),
        'pct_useful_non_oscillating': pct
    }

result = pct_non_oscillating_changes_at_iteration(
    diagnostics_df=df_50,
    iteration_of_interest=10
)

print(result)
result25 = pct_non_oscillating_changes_at_iteration(
    diagnostics_df=df_50,
    iteration_of_interest=25
)

print(result25)
# iterate over iterations
df_results_overview = []
for iteration in [5, 10, 15, 20, 25, 30, 50]:
    result_iter = pct_non_oscillating_changes_at_iteration(
        diagnostics_df=df_50,
        iteration_of_interest=iteration
    )
    df_results_overview.append(result_iter)
df_results_overview = pd.DataFrame(df_results_overview)
# save
df_results_overview.to_csv(f"{processed_spatial_data_dir}/spatial_filtering_smoothing_useful_non_"
                           f"oscillating_changes_over_iterations_{modifier}.csv", index=False)
#%%
def count_images_stable_at_iteration(
    diagnostics_df,
    iteration_of_interest=10,
    include_fully_converged=True
):
    diagnostics_df = ensure_deserialized_diagnostics(diagnostics_df)
    bouncing_df = compute_bouncing_by_iteration(diagnostics_df)

    df_iter = bouncing_df.loc[bouncing_df['iteration'] == iteration_of_interest].copy()

    if include_fully_converged:
        stable_mask = (
            (df_iter['n_changes'] == 0) |
            ((df_iter['n_changes'] > 0) & (df_iter['pct_bouncing_of_total'] == 100))
        )
    else:
        stable_mask = (
            (df_iter['n_changes'] > 0) & (df_iter['pct_bouncing_of_total'] == 100)
        )

    stable_images = sorted(df_iter.loc[stable_mask, 'image'].unique())

    result = {
        'iteration': iteration_of_interest,
        'n_images_total': df_iter['image'].nunique(),
        'n_images_stable': len(stable_images),
        'pct_images_stable': 100 * len(stable_images) / df_iter['image'].nunique()
        if df_iter['image'].nunique() > 0 else np.nan,
        'stable_images': stable_images
    }

    return result, df_iter.loc[stable_mask].copy()

result_10, stable_df_10 = count_images_stable_at_iteration(
    diagnostics_df=df_50,
    iteration_of_interest=10,
    include_fully_converged=True
)

print(result_10)