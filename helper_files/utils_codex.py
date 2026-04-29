#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : frequently used functions
# @Desc updated: Shared CODEx utility functions.
# '''=================================================
import pandas as pd
import numpy as np
import networkx as nx
import re
from scipy import stats
from tqdm import tqdm
import scanpy as sc
from collections import Counter
from sklearn.neighbors import NearestNeighbors

def adata_to_df(adata, additional_columns=None):
    """
    Transforms anndata object to a dataframe with columns 'X', 'Y' and 'name', where name
    contains cell type information.
    """
    df = pd.DataFrame()
    try:
        df['Cell type'] = adata.obs['Cell type']
    except KeyError:  # dataset not fine-tuned yet
        df['Cell type'] = adata.obs['Compartment']
    df['X'] = adata.obsm['position'].T[0].copy()
    df['Y'] = adata.obsm['position'].T[1].copy()
    if additional_columns:
        df = pd.concat([df, adata.obs[additional_columns]], axis=1)
    return df


def prune_csv_to_markers(first_line, markers_of_interest):
    """
    Prunes the first line of a csv file to only contain the markers in the marker list. Then verifies the order of the
    markers in the list.
    :param first_line: str, first line of a csv file
    :param markers_of_interest: list, list of markers to find indices off
    :return: indices of markers of interest in the csv file
    """
    first_line_list = first_line.split(',')  # make the first_line string into a list split by ','
    # keep only the entries ending in 'foreground-mean_membrane'. assert that there are any (could have been omitted in
    # preprocessing)
    assert any('foreground-mean_membrane' in entry for entry in first_line_list), \
        'No foreground-mean_membrane in first line, please pick another measurement.'
    pruned_first_line_list = [entry for entry in first_line_list if 'foreground-mean_membrane' in entry]
    # from each entry, remove '_foreground-mean_membrane' and cycle and channel info by splitting at '_'
    pruned_first_line_list = [entry.split('_')[2] for entry in pruned_first_line_list]
    # now remove the extra DAPI channels. keep only the first one.
    # index of first DAPI channel
    dapi_index = [index for index, entry in enumerate(pruned_first_line_list) if 'DAPI' in entry][0]
    # remove all entries that contain DAPI except the one at dapi_index (most likely 0)
    pruned_first_line_list = [entry for entry in pruned_first_line_list if 'DAPI' not in entry or
                              pruned_first_line_list.index(entry) == dapi_index]
    # remove the number at the end of the DAPI channel name
    pruned_first_line_list[dapi_index] = re.sub(r'\d+', '', pruned_first_line_list[dapi_index])
    # remove all entries that contain 'Blank' or 'blank'
    pruned_first_line_list_6 = [entry for entry in pruned_first_line_list if 'Blank' not in entry and
                                'blank' not in entry]
    # remove all entries that contain 'Empty' or 'empty'
    pruned_first_line_list_6 = [entry for entry in pruned_first_line_list_6 if 'Empty' not in entry and
                                'empty' not in entry]
    # replace any '/' with '-'
    pruned_first_line_list_6 = [entry.replace('/', '-') for entry in pruned_first_line_list_6]
    # now match the pruned list to the marker list. assert that markers of interest are in the csv file.
    indices = [pruned_first_line_list_6.index(marker) for marker in markers_of_interest]
    assert all(marker in pruned_first_line_list_6 for marker in markers_of_interest), \
        'Not all markers of interest are in the csv file.'
    return indices, pruned_first_line_list_6

def map_cell_type(cell_type, mapping):
    for high_level, types in mapping.items():
        if cell_type in types:
            return high_level
    return cell_type  # Return original if no match found


def graph_information(adata_for_graph):
    """
    Check a few important graph properties and return them as a pd df.
    :param adata_for_graph: anndata object with spatial graph saved in adata.obsp['spatial_connectivities'] as a sparse
    matrix according to the squidpy documentation.
    :return:
    """
    # extract a workable graph from adata
    connectivities = adata_for_graph.obsp['spatial_connectivities']
    # Create a NetworkX graph from the connectivity matrix
    graph = nx.from_scipy_sparse_array(connectivities)
    # average amount of edges per node + distribution (quantiles)
    edge_per_node = [len(list(graph.neighbors(node))) for node in graph.nodes]
    avg_edges_neighbours = np.mean(edge_per_node)
    edge_neighbours_quantiles = np.quantile(edge_per_node, [0.25, 0.5, 0.75])
    # average edge length + distribution (quantiles)
    edge_lengths = []
    for edge in graph.edges:
        edge_lengths.append(np.linalg.norm(adata_for_graph.obsm['position'][edge[0]] -
                                           adata_for_graph.obsm['position'][edge[1]]))
    edge_len_avg = np.mean(edge_lengths)
    edge_len_quantiles = np.quantile(edge_lengths, [0.25, 0.5, 0.75])
    # amount of nodes in main graph
    n_nodes = len(graph.nodes)
    # amount of single points --> to be filtered out
    single_points = [node for node in graph.nodes if len(list(graph.neighbors(node))) == 0]
    n_single_points = len(single_points)
    # amount of graphs containing less than 5% of data --> to be filtered out after inspection
    small_graphs = [c for c in nx.connected_components(graph) if len(c) < 0.05 * n_nodes]
    # transform to list of lists
    small_graphs = [[list(c)] for c in small_graphs]
    # count all nodes in the small graphs
    all_nodes_small_graphs = [node for graph in small_graphs for node in graph[0]]
    n_nodes_small_graphs = len(all_nodes_small_graphs)
    # remove the single points from the small graphs. graph is a list of a set, so access element 0
    small_graphs = [list(graph) for graph in small_graphs if len(graph[0]) > 1]
    n_small_graphs = len(small_graphs)
    # average amount of nodes of these smaller graphs (nesting, so [0] again)
    avg_small_graphs = np.mean([len(graph[0]) for graph in small_graphs])
    # % of nodes in main graph
    perc_main_graph = 100 * (n_nodes - n_single_points - n_nodes_small_graphs) / n_nodes
    # return as dict
    graph_info = {'n_nodes': n_nodes, 'n_single_points': n_single_points, 'n_small_graphs': n_small_graphs,
                  'avg_small_graphs': avg_small_graphs, 'n_nodes_small_graphs': n_nodes_small_graphs,
                  'perc_main_graph': perc_main_graph,
                  'edge_len_avg': edge_len_avg, 'edge_len_quantiles': edge_len_quantiles,
                  'avg_edges_neighbours': avg_edges_neighbours,
                  'edge_neighbours_quantiles': edge_neighbours_quantiles}
    return graph_info


def filter_spatial_graph(adata_with_graph, node_fraction_threshold=0.01):
    """
    Filter out single points and small graphs from the spatial graph.
    :param adata_with_graph: anndata object with spatial graph saved in adata.obsp['spatial_connectivities'] as a sparse
    matrix according to the squidpy documentation.
    :param node_fraction_threshold: fraction of nodes in the main graph that is considered a small graph
    :return: filtered anndata object, dict containing filtered points
    """
    # add a new index starting at 0 for the adata
    adata_with_graph.obs['index_sub'] = range(len(adata_with_graph.obs))
    # extract a workable graph from adata
    connectivities = adata_with_graph.obsp['spatial_connectivities']
    # Create a NetworkX graph from the connectivity matrix
    graph = nx.from_scipy_sparse_array(connectivities)
    n_nodes = len(graph.nodes)  # amount of nodes in main graph
    # single points
    single_points = [node for node in graph.nodes if len(list(graph.neighbors(node))) == 0]
    # amount of graphs containing less than n% of data
    small_graphs = [c for c in nx.connected_components(graph) if len(c) < node_fraction_threshold * n_nodes]
    small_graphs = [[list(c)] for c in small_graphs]  # transform to list of lists
    # remove the single points from the small graphs.
    small_graphs = [list(graph) for graph in small_graphs if len(graph[0]) > 1]
    # count all nodes in the small graphs
    all_nodes_small_graphs = [node for graph in small_graphs for node in graph[0]]
    # map single points (in index_sub indexing) to original index
    single_points_original_idx = adata_with_graph.obs.index[single_points].tolist()
    # map small graphs (in index_sub indexing) to original index
    small_graphs_original_idx = adata_with_graph.obs.index[all_nodes_small_graphs].tolist()
    # return dict with filtered points information
    return {'single_points': single_points_original_idx, 'small_graphs': small_graphs_original_idx}


def identify_disconnected_points(adata, library_key=None, node_fraction_threshold=0.01):
    """
    allow for adata to contain one or more samples in library_key column. Filter out single points and small graphs from
    the spatial graph.
    :param adata:
    :param library_key: If multiple library_id, column in anndata.AnnData.obs which stores mapping between library_id and obs.
    :param node_fraction_threshold:
    :return: adata with filtering decisions in adata.obs['disconnected_point']
    """
    # add col to adata with filtering decisions
    adata.obs['disconnected_point'] = 'connected'
    # check if multiple library_id in library_key col
    if library_key:
        if len(adata.obs[library_key].unique()) > 1:
            filtered_points = []
            for library_id in adata.obs[library_key].unique():
                adata_sub = adata[adata.obs[library_key] == library_id].copy()
                points = filter_spatial_graph(adata_sub, node_fraction_threshold)
                filtered_points.append(points)
        else:
            filtered_points = filter_spatial_graph(adata, node_fraction_threshold)
    else:
        filtered_points = filter_spatial_graph(adata, node_fraction_threshold)
    # %% add filtering decisions to adata. if in list single points, set as 'single point'
    single_point_array = adata.obs.index.isin(
        [point for sublist in filtered_points for point in sublist['single_points']])
    # add 'single point' to the disconnected_point column for the single points
    adata.obs['disconnected_point'] = adata.obs['disconnected_point'].mask(single_point_array, 'single point')
    # if in list small graphs, set as 'small graph'
    small_graph_array = adata.obs.index.isin(
        [point for sublist in filtered_points for point in sublist['small_graphs']])
    adata.obs['disconnected_point'] = adata.obs['disconnected_point'].mask(
        small_graph_array, f'small graph < {node_fraction_threshold * 100}%')
    return adata


def evaluate_kde(array, bw_adjust, step_size):
    """
    copied from codex_pointpatternanalysis -> freezing_artefact_preprocessing_functions.py
    evaluate the kde on a grid of points. can have long run times (e.g. overnight) if step_size is low.
    :param array: np.array of shape (2, n) with x and y coordinates
    :param bw_adjust: adjustment factor for the bandwidth, analog to seaborn kdeplot parameter bw_adjust.
    :param step_size: step size for the grid of points to evaluate the kde on.
    :return: kde_populated: np.array of shape (n, m) with the kde value per pixel.
    """
    kernel = stats.gaussian_kde(array)
    kernel.set_bandwidth(bw_method=kernel.factor * bw_adjust)
    # %% grid of points to evaluate the kde on:
    xmin = 0  # x min and y min at 0 to prevent shift between kde and point pattern
    xmax = array[0].max()
    ymin = 0
    ymax = array[1].max()
    X, Y = np.mgrid[xmin:xmax:step_size, ymin:ymax:step_size]
    positions = np.vstack([X.ravel(), Y.ravel()])  # shape (2, 25883)
    kde_populated = np.reshape(kernel(positions).T, X.shape)  # actual kde value per pixel. runs overnight if step=1
    return kde_populated


def spatial_metacluster_filtering(adata, metacluster_key='Metacluster',
                                  spatial_key='position', n_neighbors=10,
                                  majority_threshold=0.7, max_iterations=5,
                                  library_key='dataset_name', key_added='metacluster_filtered'):
    """
    Reassigns cells to the predominant metacluster of their spatial neighbors within each image.

    Parameters:
    - adata: AnnData object containing spatial data.
    - metacluster_key: Key in adata.obs indicating metacluster assignments.
    - spatial_key: Key in adata.obsm containing spatial coordinates.
    - n_neighbors: Number of spatial neighbors to consider.
    - majority_threshold: Proportion threshold to trigger reassignment.
    - max_iterations: Maximum number of iterations for reassignment.
    - library_key: Key in adata.obs that identifies separate spatial datasets (e.g. image/library ID).

    Returns:
    - Updated AnnData object with refined metacluster assignments.
    """
    adata_ref = adata
    adata_ref.obs[key_added] = adata_ref.obs[metacluster_key]

    for image in tqdm(adata_ref.obs[library_key].unique()):
        idx = adata_ref.obs[library_key] == image
        coords = adata_ref.obsm[spatial_key][idx]
        metaclusters = adata_ref.obs.loc[idx, metacluster_key].copy()

        nn = NearestNeighbors(n_neighbors=n_neighbors + 1)  # +1 to include self
        nn.fit(coords)
        _, indices = nn.kneighbors(coords)

        for _ in range(max_iterations):
            changes = 0
            new_metaclusters = metaclusters.copy()

            for i in range(len(metaclusters)):
                neighbor_indices = indices[i][1:]  # skip self (first index)
                neighbor_clusters = metaclusters.iloc[neighbor_indices]
                most_common, count = Counter(neighbor_clusters).most_common(1)[0]
                proportion = count / n_neighbors

                if proportion >= majority_threshold and metaclusters.iloc[i] != most_common:
                    new_metaclusters.iloc[i] = most_common
                    changes += 1

            if changes == 0:
                break
            metaclusters = new_metaclusters

        # Save back to obs
        adata_ref.obs.loc[idx, key_added] = metaclusters

def spatial_metacluster_filtering_incl_diagnostics(
    adata,
    metacluster_key='Metacluster',
    spatial_key='position',
    n_neighbors=10,
    majority_threshold=0.7,
    max_iterations=5,
    library_key='dataset_name',
    key_added='metacluster_filtered'
):
    """
    Reassign cells to the predominant metacluster of their spatial neighbors within each image,
    with per-iteration diagnostics.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing spatial data.
    metacluster_key : str
        Key in adata.obs indicating metacluster assignments.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    n_neighbors : int
        Number of spatial neighbors to consider.
    majority_threshold : float
        Proportion threshold to trigger reassignment.
    max_iterations : int
        Maximum number of reassignment iterations.
    library_key : str
        Key in adata.obs that identifies separate spatial datasets/images.
    key_added : str
        Column name in adata.obs where filtered labels will be stored.

    Returns
    -------
    adata : AnnData
        AnnData object with refined metacluster assignments in `adata.obs[key_added]`.
    diagnostics_df : pd.DataFrame
        Per-image, per-iteration diagnostics dataframe with:
            - image
            - iteration
            - n_changes
            - changed_cell_ids
            - old_labels
            - new_labels
            - transition_counts
            - converged_this_iteration
            - iterations_to_stability
            - stopped_because_max_iterations
    """
    adata_ref = adata.copy()
    adata_ref.obs[key_added] = adata_ref.obs[metacluster_key].copy()

    diagnostics_records = []

    for image in tqdm(adata_ref.obs[library_key].unique()):
        idx = adata_ref.obs[library_key] == image
        coords = adata_ref.obsm[spatial_key][idx]
        metaclusters = adata_ref.obs.loc[idx, key_added].copy()
        cell_ids = metaclusters.index.to_numpy()

        # guard against requesting more neighbors than available cells
        n_cells = coords.shape[0]
        n_neighbors_eff = min(n_neighbors, max(n_cells - 1, 0))

        if n_cells <= 1 or n_neighbors_eff == 0:
            diagnostics_records.append({
                'image': image,
                'iteration': 0,
                'n_changes': 0,
                'changed_cell_ids': [],
                'old_labels': [],
                'new_labels': [],
                'transition_counts': {},
                'converged_this_iteration': True,
                'iterations_to_stability': 0,
                'stopped_because_max_iterations': False
            })
            continue

        nn = NearestNeighbors(n_neighbors=n_neighbors_eff + 1)  # +1 includes self
        nn.fit(coords)
        _, indices = nn.kneighbors(coords)

        converged = False
        iterations_run = 0

        for iteration in range(1, max_iterations + 1):
            changed_cell_ids = []
            old_labels = []
            new_labels = []
            transition_pairs = []

            new_metaclusters = metaclusters.copy()

            for i in range(len(metaclusters)):
                neighbor_indices = indices[i][1:]  # exclude self
                neighbor_clusters = metaclusters.iloc[neighbor_indices]

                if len(neighbor_clusters) == 0:
                    continue

                most_common, count = Counter(neighbor_clusters).most_common(1)[0]
                proportion = count / len(neighbor_clusters)

                old_label = metaclusters.iloc[i]

                if proportion >= majority_threshold and old_label != most_common:
                    new_metaclusters.iloc[i] = most_common
                    changed_cell_ids.append(cell_ids[i])
                    old_labels.append(old_label)
                    new_labels.append(most_common)
                    transition_pairs.append((old_label, most_common))

            n_changes = len(changed_cell_ids)
            transition_counts = dict(Counter(transition_pairs))
            converged_this_iteration = n_changes == 0

            diagnostics_records.append({
                'image': image,
                'iteration': iteration,
                'n_changes': n_changes,
                'changed_cell_ids': changed_cell_ids,
                'old_labels': old_labels,
                'new_labels': new_labels,
                'transition_counts': transition_counts,
                'converged_this_iteration': converged_this_iteration,
                'iterations_to_stability': None,
                'stopped_because_max_iterations': False
            })

            iterations_run = iteration

            if converged_this_iteration:
                converged = True
                break

            metaclusters = new_metaclusters

        # save final assignments back
        adata_ref.obs.loc[idx, key_added] = metaclusters

        # annotate all rows for this image with the final stopping info
        image_rows = [
            i for i, rec in enumerate(diagnostics_records)
            if rec['image'] == image
        ]
        for rec_idx in image_rows:
            diagnostics_records[rec_idx]['iterations_to_stability'] = iterations_run
            diagnostics_records[rec_idx]['stopped_because_max_iterations'] = not converged

    diagnostics_df = pd.DataFrame(diagnostics_records)

    return adata_ref, diagnostics_df


def annotation_corrections(adata):
    """
    Corrects several annotation mistakes found in the data, due to labeling mistakes that were later
    clarified in discussion with experimentalists/documentation. This function deals with codex data only.
    :param adata:
    :return: adata updated
    """

    # if HOE-2257 is found in mouse_ID: drop this entire mouse (is cyclo d4, not act d3)
    adata = adata[adata.obs.query("mouse_id != 'HOE-2257'").index]
    #  BAL-3614: misnamed as TB, is cyclo_d1
    mask_3614 = (adata.obs["mouse_id"] == "BAL-3614")
    adata.obs.loc[mask_3614, "Treatment"] = "Cyclo_d1"
    # BAL-3721 brLNr --> inLNl
    mask_3721 = (adata.obs["mouse_id"] == "BAL-3721") & (adata.obs["Organ"] == "brLNr")
    adata.obs.loc[mask_3721, "Organ"] = "inLNl"
    return adata
