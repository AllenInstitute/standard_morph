import re
import math
import os
import re
from collections import defaultdict, deque
from typing import List, Dict
import pandas as pd
import warnings

# t = """    
#     8. There are no nodes that share identical x,y,z coordinates
#         -- due to the way vaa3d works
#         -- REMOVE yay        

#     1. Make sure a soma node exists
#         -- keep
        
#     2. Make sure the soma is not > distance_threshold from it's children
#         -- keep but may not be relevant, manual intervene
        
#     3. The immediate children of the soma are not furcation nodes
#         -- keep manmual interevene
        
#     4. There is only one root node (i.e. one node who has parent = -1)
#         --keep rolled inb with 1. 
        
#     5. There is only 1 place of axon origination
#         --keep, manual intervene
#     10. All axon nodes have parent of either axon, soma, or basal dendrite.
#         --keep, manual intervene
        
#     6. All nodes parents' exist in the tree (with the exception of a parent=-1, this will go towards the root node count)
#         --keep, manual intervene
        
#     9. All apical/basal dendrite nodes have parent of either soma or apical/basal dendrite respectively
#         --keep, manual intervene
    
#     11. There are no loops
#         --keep, manual intervene
    
#     12. The soma ID is 1, nodes sorted
#         --keep, programatic fix
    
#     """
def euclidean_distance(x, y):
    return sum((a - b) ** 2 for a, b in zip(x, y)) ** 0.5
     
def soma_and_soma_children_qc(df, allow_soma_children_to_branch=False, soma_child_distance_threshold=50):
    """
    Perform QC checks on the soma and its immediate children.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing neuron structure data with columns: 'compartment', 'parent', 'x', 'y', 'z', 'number_of_children'.
    allow_soma_children_to_branch : bool, optional
        When True, immediate children nodes of the soma are permitted to branch.
    soma_child_distance_threshold : float, optional
        Maximum allowed distance (in microns) for soma children from the soma. Default is 50.
    Returns
    -------
    list of dict
        List of QC test results, each as a dictionary with test metadata and node IDs with errors.
    """
    soma_df = df[(df['compartment'] == 1) & (df['parent'] == -1)]
    n_somas = soma_df.shape[0]

    soma_children_furcation_test = {
        "test": "SomaChildrenFurcation",
        "description": "Children nodes of the soma should not branch. Returned node IDs are soma children that branch.",
        "nodes_with_error": None
    }

    soma_children_distance_test = {
        "test": "SomaChildrenDistance",
        "description": f"Immediate children of the soma should be within {soma_child_distance_threshold} microns. Returned node IDs exceed that distance.",
        "nodes_with_error": None
    }

    n_soma_error = {
        "test": "NumberOfSomas",
        "description": "There should be exactly one node with type=1 and parent=-1. Returned node IDs do not meet this criterion.",
        "nodes_with_error": None
    }

    if n_somas != 1:
        warn_msg = f"{n_somas} somas detected, unable to run soma_and_soma_children_qc. Verify axon + dend. swc files merged correctly"
        warnings.warn(warn_msg, UserWarning)
        if not soma_df.empty:
            n_soma_error['nodes_with_error'] = list(soma_df[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))
        else:
            n_soma_error['nodes_with_error'] = [(1,0,0,0)]
            
    else:
        soma_node_id = soma_df['node_id'].values[0]
        soma_children_df = df[df['parent'] == soma_node_id].copy()
        problem_children_branching = soma_children_df[soma_children_df['number_of_children'] != 1]
        if not problem_children_branching.empty:
            if not allow_soma_children_to_branch:
                soma_children_furcation_test['nodes_with_error'] = list(problem_children_branching[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))

        soma_coord = [soma_df['x'].values[0], soma_df['y'].values[0], soma_df['z'].values[0]]
        soma_children_df['distance_from_soma'] = soma_children_df.apply(
            lambda rw: euclidean_distance([rw.x, rw.y, rw.z], soma_coord), axis=1
        )
        problem_children_distance = soma_children_df[
            soma_children_df['distance_from_soma'] > soma_child_distance_threshold
        ]
        if not problem_children_distance.empty:
            soma_children_distance_test['nodes_with_error'] = list(problem_children_distance[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))
    
    return [n_soma_error, soma_children_furcation_test, soma_children_distance_test]


def axon_origination_qc(df):
    """
    Check that the axon originates from an appropriate compartment (axon, soma, or basal dendrite).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing neuron structure data with axon compartment information.

    Returns
    -------
    list of dict
        QC test result for axon origination.
    """
    axon_origination_error = {
        "test": "AxonOrigins",
        "description": "Axon should originate from a single location and stem from axon, soma, or basal dendrite.",
        "nodes_with_error": None
    }

    axon_df = df[df['compartment'] == 2]
    axon_origination_df = axon_df[axon_df['parent_node_type'] != 2]
    num_origins = axon_origination_df.shape[0]

    if num_origins != 1 or int(axon_origination_df['parent_node_type'].values[0]) not in [1, 3]:
        axon_origination_error['nodes_with_error'] = list(axon_origination_df[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))

    return [axon_origination_error]


def orphan_node_check(df):
    """
    Identify nodes with missing or invalid parent assignments. This assumes that an orphaned
    node will have an empty 'parent_node_type' column value, because the parent
    node does not exist in the tree.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with neuron data including parent-child relationships.

    Returns
    -------
    list of dict
        QC test result for orphaned nodes.
    """
    orphaned_node_error = {
        "test": "OrphanedNodes",
        "description": "Nodes with missing or incorrect parent assignments.",
        "nodes_with_error": None
    }

    orphaned_nodes = df[(df['parent_node_type'].isnull()) & (df['parent'] != -1)]
    if not orphaned_nodes.empty:
        orphaned_node_error['nodes_with_error'] = list(orphaned_nodes[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))

    return [orphaned_node_error]


def dendrite_origins_qc(df):
    """
    Validate that apical and basal dendrites originate from appropriate compartments.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing dendritic compartment information.

    Returns
    -------
    list of dict
        QC test result for dendrite origin issues.
    """
    dend_df = df[df['compartment'].isin([3, 4])]
    invalid_dend_origins = []

    for comp, sub_df in dend_df.groupby("compartment"):
        origins_df = sub_df[sub_df['parent_node_type'] != comp]
        invalid_origins = origins_df[origins_df['parent_node_type'] != 1]
        invalid_dend_origins.extend(invalid_origins[['node_id', 'x', 'y', 'z']].itertuples(index=False, name=None))

    if not invalid_dend_origins:
        invalid_dend_origins = None

    return [{
        "test": "DendriteOrigins",
        "description": "Dendritic nodes should originate from soma or matching dendrite type.",
        "nodes_with_error": invalid_dend_origins
    }]


def check_cycles_and_topological_sort(df, child_dict):
    """
    Check for cycles in the graph and return a topological sort of nodes.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'node_id' and 'parent' columns.
    child_dict : dict
        Dictionary mapping each node ID to its children.

    Returns
    -------
    tuple
        ([QC report dict], node mapping dict)
    """
    report = {
        "test": "CheckForLoops",
        "description": "Check if loops exist in reconstruction.",
        "nodes_with_error": None
    }

    root_nodes = set(df.loc[df['parent'] == -1, 'node_id'])
    visited = []
    seen_ids = set()

    for root_id in root_nodes:
        stack = deque([root_id])
        while stack:
            current_node = stack.popleft()
            if current_node in seen_ids:
                report["nodes_with_error"] = [(1,0,0,0)]
                return [report], {}
            visited.append(current_node)
            seen_ids.add(current_node)
            for ch_no in child_dict.get(current_node, []):
                stack.appendleft(ch_no)

    node_mapping = {old_id: new_id for new_id, old_id in enumerate(visited, 1)}

    return [report], node_mapping


def has_valid_name(swc_file: str, name_format='AIND'):
    """
    Validate the format of an SWC filename.

    Parameters
    ----------
    swc_file : str
        Filename to validate.
    name_format : str, optional
        Expected naming convention. Default is 'AIND'.

    Returns
    -------
    list of dict
        QC report on filename formatting. If the file name is formatted appropriately, 
        the returned error report will have None for the 'nodes_with_error' field. 
        If the name is not formated correctly, the report will have [1] in the 
        'nodes_with_error' field.
    """
    if name_format == 'AIND':
        pattern = r"^N\d{1,}(-|_)\d{6}(-|_)(axon|dendrite|dendrites)(-|_)([A-Za-z]{2,3}|consensus)\.swc$"

    else:
        print("TODO: AIBS FILE NAME RE CHECK")
        pattern = ""

    is_valid = re.match(pattern, swc_file, re.IGNORECASE) is not None

    filename_error_report = {
        "test": "FileNameFormat",
        "description": "Check if file is named correctly.",
        "nodes_with_error": None
    }

    if not is_valid:
        filename_error_report["nodes_with_error"] = [(1,0,0,0)]

    return [filename_error_report]

def get_soma_mip(
        image_path: str,
        morph_df: pd.DataFrame,
        mip_path: str,
        crop_size: int = 128,
        mip_depth: int = 10,
) -> str:
    """
    Get a MIP of the soma coordinate from the Zarr image.

    Parameters
    ----------
    image_path : str
        The path to the Zarr image file.
    graph : NeuronGraph
        The neuron graph.
    mip_path : str
        The output .png file to save the mip
    crop_size : int, optional
        The size of the crop around the soma coordinate (default is 128).
    mip_depth : int, optional
        The number of slices on each side for the MIP

    Returns
    -------
    str
        The path to the saved MIP image.
    """
    try:
        import imageio
        import s3fs
        import zarr
        from skimage.exposure import rescale_intensity
    except:
        msg = "imageio, s3fs, zarr and skimage needed to run get_soma_mip"
        raise ValueError(msg)

    # Load the image
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name='us-west-2'))
    store = s3fs.S3Map(root=image_path, s3=s3, check=False)
    z = zarr.open(store, mode='r')
    arr = z['0']

    # get OME-Zarr scale metadata from root group
    metadata = z.attrs.asdict()
    scale = metadata['multiscales'][0]['datasets'][0]['coordinateTransformations'][0]['scale']

    # root = [node for node in graph.nodes if graph.in_degree(node) == 0][0]
    soma_df = morph_df[(morph_df['compartment'] == 1) & (morph_df['parent'] == -1)]
    n_somas = soma_df.shape[0]
    if n_somas >1:
        warn_msg = f"There were {n_somas} somas detected, using the first one to take soma MIP"
        warnings.warn(warn_msg, UserWarning)
        
    elif n_somas<1:
        warn_msg = "No somas found, can not take soma MIP" 
        warnings.warn(warn_msg, UserWarning)
        return None 
             
    z = int(soma_df["z"].values[0] / scale[2])
    y = int(soma_df["y"].values[0] / scale[3])
    x = int(soma_df["x"].values[0] / scale[4])

    s = math.ceil(crop_size / 2)

    mip = arr[0, 0, z - mip_depth:z + mip_depth, y - s:y + s, x - s:x + s].max(axis=0)

    mip = rescale_intensity(mip, out_range=(0, 255)).astype('uint8')

    imageio.imwrite(mip_path, mip)

    return mip_path


# @pd.api.extensions.register_dataframe_accessor("morph")
# class MorphAccessor:
#     def __init__(self, pandas_obj):
#         self._validate(pandas_obj)
#         self._obj = pandas_obj

#     @staticmethod
#     def _validate(obj):
#         # verify there is a column latitude and a column longitude
#         col_list = ['node_id', 'compartment', 'x', 'y', 'z', 'r', 'parent']
#         if not all([c in obj.columns for c in col_list]):
#             raise AttributeError("All columns must be present:\n{}".format(col_list))

#     def get_nodes_by_type(self, node_type_list):
#         "get sub dataframe for nodes of a certain type"
#         return self._obj[self._obj['compartment'].isin(node_type_list)]

    
#     def plot(self):
#         # plot this array's data on a map, e.g., using Cartopy
#         pass





# def find_duplicate_coordinates(df):
#     """
#     Identify duplicate node coordinates in the dataframe.
    
#     Parameters:
#     df (pd.DataFrame): DataFrame containing node information with columns 'x', 'y', 'z', and 'node_id'.
    
#     Returns:
#     list or None: A list of node ID lists where each inner list contains IDs of nodes sharing the same (x, y, z) coordinates.
#                   Returns None if no duplicates are found.
#     """
#     duplicates = df.groupby(['x', 'y', 'z'])['node_id'].apply(list)
#     dup_ids = [ids for ids in duplicates if len(ids) > 1]
    
#     if dup_ids:
#         return dup_ids
    
#     return None

