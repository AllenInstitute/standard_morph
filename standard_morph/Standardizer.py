import os
import pandas as pd
import re
from pathlib import Path
from collections import defaultdict
from typing import Optional
from enum import Enum
import warnings

from standard_morph.tools import (
    soma_and_soma_children_qc,
    axon_origination_qc,
    dendrite_origins_qc,
    orphan_node_check,
    check_cycles_and_topological_sort,
    has_valid_name,
    get_soma_mip,
)

SWC_COLUMN_NAMES = ['node_id', 'compartment', 'x', 'y', 'z', 'r', 'parent']
COLUMN_CASTS = {"node_id": int, "parent": int, "compartment": int}


class FilenameFormat(str, Enum):
    NONE = "None"
    AIND = "AIND"
    AIBS = "AIBS"


def apply_casts(df, casts):
    for key, typ in casts.items():
        df[key] = df[key].astype(typ)


def get_version():
    base_dir = Path(__file__).resolve().parents[1]
    pyproject_path = base_dir / "pyproject.toml"

    if not pyproject_path.exists():
        raise FileNotFoundError(f"Could not find pyproject.toml at {pyproject_path}")

    with pyproject_path.open("r", encoding="utf-8") as f:
        for line in f:
            match = re.match(r'version\s*=\s*"(.*?)"', line)
            if match:
                return match.group(1)

    raise ValueError("Version not found in pyproject.toml")


class Standardizer:

    def __init__(
        self,
        path_to_swc: Optional[str] = None,
        input_morphology_df: Optional[pd.DataFrame] = None,
        soma_children_distance_threshold: float = 50,
        swc_separator: str = ' ',
        valid_filename_format: FilenameFormat = FilenameFormat.NONE,
        soma_mip_kwargs: dict = None,
        allow_soma_children_to_branch: bool = False,
    ):
        """
        Class for running SWC standardization.

        Args:
            path_to_swc (Optional[str]): Path to SWC file. Required if no dataframe is passed.
            input_morphology_df (Optional[pd.DataFrame]): Input morphology dataframe.
            soma_children_distance_threshold (float): Distance threshold for soma children.
            swc_separator (str): Column separator in SWC file.
            valid_filename_format (FilenameFormat): Filename format validation option.
            soma_mip_kwargs (dict): Keyword arguments for get_soma_mip function.
            allow_soma_children_to_branch (bool): when True, immediate children of soma are allowed to branch 
        """
        self.path_to_swc = path_to_swc
        self.input_morphology_df = input_morphology_df
        self.separator = swc_separator
        self.valid_filename_format = valid_filename_format
        self.soma_mip_kwargs = soma_mip_kwargs or {}
        self.soma_children_distance_threshold = soma_children_distance_threshold
        self.allow_soma_children_to_branch = allow_soma_children_to_branch
        
        self._validate_inputs()
        self.load_data()
        self.extract_node_relationships()

        self.StandardizationReport = {
            "errors": [],
            "input_file": self.path_to_swc,
            "StandardMorphVersion": get_version(),
            "path_to_mip": None,
        }

    def _validate_inputs(self):
        """Validate the inputs."""
        if self.path_to_swc is not None and not isinstance(self.path_to_swc, str):
            raise ValueError(f"Expected path_to_swc to be a string, got {type(self.path_to_swc)}")

        if self.path_to_swc is None and self.input_morphology_df is None:
            raise ValueError("Either path_to_swc or input_morphology_df is required.")

        column_mapping = {"struct_type": "compartment", "parent_id": "parent", "radius":"r"}
        if self.input_morphology_df is not None:
            if not isinstance(self.input_morphology_df, pd.DataFrame):
                raise ValueError(f"Expected input_morphology_df to be a pandas DataFrame, got {type(self.input_morphology_df)}")

            self.input_morphology_df = self.input_morphology_df.rename(columns=column_mapping)
            missing = [c for c in SWC_COLUMN_NAMES if c not in self.input_morphology_df.columns]
            if missing:
                raise ValueError(f"The following columns are required and missing: {missing}")
        
        if isinstance(self.valid_filename_format, str):
            try:
                FilenameFormat(self.valid_filename_format)
            except ValueError:
                raise ValueError(f"Invalid filename format: {self.valid_filename_format}. Must be one of: {[e.value for e in FilenameFormat]}")


    def load_data(self):
        """Assign input data (from either input df or SWC file) to self._swc_df."""
        if self.path_to_swc:
            swc_df = pd.read_csv(
                self.path_to_swc,
                sep=self.separator,
                comment='#',
                header=None,
                names=SWC_COLUMN_NAMES,
            )
        else:
            swc_df = self.input_morphology_df.copy()

        if swc_df.empty:
            raise ValueError("Input dataframe is empty.")

        self._swc_df = swc_df

    def extract_node_relationships(self):
        """
        Compute structural relationships between nodes (e.g. child counts and parent node types).
        """
        swc_df = self._swc_df
        apply_casts(swc_df, casts=COLUMN_CASTS)

        node_ids = swc_df["node_id"].values
        parent_ids = swc_df["parent"].values
        self._parent_ids_dict = dict(zip(node_ids, parent_ids))

        child_ids_dict = defaultdict(list)
        for nid, pid in zip(node_ids, parent_ids):
            if pid != -1:
                child_ids_dict[pid].append(nid)
        self._child_ids_dict = child_ids_dict

        child_counts = {pid: len(children) for pid, children in child_ids_dict.items()}
        swc_df["number_of_children"] = swc_df["node_id"].map(child_counts).fillna(0).astype(int)

        swc_df = swc_df.set_index("node_id", drop=False)
        swc_df['node_id'] = swc_df.index
        swc_df["parent_node_type"] = swc_df["parent"].map(swc_df["compartment"])

        self.morph_df = swc_df

    def _append_if_error(self, report_list):
        """Append QC report if it contains errors."""
        for report in report_list:
            if report.get('nodes_with_error') is not None:
                self.StandardizationReport['errors'].append(report)

    def validate(self):
        """Run validation checks and build the report."""
        self._append_if_error(soma_and_soma_children_qc(self.morph_df, self.allow_soma_children_to_branch, self.soma_children_distance_threshold))
        self._append_if_error(axon_origination_qc(self.morph_df))
        self._append_if_error(dendrite_origins_qc(self.morph_df))
        self._append_if_error(orphan_node_check(self.morph_df))

        cycle_report, sorted_nodes = check_cycles_and_topological_sort(
            df=self.morph_df,
            child_dict=self._child_ids_dict
        )
        self._append_if_error(cycle_report)
        self.sorted_node_id_dict = sorted_nodes

        if self.valid_filename_format != FilenameFormat.NONE and self.path_to_swc:
            filename = os.path.basename(self.path_to_swc)
            self._append_if_error(has_valid_name(filename))

        if self.soma_mip_kwargs:
            mip_path = get_soma_mip(morph_df=self.morph_df, **self.soma_mip_kwargs)
            self.StandardizationReport["path_to_mip"] = mip_path


    def write_to_swc(self, output_swc_path):
        """write data to swc file

        Args:
            output_swc_path (str): path to output swc
            
        """
        node_mapping = self.sorted_node_id_dict
        output_df = self.morph_df
        apply_casts(output_df, casts=COLUMN_CASTS)
        
        if node_mapping != {}:
            # sort nodes in DFS manor
            output_df['node_id'] = output_df['node_id'].map(node_mapping)
            output_df['parent'] = output_df['parent'].map(node_mapping).fillna(-1)
        
        output_df[SWC_COLUMN_NAMES].to_csv(output_swc_path,
                                    sep=self.separator, 
                                    index=False, 
                                    header=None, 
                                    mode="a" )
          
    def report_to_html(self, report_path, node_display_mode='id'):
        """
        Generate an HTML report for a single neuron based on `self.StandardizationReport`.

        report_path: path to report html file
        node_display_mode: 'id', 'coord', or 'both' to show node ID, coordinates, or both in error list
        """
        report = self.StandardizationReport
        neuron_name = os.path.basename(report["input_file"])
        path_to_mip = report.get("path_to_mip", None)
        errors = report.get("errors", [])

        def format_node(node):
            if node_display_mode == 'id':
                return str(node[0])
            elif node_display_mode == 'coord':
                return f"({node[1]}, {node[2]}, {node[3]})"
            elif node_display_mode == 'both':
                return f"{node[0]}: ({node[1]}, {node[2]}, {node[3]})"
            else:
                warn = f"invalid node_display_mode passed, expected ['id','coord','both']. Got: {node_display_mode}"
                warnings.warn(warn, UserWarning)
                return str(node[0])

        if errors:
            error_details = "<ul>"
            for error in errors:
                error_details += f"<li><b>{error['test']}:</b> {error['description']}"
                if error["nodes_with_error"]:
                    formatted_nodes = ', '.join(format_node(n) for n in error["nodes_with_error"])
                    error_details += f" (Nodes: {formatted_nodes})"
                error_details += "</li>"
            error_details += "</ul>"
            error_class = "error"
        else:
            error_details = "<span class='no-error'>No errors found</span>"
            error_class = "no-error"

        html_content = f"""
        <html>
        <head>
            <title>Standardization Report: {neuron_name}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1 {{ color: #333; }}
                .error {{ color: red; }}
                .no-error {{ color: green; }}
                img {{ max-width: 800px; display: block; margin-top: 20px; }}
            </style>
        </head>
        <body>
            <h1>Standardization Report: {neuron_name}</h1>
            <p><strong>StandardMorph Version:</strong> {report.get("StandardMorphVersion", "Unknown")}</p>
            <p class="{error_class}"><strong>Errors:</strong><br>{error_details}</p>
            {'<img src="' + path_to_mip + '" alt="MIP Image">' if path_to_mip else ''}
        </body>
        </html>
        """

        with open(report_path, 'w') as f:
            f.write(html_content)

        print(f"Report saved to {report_path}")


def create_html_report(data: list, report_path: str, node_display_mode='id') -> None:
    """
    Create an HTML report from a list of StandardizationReports.

    node_display_mode: 'id', 'coord', or 'both' to control display of nodes_with_error
    """
    def format_node(node):
        if node_display_mode == 'id':
            return str(node[0])
        elif node_display_mode == 'coord':
            return f"({node[1]}, {node[2]}, {node[3]})"
        elif node_display_mode == 'both':
            return f"{node[0]}: ({node[1]}, {node[2]}, {node[3]})"
        else:
            warn = f"invalid node_display_mode passed, expected ['id','coord','both']. Got: {node_display_mode}"
            warnings.warn(warn, UserWarning)
            return str(node[0])

    html_content = """
    <html>
    <head>
        <title>Standardization Report Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1 { color: #333; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .error { color: red; }
            .no-error { color: green; }
            img { max-width: 200px; display: block; }
        </style>
    </head>
    <body>
        <h1>Standardization Report Summary</h1>
        <table>
            <tr>
                <th>Neuron Name</th>
                <th>StandardMorph Version</th>
                <th>Errors</th>
                <th>MIP Image</th>
            </tr>
    """

    for neuron_data in data:
        neuron_name = os.path.basename(neuron_data["input_file"])
        path_to_mip = neuron_data.get("path_to_mip", None)
        version = neuron_data.get("StandardMorphVersion", "Unknown")
        errors = neuron_data.get("errors", [])

        if errors:
            error_details = "<ul>"
            for error in errors:
                error_details += f"<li><b>{error['test']}:</b> {error['description']}"
                if error["nodes_with_error"]:
                    formatted_nodes = ', '.join(format_node(n) for n in error["nodes_with_error"])
                    error_details += f" (Nodes: {formatted_nodes})"
                error_details += "</li>"
            error_details += "</ul>"
            error_class = "error"
        else:
            error_details = "<span class='no-error'>No errors found</span>"
            error_class = "no-error"

        mip_img_tag = f'<img src="{path_to_mip}" alt="MIP Image">' if path_to_mip else "N/A"

        html_content += f"""
            <tr>
                <td>{neuron_name}</td>
                <td>{version}</td>
                <td class="{error_class}">{error_details}</td>
                <td>{mip_img_tag}</td>
            </tr>
        """

    html_content += """
        </table>
    </body>
    </html>
    """

    with open(report_path, 'w') as f:
        f.write(html_content)

    print(f"Report saved to {report_path}")
