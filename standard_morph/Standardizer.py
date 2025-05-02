import os
import pandas as pd
import re
from pathlib import Path
from collections import defaultdict
from typing import Optional
from enum import Enum

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
        """
        self.path_to_swc = path_to_swc
        self.input_morphology_df = input_morphology_df
        self.separator = swc_separator
        self.valid_filename_format = valid_filename_format
        self.soma_mip_kwargs = soma_mip_kwargs or {}
        self.soma_children_distance_threshold = soma_children_distance_threshold
        
        self._validate_inputs()
        self.load_data()
        self.extract_node_relationships()

        self.StandardizationReport = {
            "errors": [],
            "input_file": self.path_to_swc,
            "StandardMorphVersion": get_version(),
            "path_to_mip": None,
        }
        print(self.morph_df)

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
            if report.get('node_ids_with_error') is not None:
                self.StandardizationReport['errors'].append(report)

    def validate(self):
        """Run validation checks and build the report."""
        self._append_if_error(soma_and_soma_children_qc(self.morph_df, self.soma_children_distance_threshold))
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
        
        
        
    def report_to_html(self, report_path):
        """
        Generate an HTML report for a single neuron based on `self.StandardizationReport`.

        report_path: path to report html file
        """
        report = self.StandardizationReport

        # Extract details
        neuron_name = os.path.basename(report["input_file"])  # Get filename as neuron name

        path_to_mip = report.get("path_to_mip", None)

        # Handle errors
        errors = report.get("errors", [])
        if errors:
            error_details = "<ul>"
            for error in errors:
                error_details += f"<li><b>{error['test']}:</b> {error['description']}"
                if error["node_ids_with_error"]:
                    error_details += f" (Nodes: {', '.join(map(str, error['node_ids_with_error']))})"
                error_details += "</li>"
            error_details += "</ul>"
            error_class = "error"
        else:
            error_details = "<span class='no-error'>No errors found</span>"
            error_class = "no-error"

        # HTML content
        html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{neuron_name} - QC Report</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    font-size: 18px;
                }}
                table {{
                    width: 100%;
                    border-collapse: collapse;
                    font-size: 18px;
                }}
                table, th, td {{
                    border: 1px solid black;
                }}
                th, td {{
                    padding: 5px;
                    text-align: left;
                }}
                th {{
                    background-color: #f2f2f2;
                }}
                img {{
                    max-width: 256px;
                    height: auto;
                }}
                .error {{
                    color: red;
                    font-weight: bold;
                }}
                .no-error {{
                    color: green;
                    font-weight: bold;
                }}
            </style>
        </head>
        <body>
            <h1>QC Report for {neuron_name}</h1>
            <table>
                <tr>
                    <th>Neuron</th>
                    <th>QC Errors</th>
                    <th>Soma MIP</th>
                    <th>StandardMorph Version</th>
                </tr>
                <tr>
                    <td><b>{neuron_name}</b></td>
                    <td class="{error_class}">{error_details}</td>
                    <td>{f'<img src="{path_to_mip}" alt="Soma MIP">' if path_to_mip else "No Image Available"}</td>
                    <td>{report.get("StandardMorphVersion", "Unknown")}</td>
                </tr>
            </table>
        </body>
        </html>
        """

        # Save report
        with open(report_path, 'w') as f:
            f.write(html_content)

        print(f"Report saved to {report_path}")
        

def create_html_report(data: list, report_path: str) -> None:
    """
    Create an HTML report from the QC validation data. This is used when you have a list of 
    reports that you want to write to an html. Standardizer.report_to_html should be used for
    the case of writing a single cell report to an html.

    Parameters
    ----------
    data : list
        A list of dictionaries, each representing a Standardizer.StandardizationReport
    report_path : str
        The path to save the html report.
    """

    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>SWC QC Report</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                font-size: 18px;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                font-size: 18px;
            }
            table, th, td {
                border: 1px solid black;
            }
            th, td {
                padding: 5px;
                text-align: left;
            }
            th {
                background-color: #f2f2f2;
            }
            img {
                max-width: 128px;
                height: auto;
            }
            .error {
                color: red;
                font-weight: bold;
            }
            .no-error {
                color: green;
                font-weight: bold;
            }
        </style>
    </head>
    <body>
        <h1>SWC QC Report</h1>
        <table>
            <tr>
                <th>Neuron</th>
                <th>QC Errors</th>
                <th>Soma MIP</th>
                <th>StandardMorph Version</th>
            </tr>
    """

    for neuron_data in data:
        neuron_name = os.path.basename(neuron_data["input_file"])  # Extract filename as neuron name
        path_to_mip = neuron_data.get("path_to_mip", None)
        standard_morph_version = neuron_data.get("StandardMorphVersion", "Unknown")
        
        # Handle errors
        errors = neuron_data.get("errors", [])
        if errors:
            error_details = "<ul>"
            for error in errors:
                error_details += f"<li><b>{error['test']}:</b> {error['description']}"
                if error["node_ids_with_error"]:
                    error_details += f" (Nodes: {', '.join(map(str, error['node_ids_with_error']))})"
                error_details += "</li>"
            error_details += "</ul>"
            error_class = "error"
        else:
            error_details = "<span class='no-error'>No errors found</span>"
            error_class = "no-error"

        # Append the row for this neuron
        html_content += f"""
        <tr>
            <td><b>{neuron_name}</b></td>
            <td class="{error_class}">{error_details}</td>
            <td>{f'<img src="{path_to_mip}" alt="Soma MIP">' if path_to_mip else "No Image Available"}</td>
            <td>{standard_morph_version}</td>
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
