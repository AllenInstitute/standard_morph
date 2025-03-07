import os
import pandas as pd
from collections import defaultdict
from standard_morph.tools import (soma_and_soma_children_qc, 
                                  axon_origination_qc, 
                                  dendrite_origins_qc,
                                  orphan_node_check, 
                                  check_cycles_and_topological_sort,
                                  has_valid_name,
                                  get_soma_mip)

COLUMN_CASTS = {"node_id": int, "parent": int, "compartment": int}
def apply_casts(df, casts):
    for key, typ in casts.items():
        df[key] = df[key].astype(typ)


def get_version():
    with open("standard_morph/__init__.py") as f:
        for line in f:
            if line.startswith("__version__"):
                version = line.split('"')[1]
    return version


class Standardizer():
    
    def __init__(self,
                 path_to_swc: str, 
                 soma_children_distance_threshold: float = 50,
                 swc_separator: str = ' ',
                 valid_filename_format: str = 'None',
                 soma_mip_kwargs: dict = {},
                 ):
        """Class for runing swc standardization

        Args:
            path_to_swc (str): path to swc file
            
            soma_children_distance_threshold (float, optional): maximum distance a soma child can be from the soma. Defaults to 50.
            
            swc_separator (str, optional): Column separator in the swc file. Defaults to ' '.
            
            valid_filename_format (str, optional): choose from ['AIND', 'AIBS', 'None']. When 'None', will not run file name formatting chcek.
            Defaults to 'None'.
            
            soma_mip_kwargs (dict, optional): keyword arguments for running get_soma_mip. When an empty dictionary
            is passed, get_soma_mip will not be ran. See standard_morph.tools.get_soma_mip for kwargs. Defaults to {}.

        Raises:
            ValueError: _description_
        """

        
        if not isinstance(path_to_swc, str):
            type_given = type(path_to_swc)
            err_msg = f"Expected path_to_swc to be a string, but received a {type_given}"
            raise ValueError(err_msg)
        
        self.path_to_swc = path_to_swc
        self.separator = swc_separator
        self.valid_filename_format = valid_filename_format
        self.soma_mip_kwargs = soma_mip_kwargs
        self.soma_children_distance_threshold = soma_children_distance_threshold
        self._validate_inputs()
        self.load_swc()
        report = {
            "errors":[],
            "input_file":path_to_swc,
            "StandardMorphVersion":get_version(),
            "path_to_mip":None
        }
        self.StandardizationReport = report
        self.dfs_sorted_node_ids = []
        
    def _validate_inputs(self):
        #TODO
        # file_content = read_file(self.path_to_file)
        # lines = file_content.splitlines()

        # try:
        #     file_content = read_file(self.path_to_file)
        #     lines = file_content.splitlines()
        # except Exception as e:
        #     err_msg = f"Invalid SWC file: {self.path_to_file}"
        #     raise ValueError(err_msg) from e
        return True
    
    def load_swc(self):
        
        swc_df = pd.read_csv(self.path_to_swc,
                sep=self.separator,
                comment='#',
                header=None,
                names=['node_id', 'compartment', 'x', 'y', 'z', 'r', 'parent'])
        apply_casts(swc_df, casts=COLUMN_CASTS)
        
        parent_ids_dict = dict(zip(swc_df["node_id"], swc_df["parent"]))
        child_ids_dict = defaultdict(list) #{ nid:[] for nid in parent_ids_dict }
        for nid in parent_ids_dict:
            pid = parent_ids_dict[nid]
            if pid != -1:
                child_ids_dict[pid].append(nid)
        
        swc_df = swc_df.set_index("node_id")
        swc_df['node_id'] = swc_df.index
        
        swc_df["number_of_children"] = swc_df.node_id.map(lambda x: len(child_ids_dict[x]))
        

        # this assumes all parent values are in the swc_df (i.e. no orphaned nodes)        
        swc_df['parnet_node_type'] = swc_df.apply(lambda rw: swc_df.loc[rw['parent']]['compartment'] if (rw['parent']!=-1) and (rw['parent'] in swc_df.index.values) else None ,axis=1 )
        
        # create the child look up dictionary
        self._child_ids_dict = child_ids_dict
        self._parent_ids_dict = parent_ids_dict
        self.morph_df = swc_df

    def validate(self):
        
        # check soma and soma children
        soma_and_ch_errors = soma_and_soma_children_qc(self.morph_df, self.soma_children_distance_threshold)
        for s_ch_report in soma_and_ch_errors:
            if s_ch_report['node_ids_with_error'] is not None:
                self.StandardizationReport['errors'].append(s_ch_report) 
                print("Found soma error")
                
        # axon origination QC
        axon_origin_qc_report = axon_origination_qc(self.morph_df)
        if axon_origin_qc_report['node_ids_with_error'] is not None:
            self.StandardizationReport['errors'].append(axon_origin_qc_report) 
        
        # dendrite origination QC
        dendrite_origin_qc_report = dendrite_origins_qc(self.morph_df)
        if dendrite_origin_qc_report['node_ids_with_error'] is not None:
            self.StandardizationReport['errors'].append(dendrite_origin_qc_report) 
        
        # orphan node QC
        orphan_node_report = orphan_node_check(self.morph_df)
        orphaned_node_bool = False
        if orphan_node_report['node_ids_with_error'] is not None:
            orphaned_node_bool = True
            self.StandardizationReport['errors'].append(orphan_node_report) 
        
        # check for loops/cycles
        cycle_error_report, sorted_node_id_dict = check_cycles_and_topological_sort(df=self.morph_df,
                                                                                     child_dict = self._child_ids_dict)
                
        if cycle_error_report['node_ids_with_error']  is not None:
            self.StandardizationReport['errors'].append(cycle_error_report) 
        
        self.sorted_node_id_dict = sorted_node_id_dict
        if self.valid_filename_format != "None":
            filename = os.path.basename(self.path_to_swc)
            filename_report = has_valid_name(filename)
            if filename_report['node_ids_with_error'] is not None:
                self.StandardizationReport['errors'].append(filename_report) 
        
        if self.soma_mip_kwargs != {}:
            mip_kwargs = self.soma_mip_kwargs
            mip_kwargs['morph_df'] = self.morph_df
            mip_path = get_soma_mip(**mip_kwargs)   
            self.StandardizationReport['path_to_mip'] = mip_path
       
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
        
        keep_cols = ['node_id', 'compartment', 'x', 'y', 'z', 'r', 'parent']
        output_df[keep_cols].to_csv(output_swc_path,
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
