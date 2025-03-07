# Standard Morph


Standard Morph is a Python library designed for processing and validating SWC files. It provides tools for standardizing SWC structures and generating quality control (QC) reports. Future direction will include integration with backend databases. 


# Installation

coming soon:  
pip install standard-morph  

For now:  
git clone  

cd standard_morph  
pip install -r requirements.txt  
pip install .  


# Usage
## Example 1: Standardizing a Single SWC File

```python
from standard_morph.Standardizer import Standardizer

# Path to the SWC file
path_to_swc = 'Vip-IRES-Cre_Ai14_IVSCC_-259339.03.01.02_678343034_m.swc'

# Initialize the Standardizer worker
worker = Standardizer(
    path_to_swc=path_to_swc,
    swc_separator=" ",
    soma_children_distance_threshold=50,
    valid_filename_format="None",
    soma_mip_kwargs={}
)

# Validate the SWC file
worker.validate()

# Write the QC'd SWC file with sorted node IDs
worker.write_to_swc(output_swc_path="QCd_File.swc")

# Generate an HTML report for the single cell
worker.report_to_html(report_path="Cell_QC_Report.html")
```
## Example 2: Standardizing multiple SWC Files 
```python
# Define multiple SWC files
path_to_swc_1 = "Vip-IRES-Cre_Ai14_IVSCC_-259339.03.01.02_678343034_m.swc"
path_to_swc_2 = "Vip-IRES-Cre_Ai14_IVSCC_-259339.03.01.02_678343034_m.swc"

all_reports = []
for swc_path in [path_to_swc_1, path_to_swc_2]:
    worker = Standardizer(path_to_swc=swc_path,
                          swc_separator=" ",
                          soma_children_distance_threshold=50,
                          valid_filename_format="None",
                          soma_mip_kwargs={})
    worker.validate()
    this_report = worker.StandardizationReport
    all_reports.append(this_report)
    
# Generate a combined HTML report for multiple cells
create_html_report(data=all_reports, report_path="MultiCellReport.html")
``````


## Example 3: Standardizing a Single SWC File + Soma MIP QC. UNTESTED
```python
from standard_morph.Standardizer import Standardizer

# Path to the SWC file
path_to_swc = 'Vip-IRES-Cre_Ai14_IVSCC_-259339.03.01.02_678343034_m.swc'

soma_mip_kwargs = {
    "image_path":"/soma/omezarr/image/",
    "output_dir":"/where/to/save/mip.png",
    "crop_size":128,
    "mip_depth":10
}
# Initialize the Standardizer worker,
# Note the valid_filename_format will check for AIND style regular expressions
worker = Standardizer(
    path_to_swc=path_to_swc,
    swc_separator=" ",
    soma_children_distance_threshold=50,
    valid_filename_format="AIND",
    soma_mip_kwargs=soma_mip_kwargs
)

# Validate the SWC file
worker.validate()

# Generate an HTML report for the single cell
worker.report_to_html(report_path="Single_Cell_QC_Report_With_Soma_MIP.html")
```
