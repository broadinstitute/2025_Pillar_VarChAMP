import pooch
from pooch import Unzip

## For future data download from Zenodo with hash validation, currently not working due to restricted downloading access

# pillar_data_clinvar = pooch.retrieve(
#     # URL to one of Pooch's test files
#     url="https://zenodo.org/api/records/15558283/files-archive",
#     known_hash=None,
#     fname="pillar_data_15558283.zip",
#     path="raw_inputs/pillar_data", 
#     processor=Unzip(),
#     progressbar=True
# )

# pillar_data_condensed = pooch.retrieve(
#     # URL to one of Pooch's test files
#     url="https://zenodo.org/records/15558283/files/pillar_data_condensed_051325_wREVEL.csv.zip?download=1",
#     known_hash=None,
#     fname="pillar_data_condensed_051325_wREVEL.csv.zip",
#     path="raw_inputs/pillar_data", 
#     # processor=Decompress(),
#     progressbar=True
# )

# pp_readme_v2 = pooch.retrieve(
#     # URL to one of Pooch's test files
#     url="https://zenodo.org/records/15558283/files/PP_Readme_v2.md?download=1",
#     known_hash=None,
#     fname="PP_Readme_v2.md",
#     path="raw_inputs/pillar_data", 
#     # processor=Decompress(),
#     progressbar=True
# )

# file_path = pooch.retrieve(
#     # URL to one of Pooch's test files
#     url="https://zenodo.org/records/15558283/files/summary_pp_df_051325_expanded.tsv?download=1",
#     known_hash=None,
#     fname="summary_pp_df_051325_expanded.tsv",
#     path="raw_inputs/pillar_data", 
#     # processor=Decompress(),
#     progressbar=True
# )