# Import necessary libraries
import data_processing  # Custom module with parsing and utility functions
import logging
import pandas as pd
pd.set_option('display.max_rows', 1000)  # Show more rows when displaying DataFrames
pd.set_option('display.max_colwidth', None)  # Don't truncate column contents
from bson import ObjectId  # Used for MongoDB object IDs if needed
import yaml
import numpy as np
import glob
import os
from collections import defaultdict



def setup_logging(log_file):
    """
    Function to configure logging settings

    Parameters:
    ----------
    log_file : str
        The path of the log file where log messages will be written
    """
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )



def retrieve_samples(sample_sheet_path):
    """
    Reads an Excel sample sheet and builds full paths to each sample directory.

    Parameters:
    ----------
    sample_sheet_path : str
        Path to the Excel file containing sample IDs

    Returns:
    -------
    pd.Series
        Series of full paths to sample directories
    """
    df = pd.read_excel(sample_sheet_path)
    df["SampleID"] = os.path.dirname(sample_sheet_path) + '/' + df["SampleID"].astype(str)
    return df["SampleID"]



def check_samples(folder_paths):
    """
    Checks each sample directory for the expected output files from Bifrost.

    Parameters:
    ----------
    folder_paths : list of str
        List of sample directory paths

    Returns:
    -------
    dict
        A dictionary showing whether each sample directory and required files exist
    """
    bifrost_results = ["__amrfinderplus_fbi.yaml", 
                       "__ariba_mlst.yaml", 
                       "__ariba_plasmidfinder.yaml",
                       "__ariba_resfinder.yaml", 
                       "__ariba_virulencefinder.yaml", 
                       "__assemblatron.yaml",
                       "__kma_pointmutations.yaml", 
                       "__min_read_check.yaml", 
                       "__reslab_stamper.yaml",
                       "__sp_cdiff_fbi.yaml", 
                       "__sp_ecoli_fbi.yaml", 
                       "__sp_salm_fbi.yaml", 
                       "__ssi_stamper.yaml", 
                       "__whats_my_species.yaml"]

    status = {}
    for folder in folder_paths:
        file_status = {}
        sample_name = os.path.basename(folder)
        
        if os.path.isdir(folder):  # Check if folder exists
            for result in bifrost_results:
                file_name = sample_name + result
                file_path = os.path.join(folder, file_name)
                file_status[file_name] = os.path.isfile(file_path)  # Check if file exists
        else:
            logging.error(f"The folder : {folder} could not be found. Please check again")
        
        status[folder] = {
            'exists': os.path.isdir(folder),
            'files': file_status
        }
    return status



def data_collection_from_dict(sample_dict):
    """
    Gathers and processes data from YAML result files found in each sample directory.

    Parameters:
    ----------
    sample_dict : dict
        Output from check_samples() containing file presence per sample

    Returns:
    -------
    tuple of DataFrames
        Parsed results from different Bifrost analyses
    """
    analysis_files = {}

    for sample_path, sample_info in sample_dict.items():
        if not sample_info["exists"]:
            continue  # Skip if directory doesn't exist

        for file_name, is_present in sample_info["files"].items():
            if is_present:
                try:
                    analysis_name = file_name.strip(".yaml").split("__")[1]
                    analysis_files.setdefault(analysis_name, []).append(os.path.join(sample_path, file_name))
                except:
                    continue  # Ignore parsing errors

    # Parse different analysis results using the data_processing module
    mlst_df = data_processing.parse_mlst(analysis_files.get("ariba_mlst", []))

    plasmid_finder_df = data_processing.parse_finder_tools(analysis_files.get("ariba_plasmidfinder", []), "ariba_plasmidfinder")
    plasmid_finder_df = plasmid_finder_df[(plasmid_finder_df["%COVERAGE"] >= 80) & (plasmid_finder_df["%IDENTITY"] >= 80)]

    resfinder_df = data_processing.parse_finder_tools(analysis_files.get("ariba_resfinder", []), "ariba_resfinder")
    resfinder_df = resfinder_df[(resfinder_df["%COVERAGE"] >= 60) & (resfinder_df["%IDENTITY"] >= 90)]

    virulence_df = data_processing.parse_finder_tools(analysis_files.get("ariba_virulencefinder", []), "ariba_virulencefinder")
    virulence_df = virulence_df[(virulence_df["%COVERAGE"] >= 60) & (virulence_df["%IDENTITY"] >= 90)]

    assemblatron_df = data_processing.parse_assemblatron(analysis_files.get("assemblatron", []))
    kma_df = data_processing.parse_kmapointmutations(analysis_files.get("kma_pointmutations", []))
    amr_df = data_processing.parse_amrfinder(analysis_files.get("amrfinderplus_fbi", []))
    ssi_stamper_df = data_processing.check_stampers(analysis_files.get("ssi_stamper", []))
    reslab_stamper_df = data_processing.check_stampers(analysis_files.get("reslab_stamper", []))

    return (mlst_df,
            plasmid_finder_df,
            resfinder_df, 
            virulence_df,
            assemblatron_df, 
            kma_df, 
            amr_df, 
            ssi_stamper_df, 
            reslab_stamper_df
    ) 



def extract_prefix(sample_name):
    """
    Extracts the first three underscore-separated components from a sample name.

    Parameters:
    ----------
    sample_name : str

    Returns:
    -------
    str
        Prefix composed of first 3 elements (e.g., "HER_BTP_WGS")
    """
    return "_".join(sample_name.split("_")[:3])



if __name__ == '__main__':
    # Load environment config with sample paths
    config = data_processing.get_config()

    # Load and check new Illumina sample directories and files
    new_illumina = check_samples(retrieve_samples(config["Illumina"]["new"]))

    # Process the validated files from new Illumina samples
    l_new = data_collection_from_dict(new_illumina)

    # Print result for quick inspection
    print(l_new)

    


    # l_dict_dfs=[]
    
    # for df in l_new:
    #     try:
    #         df["Prefix"] = df.index.to_series().apply(extract_prefix)
    #         dict_dfs = {prefix: sub_df.drop(columns=["Prefix"]) for prefix, sub_df in df.groupby("Prefix")}
    #         l_dict_dfs.append(dict_dfs)
    #         dfs ={}
    #     except:
    #         continue
    #         #print( df.index.to_series())
    #         #print(df)
    

    # print(l_dict_dfs)    
    # mlst_dict = {key: pd.concat([df, l_og[0]]) for key, df in l_dict_dfs[0].items()}
    # plasmid_finder = {key: pd.concat([df, l_og[1]]) for key, df in l_dict_dfs[1].items()}
    

    # # Check one of the updated DataFrames
    #print(plasmid_finder["HER_BTP_WGS"])
    

