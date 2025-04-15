#!/usr/bin/env python

import importlib
import importlib.util
import logging
import pandas as pd
from bson import ObjectId
import yaml
import numpy as np
import warnings
import os
import envyaml


# This script includes helper functions to parse and normalize various YAML outputs
# from a bioinformatics pipeline such as Bifrost. These include:
# - MLST (Multi-Locus Sequence Typing)
# - Point mutation detection
# - AMR (Antimicrobial Resistance)
# - Plasmid and virulence detection
# - Species classification
# - Assembly quality metrics
# - QC stamping for pipeline success/failure


# ---------------
# CONFIG LOADING
# ---------------
# Define the core package name used for relative configuration paths
PACKAGE_NAME: str = "bifrost_reporter"

# Dynamically locate the package and its directory for loading config files
try:
    spec = importlib.util.find_spec(PACKAGE_NAME)
    if spec is None:
        raise ModuleNotFoundError(f"Package '{PACKAGE_NAME}' not found.")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    PACKAGE_DIR = os.path.dirname(module.__file__)  # Extract directory location

except ModuleNotFoundError as e:
    print(f"Error: {e}")
    PACKAGE_DIR = None
except AttributeError:
    print(f"Error: Could not determine package directory for '{PACKAGE_NAME}'.")
    PACKAGE_DIR = None
except Exception as e:
    print(f"Unexpected error: {e}")
    PACKAGE_DIR = None



def get_config(config_path: str = None):
    """
    Load specified YAML configuration file. If the path is None falls back to 
    the default config in the package directory 

    Parameters:
    ----------
    config_path : str, optional
        Path to the YAML config file. If None, defaults to 'config.default.yaml' in PACKAGE_DIR.

    Returns:
    -------
    dict
        Configuration parameters as a dictionary.
    """
    if config_path is None:
        config_path = os.path.join(PACKAGE_DIR, "config", "config.default.yaml")

    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except Exception as e:
        raise RuntimeError(f"Failed to load config file: {config_path}. Error: {str(e)}")



# ----------------------
# CUSTOM YAML HANDLERS
# ----------------------
# Enable parsing of BSON ObjectIds (common in MongoDB-based YAMLs)
def bson_objectid_constructor(loader, node):
    value = loader.construct_scalar(node)
    return ObjectId(value)

# Register the constructor for custom YAML tags
yaml.add_constructor('!bson.objectid.ObjectId', bson_objectid_constructor)





# ------------------------------------
# PARSERS FOR DIFFERENT BIFROST TOOLS
# ------------------------------------

def parse_mlst(list_files):
    """
    Parse MLST YAML files, extracting the 7 loci alleles and ST (sequence type).
    Returns a DataFrame with one row per sample.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            temp = yaml.load(f, Loader=yaml.Loader)
            if temp["status"] == "Success":
                d[temp["sample"]["name"]] = temp["summary"]["mlst_report"]
            else:
                d[temp["sample"]["name"]] = "N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A"

    # Split CSV-like MLST reports into list of 8 values (ST + 7 loci)
    for key, value in d.items():
        d[key] = value.split(',')

    df = pd.DataFrame.from_dict(d, orient='index', dtype=str)
    df.columns = range(df.shape[1])
    return df



def parse_kmapointmutations(list_files):
    """
    Extract point mutation data from KMA-based pointmutations_tsv reports.
    Returns a multi-indexed DataFrame with mutation records per sample.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("results", {}).get("pointmutations_tsv", {}).get("values", [])
                if not summary:
                    warnings.warn(f"Missing pointmutation summary for: {data['sample']['name']}")
                    continue
                df = pd.DataFrame.from_dict(summary)
                d[data["sample"]["name"]] = df
    if d:
        df = pd.concat(d, names=['Sample Name']).reset_index(level=0)
        df = df.set_index("Sample Name").drop(columns=["#Sample"])
    else:
        df = pd.DataFrame()
    return df



def check_stampers(list_files):
    """
    Evaluate the quality control (QC) summary stamps from YAML.
    If all QC checks are 'pass', mark overall as 'Pass'.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("results", {})
                df = pd.DataFrame.from_dict(summary)
                d[data["sample"]["name"]] = "Pass" if df["status"].eq("pass").all() else "Fail"
            else:
                d[data["sample"]["name"]] = "Requirement Not Met"
    return pd.DataFrame.from_dict(d, orient='index')



def parse_amrfinder(list_files):
    """
    Parse AMRFinder tool output from YAMLs, extracting AMR gene hits.
    Returns long-form DataFrame of gene hits per sample.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("summary", {}).get("output_tsv", [])
                df = pd.DataFrame(summary)
            else:
                df = pd.DataFrame(np.nan, index=range(1), columns=[
                    '% Coverage of reference sequence', '% Identity to reference sequence',
                    'Accession of closest sequence', 'Alignment length', 'Class', 'Contig id',
                    'Element subtype', 'Element type', 'Gene symbol', 'HMM description', 'HMM id',
                    'Method', 'Name of closest sequence', 'Protein identifier',
                    'Reference sequence length', 'Scope', 'Sequence name', 'Start', 'Stop',
                    'Strand', 'Subclass', 'Target length'])
            d[data["sample"]["name"]] = df
    df = pd.concat(d, names=['Sample Name']).reset_index(level=0).set_index("Sample Name")
    return df



def parse_species(list_files):
    """
    Parse Kraken-style species classification output.
    Adds custom column to compute unclassified + top-species proportion.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            if data["status"] == "Success":
                summary = data.get("summary", {})
                d[data["sample"]["name"]] = summary
    df = pd.DataFrame.from_dict(d, orient='index')
    df['sum_unclassified_species1'] = df['percent_unclassified'] + df['percent_classified_species_1']
    return df[[
        "name_classified_species_1", "percent_classified_species_1",
        "name_classified_species_2", "percent_classified_species_2",
        'percent_unclassified', 'sum_unclassified_species1']]



def parse_finder_tools(list_files, ariba_type):
    """
    Generic parser for ARIBA-style outputs (e.g., virulencefinder, plasmidfinder).
    Dynamically adapts to the provided ARIBA result type.
    """
    data_df = {}

    def extract_data(info):
        extracted_data = []
        if info:
            for entry in info:
                extracted_data.append({
                    'GENE': entry.get('GENE', np.nan),
                    '%COVERAGE': entry.get('%COVERAGE', np.nan),
                    '%IDENTITY': entry.get('%IDENTITY', np.nan),
                    'SEQUENCE': entry.get('SEQUENCE', np.nan),
                    'START': entry.get('START', np.nan),
                    'END': entry.get('END', np.nan),
                    'DATABASE': entry.get('DATABASE', np.nan),
                    'ACCESSION': entry.get('ACCESSION', np.nan)
                })
        else:
            extracted_data.append({
                'GENE': np.nan, '%COVERAGE': np.nan, '%IDENTITY': np.nan,
                'SEQUENCE': np.nan, 'START': np.nan, 'END': np.nan,
                'DATABASE': np.nan, 'ACCESSION': np.nan})
        return extracted_data

    for file in list_files:
        with open(file) as f:
            data = yaml.load(f, Loader=yaml.Loader)
            sample_name = data["sample"]["name"]
            if data.get("status") == "Success":
                info = data["summary"].get(ariba_type, [])
                data_df[sample_name] = extract_data(info)
            else:
                data_df[sample_name] = extract_data(None)

    rows = []
    for sample, entries in data_df.items():
        for entry in entries:
            rows.append((sample, entry['GENE'], entry['%COVERAGE'], entry['%IDENTITY'],
                         entry['SEQUENCE'], entry['START'], entry['END'],
                         entry['DATABASE'], entry['ACCESSION']))

    df = pd.DataFrame(rows, columns=[
        'Sample', 'GENE', '%COVERAGE', '%IDENTITY', 'SEQUENCE',
        'START', 'END', 'DATABASE', 'ACCESSION'])

    df[['%COVERAGE', '%IDENTITY', 'START', 'END']] = df[[
        '%COVERAGE', '%IDENTITY', 'START', 'END']].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    df.set_index('Sample', inplace=True)
    return df

def parse_assemblatron(list_files):
    """
    Parse genome assembly metrics from Assemblatron output YAMLs.
    Includes GC content, N50, contig counts, and genome sizes at depth.
    """
    d = {}
    for file in list_files:
        with open(file) as f:
            temp = yaml.load(f, Loader=yaml.Loader)
            if temp["status"] == "Success":
                d[temp["sample"]["name"]] = [
                    temp["summary"]["GC"], temp["summary"]["N50"],
                    temp["summary"]["bin_contigs_at_1x"],
                    temp["summary"]["bin_contigs_at_10x"],
                    temp["summary"]["bin_coverage_at_1x"],
                    temp["summary"]["bin_length_at_1x"],
                    temp["summary"]["bin_length_at_10x"],
                    temp["summary"]["bin_length_at_25x"],
                    temp["summary"]["snp_filter_10x_10%"]]

    df = pd.DataFrame.from_dict(d, orient='index', columns=[
        "GC %", "N50", "Number of contigs (1x cov.)", "Number of contigs (10x cov.)",
        "Average coverage (1x)", "Genome size at 1x depth",
        "Genome size at 10x depth", "Genome size at 25x depth", "Ambiguous sites"])
    return df
