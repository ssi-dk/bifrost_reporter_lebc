# Bifrost YAML Data Collection and Processing

This project provides a modular Python pipeline to **collect**, **parse**, and **filter** genomic analysis results from Bifrost-generated YAML files. It extracts data related to:

- **MLST**
- **AMR genes**
- **Plasmids**
- **Virulence genes**
- **Point mutations**
- **QC stampers**




##  Package Structure

### `data_processing.py`

Contains utility functions to:
- Parse YAML outputs from Bifrost
- Filter based on coverage and identity thresholds
- Handle multiple analysis tools: `ariba`, `amrfinderplus`, `assemblatron`, `kma_pointmutations`, etc.
- Return cleaned results as pandas DataFrames



### `data_collection.py`

Handles:
- Reading a sample sheet (`.xlsx`)
- Validating sample directories and expected result files
- Calling parsing functions from `data_processing`
- Aggregating and returning dataframes
- Optional grouping by sample prefixes (e.g., batch/project IDs)



## Input Requirements

### Config file 

Configuration file in YAML format like below:

```yaml
project:
    name: "Benchmark of new samples against old one"

Illumina:
    new:  "xxx/input/new_samples/sample_sheet.xlsx"
    original:  "xxx/original_samples/sample_sheet.xlsx"
``` 

### Sample Sheet

An Excel file with at least one column named `SampleID` listing the Bifrost result folder names.

### Expected YAML Files per Sample

Each sample directory must contain:

- `__amrfinderplus_fbi.yaml`
- `__ariba_resfinder.yaml`
- `__ariba_virulencefinder.yaml`
- `__ariba_plasmidfinder.yaml`
- `__ariba_mlst.yaml`
- `__assemblatron.yaml`
- `__kma_pointmutations.yaml`   
- `__ssi_stamper.yaml`
- `__reslab_stamper.yaml`
- (others like `__sp_ecoli_fbi.yaml`, `__whats_my_species.yaml` are optional)



## Usage

### 1. Prepare Configuration

Create a config (e.g., YAML file or hardcoded dictionary) like below:

```yaml
Illumina:
  new: /path/to/sample_sheet.xlsx

