# Bifrost YAML Data Collection and Processing

This Python package was developed pipeline to **collect**, **parse**, and **filter** genomic analysis results from Bifrost-generated YAML files. It extracts data related to:

- **MLST**
- **AMR genes**
- **Plasmids**
- **Virulence genes**
- **Point mutations**
- **QC stampers**



## Package Structure

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

---

## Input Requirements

###  Configuration file

A configuration file in YAML format able to point to the input file Excel sheets listing the Bifrost sample names.
The package will scan the dictory where the Excel sheet is stored for 

An example of how it is formatted :

```
project:
    name: "EQA"

Illumina:
    new:  "/home/simonescrima/Desktop/bifrost_reporter/input/bifrost_samples/sample_sheet.xlsx"
    original:  "/dpssi/data/Projects/rar_research/proj/Shared_projects/ms_WP9_benchmark/Illumina_original/sample_sheet.xlsx"
```

Where in this case `new` are the sample to be compared again `original`.

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

---

## Usage

Run the following command:
```
python bifrost_reporter/data_collection.py
```



