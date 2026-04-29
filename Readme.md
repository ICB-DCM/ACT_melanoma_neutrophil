This is the reproducible code for the manuscript [Longitudinal spatial neutrophil profiling during ACT in murine melanoma reveals distinct lymph node infiltration patterns](https://doi.org/10.64898/2026.02.24.707484), which explores neutrophil and CD8 T cell kinetics and spatial organisation in tumor-draining vs non-draining lymph nodes in murine melanoma under adoptive T cell therapy. The manuscript is available as a preprint on [Bioarxiv](https://doi.org/10.64898/2026.02.24.707484).

# General structure
- Scripts are sorted by data modality
- Per modality, scripts are numbered by the order they should be executed in:
    - raw --> processed (script 0)
        - processed data files can be found in the data repository under 'processed'
    - processed --> postprocessing/analysis (script 1..n)
        - postprocessed/analysed files ready for the visualisation scripts can be found in the data repository under 'processed'
    - postprocessing/analysis --> visualisation for manuscript figures (script n+1).
        - figures as .pdf files can be found in the data repository under 'figures'
- If two scripts/folders have the same number with a different letter (e.g. 2a, 2b), they can be executed in any order. If none of the scripts in a folder have a number, they can be executed in any order.
- Where applicable, the (sub)figure number in the manuscript has been appended to the file name.

# Running the code
- Use either uv sync (preferred) or pip install -r requirements.txt
- In `helper_files/paths_parameters`, update `data_repo_path` to the saved location of the data repository
- Then you can start rerunning from anywhere, the data repository has all intermediate files saved. Per analysis, we recommend running it in the numbered order from preprocessing --> analysis --> visualisation. 


# Data repository 
The data can be found at https://doi.org/10.5281/zenodo.18712088 and is formatted as follows:
act_neutrophil_data_repository
```
├── figures
│   ├── codex
│   │   ├── manuscipt_plots
│   │   ├── overview_intensity_images
│   │   ├── patch_visualisations
│   │   ├── phenotyping
│   │   ├── spatial_analysis
│   │   └── whole_LN_visualisations
│   ├── flow_codex_comparison
│   ├── flow_cytometry
│   │   ├── lymph_node_sizes_all_conditions
│   │   ├── main_conditions
│   │   ├── no_cpg_conditions
│   │   └── preprocessing
│   ├── hematology
│   └── tumor_growth_curves
│       ├── merged
│       ├── per_condition
│       └── revision
├── overview_of_experiments
├── processed
│   ├── codex
│   │   ├── overview_intensity_images
│   │   ├── revision
│   │   └── statistics
│   ├── flow_cytometry
│   │   └── statistics
│   ├── hematology
│   │   └── statistics
│   └── tumor_growth
│       ├── metadata
│       ├── per_experiment
│       └── per_treatment
└── raw
    ├── codex
    │   ├── data.... (folder per image)
    ├── flow_cytometry
    │   ├── cd8_panel
    │   ├── cd8_panel_prism_export
    │   ├── early_gates
    │   ├── pan_immune
    │   └── pmels_in_blood
    ├── hematology
    │   └── per_experiment
    └── tumor_growth
└── Readme.md
```
