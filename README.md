# Interpretable LCA

Public release of the core code for the sparse latent class analysis project.

## Contents

- `main/`: core R code for the method, simulation workflow, and real-data pipeline.
- `promis/`: placeholders describing data dependencies.

## What is included

- Core method implementation in `main/functions.R`
- Simulation driver in `main/simulation.R`
- Real-data analysis pipeline in `main/real_data_pipeline.R`


## Requirements

- R
- The `poLCA` package

## Basic usage

Run the simulation workflow:

```bash
Rscript main/simulation.R
```

Run the real-data pipeline after placing the expected PROMIS input file in `promis/`:

```bash
Rscript main/real_data_pipeline.R
```

## Notes

- The files use project-relative paths and are intended to be run from the repository root.
