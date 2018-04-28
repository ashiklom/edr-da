# ED/EDR simulation analysis

Author: Alexey Shiklomanov

The `redr` package contains functions to facilitate setting up and running EDR test cases and analyses.
The R scripts in other directories contain examples of different EDR-related analyses.

## Multi-side parameter data assimilation

Use scripts in the `scripts/multi_site_pda` directory.

## ED time series ensemble simulation

Use scripts in the `scripts/ensemble` directory.

1. `scripts/ensemble/01_prepare_ensemble.R` -- This script prepares an ensemble for execution, including performing meteorological data processing and setting up ED parameter and namelist files within an intelligent directory structure. Outputs are stored in `ensemble_outputs/<prefix>/`
    - `qsub_ed_<site>.sh` -- Qsub script for an individual site. These are generated automatically, and can be submitted to the cluster with `qsub`.
2. `scripts/ensemble/02_run_edr.R` -- This script runs EDR across the time series generated from the first ED run.
    - `qsub_run_edr.sh` -- Qsub script for running an individual site
    - `submit_run_edr.sh` -- Script to submit all sites.
3. `scripts/ensemble/03_read_all_ts.R` -- Read EDR spectra time series into a single long data frame and save to `sync_data`.
4. `scripts/ensemble/04_landsat_ts.R` -- Convert hyperspectral time series from (3) to Landsat resolution and save.
5. `scripts/ensemble/05_read_history_ts.R` -- Read `history-S` files for each ensemble, and save into a single data frame for each ensemble.
6. `scripts/ensemble/99_prepare_restart.R` -- For ensembles that fail early, prepare a `ED2IN_restart` file that will restart the run from the most recent history file.

## Reload priors

To refresh the internal package priors, re-run the `data-raw/generate_priors.R` script and re-install the package.
