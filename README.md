# ED and EDR test suite

Author: Alexey Shiklomanov

Creates a set of single-cohort input files and executes a short (2-week) run of ED for each. 
Outputs can then be analyzed with EDR (ED radiative transfer module).

## Setup

Edit the `paths.mk` file with paths to the ED and EDR executables.
These will be used as targets for symlinks.

Then, just run `make` to create the inputs and do the ED runs (will take ~15-60 minutes total, depending on your computer speed).

## EDR analysis

All these scripts are in the `inversion` directory.

Test that EDR works with the `quicktest.R` script.

Perform a sensitivity analysis on the parameters of your choice with `sensitivity.R`.

Pseudodata inversion with `inversion/inversion_ps.R`
