# Cutting out the middleman: Calibrating and validating a vegetation model using remotely sensed surface reflectance

**Source code for the analysis and figures.**

To reproduce the analysis for the GMD paper:

1. Clone this repository to your local machine.
2. Download the `gmd-data.tar.gz` file from OSF (https://osf.io/b6umf/) and extract its contents to the root of the `edr-da` repository. This should create the following folders: `aviris`, `other_site_data`.
3. Run the script `scripts/inversion_r_multi_site.R` to perform the multi-site calibration. Note that the exact inversion configuration is controlled by command line flags; the exact version requires the additional argument `hetero`. This script takes a long time, so you may want to run it on an HPC cluster.
4. Once (3) has finished, you can post-process the results by running `make all`. This will run the post-processing R `drake` workflow, and then build the manuscript PDF.
