# Cutting out the middleman: Calibrating and validating a vegetation model using remotely sensed surface reflectance

**Source code for the analysis and figures.**

To reproduce the analysis for the GMD paper:

1. Clone this repository to your local machine.
2. Download the `gmd-data.tar.gz` file from OSF (https://osf.io/b6umf/) and extract its contents to the root of the `edr-da` repository. This should create the following folders: `aviris`, `other_site_data`.
3. Run the numbered scripts in the `scripts/` directory to setup and run the multi-site calibration.
   - The last script takes a long time, so you may want to run it on an HPC cluster. Note that the third script creates an `RData` file with all the data needed for the inversion, so you can run the first three scripts locally and then upload the `RData` file to your HPC and only run the 4th script on the HPC.
   - Alternatively, to get the exact MCMC samples used in the manuscript, download the `results.tar.gz` file from the same OSF repository. Once this file is unzipped, you should be able to run step 4 below and exactly reproduce the results of the manuscript.
4. Once (3) has finished, you can post-process the results by running `make all`. This will run the post-processing R `drake` workflow, and then build the manuscript PDF.
