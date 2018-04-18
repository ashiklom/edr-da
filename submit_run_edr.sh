#!/bin/bash

while read site; do
    #echo "Running $site"
    qsub -N $site qsub_run_edr.sh --site=$site
done <other_site_data/selected_sites
