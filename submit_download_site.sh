#!/bin/bash

#qsub -N BH02_site_1-25665 qsub_download_site.sh --site=BH02_site_1-25665
qsub -N BI01_site_1-25676 qsub_download_site.sh --site=BI01_site_1-25676

#for s in sites/*_site_*; do
    #site=${s#sites/}
    #qsub -N $site qsub_download_site.sh --site=$site
#done
