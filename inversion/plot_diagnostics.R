#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#--------------------------------------------------------------------------------------------------#


setwd('')  # Set this to your main output dir
pda_output <- readRDS('') # e.g. 'prog_samples.dbh40.temperate.Late_Hardwood.rds'

library(coda)
library(PEcAn.utils)
library(PEcAn.assim.batch)
samples.bt <- PEcAn.assim.batch::autoburnin(input.pda.data$samples)
samples.bt <- PEcAn.assim.batch::makeMCMCList(samples.bt)

par(mfrow=c(1,1), mar=c(2,2,0.3,0.4), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
png(filename ="trace_plot.png", width = 1500, height = 1600, res=150)
plot(samples.bt)
dev.off()

rawsamps <- do.call(rbind, samples.bt)
png(filename="pairs.png", width = 1500, height = 1600, res=150)
pairs(rawsamps)
dev.off()

png(filename="deviance.png", width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$deviance))
dev.off()

png(filename="n_eff.png", width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(input.pda.data$n_eff_list))
dev.off()

#--------------------------------------------------------------------------------------------------#
### EOF