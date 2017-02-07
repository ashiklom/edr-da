#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

library(coda)
library(PEcAn.utils)
library(PEcAn.assim.batch)
#--------------------------------------------------------------------------------------------------#


setwd(getwd())  # Set this to your main output dir
pda_output <- readRDS(list.files()[grep("prog_samples.*",list.files())]) # e.g. 'prog_samples.dbh40.temperate.Late_Hardwood.rds'
samples.bt <- PEcAn.assim.batch::autoburnin(pda_output$samples)
samples.bt <- PEcAn.assim.batch::makeMCMCList(samples.bt)

# subset MCMC object
niter(samples.bt)
nvar(samples.bt)
nchain(samples.bt)
#samples.bt[1:5,1:3]
#temp <- samples.bt[1:10,1:3]
start <- niter(samples.bt)*0.35
stop <- niter(samples.bt)
#temp <- PEcAn.assim.batch::makeMCMCList(samples.bt[seq_along(niter(samples.bt)*0.3:niter(samples.bt)),]) # drop the first 30% of the samples
temp <- PEcAn.assim.batch::makeMCMCList(samples.bt[start:stop,1:nvar(samples.bt),drop=TRUE])
niter(temp)

# Plot diagnostics
pdf(file ="trace_plot.pdf",width=8,height=6,onefile=T)
par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
#png(filename ="trace_plot.png", width = 1500, height = 1600, res=150)
plot(samples.bt)
dev.off()

rawsamps <-do.call(rbind, samples.bt)
png(filename="pairs.png", width = 1500, height = 1600, res=150)
pairs(rawsamps)
dev.off()

png(filename="deviance.png", width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(pda_output$deviance))
dev.off()

png(filename="n_eff.png", width = 1500, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(pda_output$n_eff_list))
dev.off()

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot subset
pdf(file ="subset.trace_plot.pdf",width=8,height=6,onefile=T)
par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
plot(temp)
dev.off()

rawsamps <-do.call(rbind, temp)
png(filename="subset.pairs.png", width = 1500, height = 1600, res=150)
pairs(rawsamps)
dev.off()
#--------------------------------------------------------------------------------------------------#
### EOF