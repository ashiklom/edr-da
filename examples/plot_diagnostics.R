#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

library(coda)
library(PEcAn.utils)
library(PEcAn.assim.batch)
#--------------------------------------------------------------------------------------------------#

setwd(getwd())
pda_output <- readRDS(list.files()[grep("inversion_samples_inprogress*",list.files())]) 
samples.bt <- PEcAn.assim.batch::autoburnin(pda_output$samples)
samples.bt <- PEcAn.assim.batch::makeMCMCList(samples.bt)

# Info
niter(samples.bt)
nvar(samples.bt)
nchain(samples.bt)

# Plot diagnostics
pdf(file ="trace_plot.pdf",width=8,height=6,onefile=T)
par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
#png(filename ="trace_plot.png", width = 1500, height = 1600, res=150)
plot(samples.bt)
dev.off()

#rawsamps <-do.call(rbind, samples.bt)
#png(filename="pairs.png", width = 2000, height = 1600, res=150)
#pairs(rawsamps)
#dev.off()

png(filename="deviance.png", width = 2000, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(pda_output$deviance))
dev.off()

png(filename="n_eff.png", width = 2000, height = 1600, res=150)
plot(PEcAn.assim.batch::makeMCMCList(pda_output$n_eff_list))
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
