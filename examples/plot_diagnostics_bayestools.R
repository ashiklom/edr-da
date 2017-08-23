#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

library(coda)
#library(PEcAn.utils)
#library(PEcAn.assim.batch)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
setwd(getwd())
#pda_output <- readRDS(list.files()[grep("inversion_samples_inprogress*",list.files())]) 
pda_output <- readRDS(list.files()[grep("inversion_samples_progress*",list.files())])  # to remove
#pda_output <- readRDS(list.files()[grep("inversion_samples_inprogress*",list.files())]) 

samples_mcmc <- BayesianTools::getSample(pda_output, coda = TRUE)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Info
niter(samples_mcmc)
nvar(samples_mcmc)
nchain(samples_mcmc)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Plot diagnostics
pdf(file ="trace_plot.pdf",width=8,height=6,onefile=T)
par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
#png(filename ="trace_plot.png", width = 1500, height = 1600, res=150)
plot(samples_mcmc)
dev.off()


#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF





