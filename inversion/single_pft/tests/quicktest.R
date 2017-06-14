#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

source('config.R')
library(redr)

save_plot <- FALSE #TRUE/FALSE
hidden <- TRUE #TRUE/FALSE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
if (hidden) {
  outdir <- '.quicktest_outdir'
} else {
  outdir <- 'quicktest_outdir'
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
css_df <- data.frame(year = 2000, patch = 1, cohort = 1, dbh = 20, ht = 0, pft = 9,
                     den = 0.05, bdead = 0, balive = 0, lai = -999)
data(pss_ex1)
data(site_ex1)

ed2in_changes <- list(IMONTHA = 06, IDATEA = 29, IYEARA = 2004,
                      IMONTHZ = 07, IDATEZ = 04, IYEARZ = 2004)

try_genrun <- generate_run(prefix = outdir,
                           site_lat = 45.5,
                           site_lon = -90.5,
                           css_df = css_df, 
                           pss_df = pss_df, 
                           site_df = site_df,
                           common_inputs_dir = common_inputs_dir,
                           site_met_dir = site_met_dir,
                           ed_exe_path = ed_exe_path,
                           ed2in_changes = ed2in_changes,
                           # Purge prefix directory if already exists (useful for debugging)
                           RMDIR = TRUE)

message('Running ED...')
test_run <- run_ed(outdir)
message('Done!')
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Test EDR
dir <- paste(outdir)
edr_args <- list(spectra_list = list(temperate.Early_Hardwood = 
                                       PEcAnRTM::prospect(c(1.4, 30, 0.01, 0.01), 4, TRUE)),
                 trait.values = list(temperate.Early_Hardwood = list()),
                 par.wl = 400:2499,
                 nir.wl = 2500,
                 datetime = as.POSIXlt('2004-06-29 12:00:00'))

setup <- setup_edr(dir, edr_exe_path)
message('Running EDR...')
albedo <- run_edr(dir, edr_args)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
message('Done! Albedo summary:')
print(summary(albedo))
print(head(albedo))
print(tail(albedo))
print(range(albedo))
print(quantile(albedo, c(0.025, 0.5, 0.975)))

waves <- seq(400,2500,1)
edr_path <- paste0(outdir,"/edr/")
png(paste(normalizePath(edr_path),'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
plot(waves,as.vector(unlist(albedo))*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
     cex.axis=1.5, cex.lab=1.7,col="black")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "grey80")
lines(waves,as.vector(unlist(albedo))*100,lwd=3, col="green4")
dev.off()
#--------------------------------------------------------------------------------------------------#
### EOF