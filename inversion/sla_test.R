#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
source('common.R')
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
outdir <- paste("SLA_Test", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(outdir)) dir.create(outdir,recursive=TRUE)

pft <- "temperate.Late_Hardwood"
css_df <- data.frame(year = 2000, patch = 1, cohort = 1, dbh = 20, ht = 0, pft = 11,
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
                           common_inputs_dir = file.path(edr_da_dir, 'ed-inputs/EDI'),
                           site_met_dir = file.path(edr_da_dir, 'ed-inputs/met3/US-WCr/'),
                           ed_exe_path = ed_exe_path,
                           ed2in_changes = ed2in_changes,
                           RMDIR = TRUE)

message('Running ED...')
test_run <- run_ed(outdir)
message('Done!')
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## EDR
dir <- paste(outdir)

pp <- c(1.4, 30, 8, 0.01, 0.01)
spectra_list <- list(temperate.Late_Hardwood = prospect(pp, 5, TRUE))
tv1 <- list(temperate.Late_Hardwood = list(SLA = 5))
tv2 <- list(temperate.Late_Hardwood = list(SLA = 25))
tv3 <- list(temperate.Late_Hardwood = list(SLA = 50))

edr_args <- list(spectra_list = spectra_list,
                 trait.values = tv1,
                 par.wl = 400:2499,
                 nir.wl = 2500,
                 datetime = as.POSIXlt('2004-06-29 12:00:00'))
setup <- setup_edr(dir, edr_exe_path)
message('Running EDR...')
albedo.1 <- run_edr(dir, edr_args)

edr_args <- list(spectra_list = spectra_list,
                 trait.values = tv2,
                 par.wl = 400:2499,
                 nir.wl = 2500,
                 datetime = as.POSIXlt('2004-06-29 12:00:00'))
message('Running EDR...')
albedo.2 <- run_edr(dir, edr_args)

edr_args <- list(spectra_list = spectra_list,
                 trait.values = tv3,
                 par.wl = 400:2499,
                 nir.wl = 2500,
                 datetime = as.POSIXlt('2004-06-29 12:00:00'))
message('Running EDR...')
albedo.3 <- run_edr(dir, edr_args)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
waves <- seq(400,2500,1)
edr_path <- paste0(outdir,"/edr/")
png(paste(normalizePath(edr_path),'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
plot(waves,unlist(albedo.1)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
     cex.axis=1.5, cex.lab=1.7,col="black")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "grey80")
lines(waves,unlist(albedo.1)*100,lwd=3, col="green4")
lines(waves,unlist(albedo.2)*100,lwd=3, col="green3")
lines(waves,unlist(albedo.3)*100,lwd=3, col="green2")
legend("topright",legend=c(paste0("SLA: ",tv1$temperate.Late_Hardwood$SLA),
                           paste0("SLA: ",tv2$temperate.Late_Hardwood$SLA),
                           paste0("SLA: ",tv3$temperate.Late_Hardwood$SLA)),bty="n",
       col=c("green4","green3","green2"),lty=1,lwd=4)
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
