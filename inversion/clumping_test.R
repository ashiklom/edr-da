#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

source("common.R")
#--------------------------------------------------------------------------------------------------#


outdir <- paste("Clumping_Test", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(outdir)) dir.create(outdir,recursive=TRUE)
link_ed(outdir)

pft <- "temperate.Late_Hardwood"
dbh <- 40
getvar("LAI_CO", dbh, pft)

data_dir <- normalizePath(paste0('../run-ed/1cohort/dens0.05/dbh',dbh,'/',pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))
file.copy(from = edr_exe_path,
          to = file.path(outdir, 'ed_2.1-opt'), 
          overwrite = TRUE)

pp <- c(1.4, 30, 8, 0.01, 0.01)
spectra_list <- list(temperate.Late_Hardwood = prospect(pp, 5, TRUE))
tv1 <- list(temperate.Late_Hardwood = list(clumping_factor = 0.60))
tv2 <- list(temperate.Late_Hardwood = list(clumping_factor = 0.75))
tv3 <- list(temperate.Late_Hardwood = list(clumping_factor = 0.90))

alb1 <- EDR(paths = paths,
            spectra_list = spectra_list,
            par.wl = 400:2499,
            nir.wl = 2500,
            datetime = ISOdate(2004, 07, 01, 16, 00, 00),
            trait.values = tv1,
            output.path = outdir)

alb2 <- EDR(paths = paths,
            spectra_list = spectra_list,
            par.wl = 400:2499,
            nir.wl = 2500,
            datetime = ISOdate(2004, 07, 01, 16, 00, 00),
            trait.values = tv2,
            output.path = outdir)

alb3 <- EDR(paths = paths,
            spectra_list = spectra_list,
            par.wl = 400:2499,
            nir.wl = 2500,
            datetime = ISOdate(2004, 07, 01, 16, 00, 00),
            trait.values = tv3,
            output.path = outdir)
waves <- seq(400,2500,1)
png(paste(outdir,"/",'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
plot(waves,unlist(alb1)*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
     cex.axis=1.5, cex.lab=1.7,col="black")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "grey80")
lines(waves,unlist(alb1)*100,lwd=3, col="green4")
lines(waves,unlist(alb2)*100,lwd=3, col="green3")
lines(waves,unlist(alb3)*100,lwd=3, col="green2")
legend("topright",legend=c(paste0("CF: ",tv1$temperate.Late_Hardwood$clumping_factor),
                           paste0("CF: ",tv2$temperate.Late_Hardwood$clumping_factor),
                           paste0("CF: ",tv3$temperate.Late_Hardwood$clumping_factor)),bty="n",
       col=c("green4","green3","green2"),lty=1,lwd=4)
dev.off()

