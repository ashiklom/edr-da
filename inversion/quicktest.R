#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
source("common.R")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
outdir <- "testdir"
dir.create(outdir, showWarnings = FALSE)
link_ed(outdir)

pft <- "temperate.Late_Hardwood"
dbh <- 20
getvar("LAI_CO", dbh, pft)

data_dir <- normalizePath(paste0('../run-ed/1cohort/dens0.05/dbh',dbh,'/',pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))
file.copy(from = edr_exe_path,
          to = file.path(outdir, 'ed_2.1-opt'), 
          overwrite = TRUE)

# One PFTs: spectra and traits for both
spectra_list <- list(temperate.Late_Hardwood = prospect(c(1.4, 40, 0.01, 0.01), 4, TRUE))
trait.values <- list(temperate.Late_Hardwood = list(clumping_factor = 0.9, orient_factor = -0.2, 
                                                    SLA = 10000,
                                                    b1Bl_large = 0.01, b2Bl_large = 1.488))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Generate ED2 albedo
alb <- EDR(paths = paths,
           spectra_list = spectra_list,
           par.wl = 400:2499,
           nir.wl = 2500,
           datetime = ISOdate(2004, 07, 01, 16, 00, 00),
           trait.values = trait.values,
           output.path = outdir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print(head(alb))
print(tail(alb))
print(range(alb))
print(quantile(alb, c(0.025, 0.5, 0.975)))

waves <- seq(400,2500,1)
png(paste(normalizePath(outdir),'simulated_albedo.png',sep="/"),width=4900, height =2700,res=400)
par(mfrow=c(1,1), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
plot(waves,as.vector(unlist(alb))*100,type="l",lwd=3,ylim=c(0,60),xlab="Wavelength (nm)",ylab="Reflectance (%)",
     cex.axis=1.5, cex.lab=1.7,col="black")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
       "grey80")
lines(waves,as.vector(unlist(alb))*100,lwd=3, col="green4")
dev.off()
#--------------------------------------------------------------------------------------------------#
### EOF