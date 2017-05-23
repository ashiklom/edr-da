#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

source("common.R")
#--------------------------------------------------------------------------------------------------#

outdir <- paste("Sensitivity_Test", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(outdir)) dir.create(outdir,recursive=TRUE)
link_ed(outdir)

pft <- "temperate.Late_Hardwood"
pft_pecan <- "temperate.Late_Hardwood"
dbh <- 20 # 20, 30 or 40
lai <- getvar("LAI_CO", dbh, pft)

data_dir <- normalizePath(paste0('../run-ed/1cohort/dens0.05/dbh',dbh,'/',pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))
file.copy(from = edr_exe_path,
          to = file.path(outdir, 'ed_2.1-opt'), 
          overwrite = TRUE)

pp <- c(1.4, 30, 8, 0.01, 0.01)
spectra_list <- list(temperate.Late_Hardwood = prospect(pp, 5, TRUE))

edr_sensitivity <- function(param_name, minval, maxval, steps = 7) {
    out <- matrix(numeric(), 2101, steps)
    param_seq <- seq(minval, maxval, length.out = steps)
    for (i in seq_along(param_seq)) {
        val <- param_seq[i]
        trait.values <- list()
        trait.values[[pft_pecan]] <- list()
        trait.values[[pft_pecan]][[param_name]] <- val
        alb <- EDR(paths = paths,
                    spectra_list = spectra_list,
                    par.wl = 400:2499,
                    nir.wl = 2500,
                    datetime = ISOdate(2004, 07, 01, 16, 00, 00),
                    trait.values = trait.values,
                    output.path = outdir)
        out[,i] <- alb
    }
    return(out)
}

## Results
clumping_factor <- edr_sensitivity("clumping_factor", 0, 1)
png(paste0(outdir,"/","clumping.png"))
matplot(clumping_factor, type='l', lty = 1, main = "Clumping factor")
dev.off()

orient_factor <- edr_sensitivity("orient_factor", -0.5, 0.5)
png(paste0(outdir,"/","orient.png"))
matplot(orient_factor, type='l', lty = 1, main = "orient factor")
dev.off()

sla <- edr_sensitivity("SLA", 10, 100)
png(paste0(outdir,"/","sla.png"))
matplot(sla, type='l', lty = 1, main = "SLA")
dev.off()

height_allom_1 <- edr_sensitivity("b1Ht", -1, 1)
png(paste0(outdir,"/","b1Ht.png"))
matplot(height_allom_1, type='l', lty = 1, main = "height_allom_1")
dev.off()

height_allom_2 <- edr_sensitivity("b2Ht", -1, 1)
png(paste0(outdir,"/","b2Ht.png"))
matplot(height_allom_2, type='l', lty = 1, main = "height_allom_2")
dev.off()

leaf_allom_1 <- edr_sensitivity("b1Bl_large", 0.1, 0.9)
png(paste0(outdir,"/","b1Bl_large.png"))
matplot(leaf_allom_1, type = 'l', lty = 1, main = "leaf_allom_1")
dev.off()

leaf_allom_2 <- edr_sensitivity("b2Bl_large", 0.7, 1.5)
png(paste0(outdir,"/","b2Bl_large.png"))
matplot(leaf_allom_2, type = 'l', lty = 1, main = "leaf_allom_2")
dev.off()


## multi-panel
waves <- seq(400,2500,1)
png(paste(outdir,"/",'sensitivity_test.png',sep="/"),width=4900, height =2700,res=400)
par(mfrow=c(3,3), mar=c(4.3,4.5,1.0,1), oma=c(0.1,0.1,0.1,0.1)) # B L T R
matplot(waves,clumping_factor, type='l', lty = 1, main = paste0("DBH: ",dbh," Param: Clumping factor"))
matplot(waves,orient_factor, type='l', lty = 1, main = paste0("DBH: ",dbh," Param: orient factor"))
matplot(waves,sla, type='l', lty = 1, main = paste0("DBH: ",dbh," Param: SLA"))
matplot(waves,height_allom_1, type='l', lty = 1, main = paste0("DBH: ",dbh," Param: height_allom_1"))
matplot(waves,height_allom_2, type='l', lty = 1, main = paste0("DBH: ",dbh," Param: height_allom_2"))
matplot(waves,leaf_allom_1, type = 'l', lty = 1, main = paste0("DBH: ",dbh," Param: leaf_allom_1"))
matplot(waves,leaf_allom_2, type = 'l', lty = 1, main = paste0("DBH: ",dbh," Param: leaf_allom_2"))
dev.off()
