source("common.R")

outdir <- paste("sensitivity", format(Sys.time(), format="%Y%m%d_%H%M%S"), sep = "_")
if (! file.exists(outdir)) dir.create(outdir,recursive=TRUE)
link_ed(outdir)

pft <- "temperate.Late_Hardwood"
pft_pecan <- "temperate.Late_Hardwood"
dbh <- 40
lai <- getvar("LAI_CO", dbh, pft)
prospect.param <- c('N' = 1.4,
                    'Cab' = 40,
                    'Car' = 10,
                    'Cw' = 0.01,
                    'Cm' = 0.01)

edr_sensitivity <- function(param_name, minval, maxval, steps = 7) {
    out <- matrix(numeric(), 2101, steps)
    param_seq <- seq(minval, maxval, length.out = steps)
    for (i in seq_along(param_seq)) {
        val <- param_seq[i]
        trait.values <- list()
        trait.values[[pft_pecan]] <- list()
        trait.values[[pft_pecan]][[param_name]] <- val
        alb <- EDR.run(prospect.param = prospect.param,
                       trait.values = trait.values,
                       output.path = outdir,
                       dbh = 40,
                       pft = paste0(pft))
        out[,i] <- alb
    }
    return(out)
}


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

