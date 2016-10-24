source("common.R")

outdir <- "sensitivity"
dir.create(outdir, showWarnings = FALSE)
link_ed(outdir)

pft <- "late_hardwood"
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
                       pft = "late_hardwood")
        out[,i] <- alb
    }
    return(out)
}


clumping_factor <- edr_sensitivity("clumping_factor", 0, 1)
png("clumping.png")
matplot(clumping_factor, type='l', lty = 1, main = "Clumping factor")
dev.off()

orient_factor <- edr_sensitivity("orient_factor", -0.5, 0.5)
png("orient.png")
matplot(orient_factor, type='l', lty = 1, main = "orient factor")
dev.off()

sla <- edr_sensitivity("SLA", 10, 100)
png("sla.png")
matplot(sla, type='l', lty = 1, main = "SLA")
dev.off()

height_allom_1 <- edr_sensitivity("b1Ht", -1, 1)
matplot(height_allom_1, type='l', lty = 1, main = "height_allom_1")

height_allom_2 <- edr_sensitivity("b2Ht", -1, 1)
matplot(height_allom_2, type='l', lty = 1, main = "height_allom_2")

leaf_allom_1 <- edr_sensitivity("b1Bl_large", 0.1, 0.9)
png("b1Bl_large.png")
matplot(leaf_allom_1, type = 'l', lty = 1, main = "leaf_allom_1")
dev.off()

leaf_allom_2 <- edr_sensitivity("b2Bl_large", 0.7, 1.5)
#png("b2Bl_large.png")
matplot(leaf_allom_2, type = 'l', lty = 1, main = "leaf_allom_2")
#dev.off()

