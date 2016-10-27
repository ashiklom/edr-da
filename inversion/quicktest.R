source("common.R")

outdir <- "testdir"
dir.create(outdir, showWarnings = FALSE)
link_ed(outdir)

pft <- "temperate.Late_Hardwood"
dbh <- 20
getvar("LAI_CO", dbh, pft)

pp <- c(1.4, 30, 8, 0.01, 0.01)
tv <- list(temperate.Late_Hardwood = 
           list(clumping_factor = 0.9,
                orient_factor = -0.2, 
                SLA = 10000,
                b1Bl_large = 0.01,
                b2Bl_large = 1.488))
alb <- EDR.run(prospect.param = pp,
               trait.values = tv,
               output.path = outdir,
               dbh = dbh,
               pft = pft)
print(head(alb))
print(tail(alb))
print(range(alb))
print(quantile(alb, c(0.025, 0.5, 0.975)))
