source("common.R")

outdir <- "testdir"
dir.create(outdir, showWarnings = FALSE)
link_ed(outdir)

pft <- "late_hardwood"
dbh <- 40
getvar("LAI_CO", dbh, pft)

pp <- c(1.4, 30, 8, 0.01, 0.01)
tv1 <- list(temperate.Late_Hardwood = list(SLA = 50))
tv2 <- list(temperate.Late_Hardwood = list(SLA = 100))
tv3 <- list(temperate.Late_Hardwood = list(SLA = 200))

alb1 <- EDR.run(prospect.param = pp,
               trait.values = tv1,
               output.path = outdir,
               dbh = dbh,
               pft = pft)

alb2 <- EDR.run(prospect.param = pp,
               trait.values = tv2,
               output.path = outdir,
               dbh = dbh,
               pft = pft)

alb3 <- EDR.run(prospect.param = pp,
               trait.values = tv3,
               output.path = outdir,
               dbh = dbh,
               pft = pft)

albmat <- do.call(cbind, list(alb1, alb2, alb3))
matplot(albmat, type='l', col = c("green4", "green3", "green2"))
