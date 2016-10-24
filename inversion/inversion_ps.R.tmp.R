source("common.R")

pft <- "late_hardwood"
pft_pecan <- "temperate.Late_Hardwood"
dbh <- 40
lai <- getvar("LAI_CO", dbh, pft)

invert_model <- function(param, seed = 0) {
    outdir <- paste("inversion", seed, sep = ".")
    dir.create(outdir, showWarnings = FALSE)
    try_link <- link_ed(outdir)

    prospect.param <- c('N' = 1.4,
                        'Cab' = 40,
                        'Car' = 10,
                        'Cw' = 0.01,
                        'Cm' = 0.01)

    orient_factor <- param[1]
    clumping_factor <- param[2]
    sla <- param[3]
    b1Bl_large <- param[4]
    b2Bl_large <- param[5]

    trait.values <- list()
    trait.values[[pft_pecan]] <- list(orient_factor = orient_factor,
                                      clumping_factor = clumping_factor,
                                      sla = sla,
                                      b1Bl_large = b1Bl_large,
                                      b2Bl_large = b2Bl_large)

    albedo <- EDR.run(prospect.param = prospect.param,
                      trait.values = trait.values,
                      output.path = outdir,
                      dbh = dbh,
                      pft = pft)

    return(albedo)
}

# Simulate observation
inits <- c("orient_factor" = 0,
           "clumping_factor" = 0.5,
           "sla" = 40,
           "b1Bl_large" = 0.05,
           "b2Bl_large" = 1.45)

obs <- invert_model(inits)# + generate.noise()

prior_def <- list(orient_factor = list("unif", list(-0.5, 0.5)),
                  clumping_factor = list("unif", list(0, 1)),
                  sla = list("lnorm", list(log(20), 1)),
                  b1Bl_large = list("lnorm", list(log(0.25), 0.3)),
                  b2Bl_large = list("lnorm", list(log(0.95), 0.1)))

prior <- function(params) {
    out <- 0
    for (p in seq_along(prior_def)) {
        func <- get(paste0("d", prior_def[[p]][[1]]))
        out <- out + do.call(func, c(unname(params[p]), 
                                     prior_def[[p]][[2]],
                                     TRUE))
    }
    return(out)
}

init_function <- function() {
    out <- sapply(prior_def, 
                  function(x) do.call(get(paste0("r", x[[1]])),
                                      c(1, x[[2]])))
    return(out)
}

param.mins <- c(orient_factor = -0.5,
                clumping_factor = 0,
                sla = 0,
                b1Bl_large = 0,
                b2Bl_large = 0)

param.maxs <- c(orient_factor = 0.5,
                clumping_factor = 1,
                sla = Inf,
                b1Bl_large = Inf,
                b2Bl_large = 1.488) # Hangs if any higher


invert.options <- list(model = invert_model, 
                       nchains = 3,
                       inits.function = init_function,
                       prior.function = prior,
                       ngibbs.max = 100000,
                       ngibbs.min = 500,
                       ngibbs.step = 1000,
                       param.mins = param.mins,
                       param.maxs = param.maxs,
                       adapt = 100,
                       adj_min = 0.1,
                       target = 0.234)

samples <- invert.auto(observed = obs, 
                       invert.options = invert.options,
                       parallel = TRUE)
samples.bt <- PEcAn.assim.batch::autoburnin(samples$samples)
png("trace.noise.png")
plot(samples.bt)
dev.off()

rawsamps <- do.call(rbind, samples.bt)
png("pairs.noise.png")
pairs(rawsamps)
dev.off()
