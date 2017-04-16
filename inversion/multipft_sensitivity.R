edr.exe.name <- 'ed_2.1'
pft <- list("temperate.Early_Hardwood","temperate.North_Mid_Hardwood","temperate.Late_Hardwood")
dens <- 0.015
dbh <- 20 # 20, 30 or 40
lai <- getvar("LAI_CO", dbh, pft)
num.cohorts <- 3
multi.pft <- "EMLH"

EDR.multirun <- function(prospect.param, trait.values, paths, output.path) {
    spectra_list <- lapply(prospect.param, prospect, version = 5, include.wl = TRUE)
    out <- EDR(spectra_list = spectra_list,
               trait.values = trait.values,
               paths = paths, 
               par.wl = 400:800,
               nir.wl = 801:2500,
               datetime = ISOdate(2004, 07, 01, 16, 00, 00),
               edr.exe.name = "ed_2.1",
               output.path = output.path)
    return(out)
}


data_dir <- normalizePath(paste0('../run-ed/',num.cohorts,'cohort/dens',dens,'/dbh',dbh,'/',multi.pft))
paths <- list(ed2in = file.path(data_dir, 'ED2IN'),
              history = file.path(data_dir, 'outputs'))

outdir <- 'edr_test'
dir.create(outdir, showWarnings = FALSE)
try_link <- link_ed(outdir)

# Setup PROSPECT parameters
prospect_ver <- 5
pp <- list()
pp[["temperate.Early_Hardwood"]] <- c("N" = 1.4, "Cab" = 30, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.Late_Hardwood"]] <- c("N" = 1.95, "Cab" = 65, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)
pp[["temperate.North_Mid_Hardwood"]] <- c("N" = 1.8, "Cab" = 45, "Car" = 8, "Cw" = 0.01, "Cm" = 0.01)

# Set trait values for each PFT
trait.values <- list()
trait.values[['temperate.Early_Hardwood']] <- list()
trait.values[['temperate.Late_Hardwood']] <- list()
trait.values[['temperate.North_Mid_Hardwood']] <- list()

# Single run to make sure everything works
albedo <- EDR.multirun(prospect.param = pp,
                       trait.values = trait.values,
                       output.path = outdir,
                       paths = paths)

# EDR sensitivity analysis
EDR_sensitivity_trait <- function(sens_values, sens_trait, sens_pft){
    nsens <- length(sens_values)
    albedo_mat <- matrix(NA, 2101, nsens)
    for (i in seq_len(nsens)) {
        print(sens_values[i])
        trait.values[[sens_pft]][[sens_trait]] <- sens_values[i]
        albedo_mat[,i] <- EDR.multirun(prospect.param = pp,
                                       trait.values = trait.values,
                                       output.path = outdir,
                                       paths = paths)
    }
    return(albedo_mat)
}

sensitivity_plot <- function(results_list, values_list, pft_list) {
    nparam <- length(values_list)
    npft <- length(pft_list)
    par(mfrow = c(nparam, npft))
    for (p in seq_len(nparam)) {
        for (f in seq_len(npft)) {
            param <- names(values_list)[p]
            pftf <- pft_list[[f]]
            matplot(results_list[[param]][[pftf]], type='l', main = paste(param, pftf))
        }
    }
}

nsens <- 10
trait_values_list <- list(clumping_factor = seq(0, 1, length.out = nsens),
                          orient_factor = seq(-1, 1, length.out = nsens))

trait_sensitivity_results <- list()
for (p in seq_along(trait_values_list)) {
    for (f in seq_along(pft)) {
        sens_values <- trait_values_list[[p]]
        sens_param <- names(trait_values_list)[p]
        print(sens_param)
        sens_pft <- pft[[f]]
        print(sens_pft)
        trait_sensitivity_results[[sens_param]][[sens_pft]] <- 
            EDR_sensitivity_trait(sens_values, sens_param, sens_pft)
    }
}

sensitivity_plot(trait_sensitivity_results, trait_values_list, pft)

# PROSPECT parameter sensitivity analysis
EDR_sensitivity_prospect <- function(sens_values, sens_param, sens_pft) {
    nsens <- length(sens_values)
    albedo_mat <- matrix(NA, 2101, nsens)
    for (i in seq_len(nsens)) {
        print(sens_values[i])
        pp[[sens_pft]][sens_param] <- sens_values[i]
        albedo_mat[,i] <- EDR.multirun(prospect.param = pp,
                                       trait.values = trait.values,
                                       output.path = outdir,
                                       paths = paths)
    }
    return(albedo_mat)
}

prospect_values_list <- list(N = seq(1, 4, length.out = nsens),
                             Cab = seq(0, 200, length.out = nsens),
                             Car = seq(0, 50, length.out = nsens),
                             Cw = seq(0, 0.05, length.out = nsens),
                             Cm = seq(0, 0.05, length.out = nsens))

prospect_sensitivity_results <- list()
for (p in seq_along(prospect_values_list)) {
    for (f in seq_along(pft)) {
        sens_values <- prospect_values_list[[p]]
        sens_param <- names(prospect_values_list)[p]
        print(sens_param)
        sens_pft <- pft[[f]]
        print(sens_pft)
        prospect_sensitivity_results[[sens_param]][[sens_pft]] <- 
            EDR_sensitivity_prospect(sens_values, sens_param, sens_pft)
    }
}

sensitivity_plot(prospect_sensitivity_results, prospect_values_list, pft)
