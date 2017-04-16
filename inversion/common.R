library(PEcAnRTM)
library(PEcAn.ED2)

edr_exe_path <- readLines("edr_path")

runpath <- function(dbh, pft, dens = 0.05) {
    file.path("../run-ed/1cohort/", 
              paste0("dens", dens),
              paste0("dbh", dbh), 
              pft)
}


getpaths <- function(dbh, pft) {
    run_path <- runpath(dbh, pft)
    paths <- list(ed2in = file.path(run_path, "ED2IN"),
                  history = file.path(run_path, "outputs"))
    return(paths)
}

EDR.run <- function(prospect.param, trait.values, output.path, 
                    dbh, pft) {
    spectra_list <- lapply(prospect.param, prospect, version = 5, include.wl = TRUE)
    out <- EDR(spectra_list = spectra_list,
               trait.values = trait.values,
               paths = getpaths(dbh, pft), 
               par.wl = 400:800,
               nir.wl = 801:2500,
               datetime = ISOdate(2004, 07, 01, 16, 00, 00),
               edr.exe.name = "ed_2.1",
               output.path = output.path)
    return(out)
}

link_ed <- function(output.path, 
                    output.name = "ed_2.1",
                    exe.path = edr_exe_path) {
    file.symlink(exe.path, file.path(output.path, output.name))
}

getvar <- function(varname, dbh, pft) {
    path <- file.path(runpath(dbh, pft), "outputs")
    nc <- ncdf4::nc_open(list.files(path, "history-S", full.names = TRUE))
    out <- ncdf4::ncvar_get(nc, varname)
    ncdf4::nc_close(nc)
    return(out)
}

density_line <- function(dfun, param, p_min, p_max, n = 100) {
    x <- seq(p_min, p_max, length.out = n)
    y <- do.call(dfun, c(list(x=x), as.list(param)))
    plot(x, y, type='l')
}
