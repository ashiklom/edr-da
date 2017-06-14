library(ncdf4)

nc <- nc_open("testdir/history-S-2004-07-01-120000-g01.h5")

ncvar_get(nc, "LAI_CO")
