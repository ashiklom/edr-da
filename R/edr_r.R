#' Run ED two-stream radiative transfer model
#'
#' Wrapper around [sw_two_stream] that also calls PROSPECT 5 (see
#' [PEcAnRTM::prospect()]) to generate the leaf reflectance and
#' transmittance spectra, [hapke_soil()] to generate the soil spectrum
#' given `soil_moisture`, and a "flat" incident solar spectra for
#' direct and diffuse light based on `direct_sky_frac`.
#'
#' @inherit sw_two_stream params return
#' @inherit hapke_soil params
#' @param N PROSPECT 5 effective number of leaf mesophyll layers (>1; npft)
#' @param Cab PROSPECT 5 leaf chlorophyll content (ug cm-2; npft)
#' @param Car PROSPECT 5 leaf carotenoid content (ug cm-2; npft)
#' @param Cw PROSPECT 5 leaf water content (g cm-2; npft)
#' @param Cm PROSPECT 5 leaf dry matter content (g cm-2; npft)
#' @param direct_sky_frac Fraction of incident solar radiation that is
#'   direct (0-1; 0 = all diffuse radiation)
#' @param wood_reflect Wood reflectance spectrum. Default [wood_spectrum()]
#' @return 
#' @author Alexey Shiklomanov
#' @export
edr_r <- function(pft, lai, wai, cai,
                  N, Cab, Car, Cw, Cm,
                  orient_factor, clumping_factor,
                  soil_moisture,
                  direct_sky_frac,
                  czen,
                  wood_reflect = matrix(rep(wood_spec, length(pft)), 2101),
                  wavelengths = seq(400, 2500)) {
  ncohort <- length(pft)
  npft <- length(N)
  nwl <- length(wavelengths)
  stopifnot(
    length(pft) == ncohort,
    length(lai) == ncohort, all(lai >= 0),
    length(wai) == ncohort, all(wai >= 0),
    length(cai) == ncohort,
    all(cai >= 0), all(cai <= 1),
    length(N) == npft, all(N >= 1),
    length(Cab) == npft, all(Cab > 0),
    length(Car) == npft, all(Car > 0),
    length(Cw) == npft, all(Cw > 0),
    length(Cm) == npft, all(Cm > 0),
    length(orient_factor) == npft,
    all(orient_factor < 1), all(orient_factor > -1),
    length(clumping_factor) == npft,
    all(clumping_factor > 0), all(clumping_factor <= 1),
    length(soil_moisture) == 1,
    soil_moisture <= 1, soil_moisture >= 0,
    length(direct_sky_frac) == 1,
    direct_sky_frac >= 0, direct_sky_frac <= 1,
    length(czen) == 1,
    NROW(wood_reflect) %in% c(2101, nwl)
  )

  # Wavelength indices -- everything relative to 400:2500 (so 400nm is
  # index 1, 2500 is index 2101)
  wli <- wavelengths - 399

  # If using full wood reflectance spectrum, subset to only used
  # wavelengths
  if (nwl != NROW(wood_reflect)) wood_reflect <- wood_reflect[wli, ]

  leaf_spectra <- Map(rrtm::prospect5, N, Cab, Car, Cw, Cm)
  leaf_reflect <- Reduce(
    cbind,
    Map(function(x) x[["reflectance"]], leaf_spectra)
  )
  leaf_reflect <- leaf_reflect[wli,]
  leaf_trans <- Reduce(
    cbind,
    Map(function(x) x[["transmittance"]], leaf_spectra)
  )
  leaf_trans <- leaf_trans[wli,]

  # Soil reflectance as a function of soil moisture
  soil_reflect <- hapke_soil(soil_moisture)[wli]

  # "Flat" spectra of incident solar radiation
  down0_sky <- rep(direct_sky_frac, nwl)
  down_sky <- rep(1 - direct_sky_frac, nwl)

  # Wood does not transmit in the VIS or NIR
  wood_trans <- wood_reflect
  wood_trans[] <- 0

  sw_two_stream(
    czen = czen,
    iota_g = soil_reflect,
    pft = pft,
    lai = lai,
    wai = wai,
    cai = cai,
    orient_factor = orient_factor,
    clumping_factor = clumping_factor,
    leaf_reflect = leaf_reflect,
    leaf_trans = leaf_trans,
    wood_reflect = wood_reflect,
    wood_trans = wood_trans,
    down_sky = down_sky,
    down0_sky = down0_sky,
    wavelengths = wavelengths
  )
}
