context("Load observations")

import::from(here, here)

test_that(
  "Loading observations works",
  {
    sites <- readLines(here("other_site_data/site_list"))
    obs <- load_observations(sites)
    expect_true(PEcAnRTM::is_spectra(obs))
    expect_equal(PEcAnRTM::wavelengths(obs), seq(400, 1300, by = 10))
  }
)

