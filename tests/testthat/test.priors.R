context("Priors work")

p <- create_prior()

test_that(
  "Priors don't throw any errors",
  {
    for (i in 1:50) {
      z <- p$density(p$sampler())
      expect_true(is.finite(z))
      expect_true(is.numeric(z))
    }
  }
)
