context("Priors work")

p <- create_prior()

for (i in 1:500) {
  z <- p$density(p$sampler())
  stopifnot(is.finite(z))
}
