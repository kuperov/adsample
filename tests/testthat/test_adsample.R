set.seed(123)

test_that("raw ad_sample can generate variates", {
  # use example of the Gamma(3/2, 2) density
  h <- function(tau) {
    y <- log(4*sqrt(2)*sqrt(tau)*exp(-1/2*tau*(3 - 1)^2)/sqrt(pi))
    yprime <- 1/2*(-tau*(3 - 1)^2 + 1)/tau
    c(y, yprime)
  }

  xs <- raw_ad_sample(n = 1e3, log_dens = h,
                      initialPoints = c(0.5, 1),
                      minRange = 0, maxRange = Inf)

  expect_equal(class(xs), 'numeric')
  expect_equal(length(xs), 1e3)

  # compare with true analytical CDF
  f <- function(tau) pgamma(tau, shape = 3/2, rate = 2)
  kt <- ks.test(xs, 'f')
  expect_gt(kt$p.value, 0.01)
})

test_that("raw_ad_sample_debug can generate variates and sensible state variables", {
  set.seed(123)

  # use example of the Gamma(3/2, 2) density
  h.gam <- function(tau) {
    y <- log(4*sqrt(2)*sqrt(tau)*exp(-1/2*tau*(3 - 1)^2)/sqrt(pi))
    yprime <- 1/2*(-tau*(3 - 1)^2 + 1)/tau
    c(y, yprime)
  }

  samp <- raw_ad_sample_debug(n = 1e3, log_dens = h.gam,
                              initialPoints = c(0.5, 1),
                              minRange = 0, maxRange = Inf,
                              maxiter=0)

  expect_equal(length(samp$samples), 1e3)

  expect_equal(length(samp$`T`), samp$k)
  expect_equal(length(samp$Q), samp$k)
  expect_equal(length(samp$H), samp$k)
  expect_equal(length(samp$Hprime), samp$k)
  expect_equal(length(samp$Z), samp$k+1)

  # compare with true analytical CDF
  f <- function(tau) pgamma(tau, shape = 3/2, rate = 2)
  kt <- ks.test(samp$samples, 'f')
  expect_gt(kt$p.value, 0.01)
})

test_that("adsample interface draws samples", {
  set.seed(123)

  h <- function(tau) {
    y <- log(4*sqrt(2)*sqrt(tau)*exp(-1/2*tau*(3 - 1)^2)/sqrt(pi))
    yprime <- 1/2*(-tau*(3 - 1)^2 + 1)/tau
    c(y, yprime)
  }

  xs <- adsample(n = 10, log_dens = h, initialPoints = c(0.5, 1),
                minRange = 0, maxRange = Inf)
  expect_equal(length(xs), 10)
  expect_equal('numeric', class(xs))

  x.debug <- adsample(n = 10, log_dens = h, initialPoints = c(0.5, 1),
                        minRange = 0, maxRange = Inf, debug = TRUE)

  expect_equal('adsample', class(x.debug))
  expect_equal(length(x.debug$samples), 10)
})

test_that("diagnostic plots for adsample objects work without error", {
  # use example of the Gamma(3/2, 2) density
  h.gam <- function(tau) {
    y <- log(4*sqrt(2)*sqrt(tau)*exp(-1/2*tau*(3 - 1)^2)/sqrt(pi))
    yprime <- 1/2*(-tau*(3 - 1)^2 + 1)/tau
    c(y, yprime)
  }
  samp <- adsample(100, h.gam, initialPoints = c(1/2, 1),
                   minRange = 0, maxRange = Inf, debug = TRUE)

  # perform smoke test
  tf <- tempfile(fileext = 'png')
  png(tf)
  plot(samp)
  dev.off()
  unlink(tf)
})

test_that("symbolic differentiation works", {
  # automatic
  auto <- mklogdensf(exp(-1/18*(x-2)^2), x)
  # manual
  man <- function(x) {
    y <- - 1/18*(x-2)^2
    yprime <- -(x - 2)/9
    c(y, yprime)
  }
  xs <- c(-100, -0.1, 0, 2, 8)
  for (x in xs) {
    y.auto <- auto(x)
    y.man  <- man(x)
    expect_equal(y.auto, y.man)
  }
})
