library(australian)

context("Test Copula Fit")

test_that('copula_fit_theta', {
  r <- copula_fit_theta(2, 'norm', T, F)
  expect_length(r, 3)
  expect_null(r$pars)

  r <- copula_fit_theta(2, 't', T, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('distribution.nu'))
  expect_length(r$ci, 2)
  expect_equal(dim(r$ui), c(2, 1))

  r <- copula_fit_theta(2, 'ghst', T, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('distribution.nu',
                                'distribution.gamma1', 'distribution.gamma2'))
  expect_length(r$ci, 2 + 2 * 2)
  expect_equal(dim(r$ui), c(2 + 2 * 2, 3))

  r <- copula_fit_theta(2, 'norm', F, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('alphabeta.alpha', 'alphabeta.beta'))
  expect_length(r$ci, 3)
  expect_equal(dim(r$ui), c(3, 2))

  r <- copula_fit_theta(2, 'norm', F, T)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('alphabeta.alpha', 'alphabeta.beta',
                                'upsilon.phi', 'upsilon.theta'))
  expect_length(r$ci, 3 + 2)
  expect_equal(dim(r$ui), c(3 + 2, 4))
})
