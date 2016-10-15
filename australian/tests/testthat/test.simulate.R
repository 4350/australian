library(australian)
library(foreach)
library(mvtnorm)

context("Test Simulation")

copula_spec <- function(dynamics, distribution) {
  CopulaSpecification(
    dynamics = do.call(CopulaDynamics, dynamics),
    distribution = do.call(CopulaDistribution, distribution)
  )
}

test_that('Simulate(constant, normal)', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0, beta = 0, Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  registerDoSEQ()
  simulations <- copula_simulate(spec, 100, 50)
  expect_length(simulations, 50)
  expect_equal(dim(simulations[[1]]), c(100, 2))
})

test_that('Simulate(dynamic, normal)', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0.06, beta = 0.91, Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  registerDoSEQ()
  simulations <- copula_simulate(spec, 100, 5)
  expect_length(simulations, 5)
  expect_equal(dim(simulations[[1]]), c(100, 2))
})
