library(australian)
context("Test Testing")

test_that("1 is equal to 1", {
  expect_equal(1, 1)
})

test_that("copula is good", {
  expect_equal(copula(2), 6)
})
