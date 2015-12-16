context('sampling tests')

test_that("evaluateIntegral() returns the correct output when abs(h'(x)) > 0 ", {
  z <- seq(2) # mock intersection points
  u <- function (x) { log(x) } # mock envelope function
  dhx <- 2 # mock h'(x)
  expect_equal(evaluateIntegral(z[1], z[2], u, dhx), 0.5)
})

test_that("evaluateIntegral() returns the correct output when h'(x) = 0", {
  z <- seq(2) # mock intersection points
  u <- function (x) { log(x) } # mock envelope function
  dhx <- 0 # mock h'(x)
  expect_equal(evaluateIntegral(z[1], z[2], u, dhx), 2)
})
