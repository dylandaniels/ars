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

test_that("updateIntegrals() picks correct integrals to update", {
  dummyValue <- 1000
  oldAbscissae <- seq(1, 5)
  oldIntegrals <- rep(dummyValue, 5)
  newZ <- seq(0.5, 5.5)
  u <- function (x) { log(x) }
  dhx <- rep(2, 6)

  # Interior point
  xStar <- 1.5
  # Should update integrals at indicies 1:3
  updatedIntegrals <- updateIntegrals(xStar, oldAbscissae, oldIntegrals, newZ, u, dhx)
  expect_equal(sum(updatedIntegrals[1:3] != dummyValue), 3)
  expect_equal(sum(updatedIntegrals[4:6] == dummyValue), 3)

  # Left border point
  xStar <- 0.5
  # Should update integrals at indicies 1:2
  updatedIntegrals <- updateIntegrals(xStar, oldAbscissae, oldIntegrals, newZ, u, dhx)
  expect_equal(sum(updatedIntegrals[1:2] != dummyValue), 2)
  expect_equal(sum(updatedIntegrals[3:6] == dummyValue), 4)

  # Right border point
  xStar <- 5.2
  # Should update integrals at indicies 5:6
  updatedIntegrals <- updateIntegrals(xStar, oldAbscissae, oldIntegrals, newZ, u, dhx)
  expect_equal(sum(updatedIntegrals[5:6] != dummyValue), 2)
  expect_equal(sum(updatedIntegrals[1:4] == dummyValue), 4)
})
