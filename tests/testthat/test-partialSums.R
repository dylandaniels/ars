context('partial sums for s_k(x) integral')

test_that('partialSums() returns the correct output', {
  z <- seq(3) # mock intersection points
  u <- function (x) { log(x) } # mock envelope function
  dhx <- c(2, 3)  # mock h'(x) evaluated at abscissae
  expect_equal(partialSums(z, u, dhx), cumsum(c(1/2, 1/3)))
})
