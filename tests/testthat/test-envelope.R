#to test the envelope() function
library(testthat)
context('envelope returns the correct u_k(x)')

#Test the regular cases
test_that('envelope() returns the correct output', {

  z <- c(0,1,2,3)
  abscissae<- c(.5,1.5,2.5)
  x<- 0.4
  hx <- c(10,20,8)
  dhx <- c(1,.8,.6)
  expect_equal(envelope(z, abscissae, x,hx,dhx), 9.9)
})

test_that('envelope() returns the correct output', {

  z <- c(0,1,2,3)
  abscissae <- c(.5,1.5,2.5)
  x <- 2.6
  hx <- c(10,20,8)
  dhx <- c(1,.8,.6)
  expect_equal(envelope(z, abscissae, x, hx, dhx) , 8.06)
})

#Test the corner case where 'x' equals to the upper bound and lower bound
test_that('envelope() returns the correct output', {

  z <- c(0,1,2,3)
  abscissae<- c(.5,1.5,2.5)
  x<- 3
  hx <- c(10,20,8)
  dhx <- c(1,.8,.6)
  expect_equal(envelope(z, abscissae, x,hx,dhx), 8.3)
})

test_that('envelope() returns the correct output', {

  z <- c(0,1,2,3)
  abscissae<- c(.5,1.5,2.5)
  x<- 0
  hx <- c(10,20,8)
  dhx <- c(1,.8,.6)
  expect_equal(envelope(z, abscissae, x,hx,dhx), 9.5)
})

