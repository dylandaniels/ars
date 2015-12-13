# This is equation (2) on page 340
# x is the point to evaluate
# z is the vector of envelop intersection points
# sampleFunc is h(x)
# derivSampleFunc is h'(x)

#take the input vector of z(sorted), one value of xstar, vector of h(x) and vector of h'(x)
#and output the value of u_k(xstar)
# TODO: rename x to abscissae
# TODO: rename xstar to x
envelope <- function(z, abscissae, x, hx, dhx) {
  #print(paste0('x=',x))
  #print(paste0('z=',z))

  if (x == z[length(z)]) {
    index = length(z) - 1
  } else {
    index <- sum(z <= x)
  }
  ux <- hx[index] + (x - abscissae[index]) * dhx[index]
  return(ux)
}

#to test the function
##testtaht
# context('envelope returns the correct u_k(xstar)')
#
# test_that('envelope() returns the correct output', {
#
#   z <- c(1,2,3)
#   x <- c(.5,1.5,2.5)
#   xstar <- 0.4
#   sampleFunc <- c(10,20,8)
#   derivX <- c(1,.8,.6)
#   expect_equal(envelope(z, x,xstar,sampleFunc,derivX), 9.9)
# })
#
# test_that('envelope() returns the correct output', {
#
#   z <- c(1,2,3)
#   x <- c(.5,1.5,2.5)
#   xstar <- 2.6
#   sampleFunc <- c(10,20,8)
#   derivX <- c(1,.8,.6)
#   expect_equal(envelope(z, x,xstar,sampleFunc,derivX), 8.06)
# })




