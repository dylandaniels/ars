# This is equation (2) on page 340
# x is the point to evaluate
# z is the vector of envelop intersection points
# sampleFunc is h(x)
# derivSampleFunc is h'(x)

#take the input vector of z(sorted), one value of xstar, vector of h(x) and vector of h'(x)
#and output the value of u_k(xstar)
envelope <- function(z,x,xstar,sampleFunc,derivX){
  #calculate the vector of absolute difference between z and xstar
  diffs <- abs(xstar-z)
  #calcuate the index of the min of the differences
  index <- which.min(diffs)

  #for the z_index which has the min absolute difference with xstar,
  #either it's the upper bound or lower bound of the interval in which
  #xstar falls
  if (index==1){
    #if xstar is closest to z[1],we have two cases
    if(xstar<=z[1]){
      ustar <- sampleFunc[1]+(xstar-x[1])*derivX[1]
    }else{
      ustar <- sampleFunc[2]+(xstar-x[2])*derivX[2]
    }
  } else{
    if(xstar>=z[index-1] && xstar<=z[index]){
      ustar <- sampleFunc[index]+(xstar-x[index])*derivX[index]
    }else if (xstar>z[index] && xstar<=z[index+1]){
      ustar <- sampleFunc[index+1]+(xstar-x[index+1])*derivX[index+1]
    }

  }


  #return ustar
  return(ustar)
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




