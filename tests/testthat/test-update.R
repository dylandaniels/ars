context('test updateDisVals function')

test_that('test the derivative vector is not non-increasing',{
  abscissae <- c(6,7,8,9,10)
  hx <- c(1,2,3,2,1)
  dhx <- c(-9,-8,-7,-6,-5)
  xStar <- 8.1
  hxStar <- 2.7
  dhxStar <--1
  expect_error(updateDistVals(abscissae,hx,dhx,xStar,hxStar,dhxStar),
               'The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold')
})

test_that('test if xstar is in the middle of abscissae',{
  abscissae <- c(6,7,8,9,10)
  hx <- c(1,2,3,2,1)
  dhx <- c(6,1,0,-1,-2)
  xStar <- 8.1
  hxStar <- 2.7
  dhxStar <--0.5
  expect_equal(updateDistVals(abscissae, hx, dhx, xStar, hxStar, dhxStar)$hx,c(1,2,3,2.7,2,1))
})

test_that('test if xstar is in the beginning of abscissae',{
  abscissae <- c(1,2,3,4,5)
  hx <- c(6,7,8,7,6)
  dhx <- c(3,2,-1,-2,-3)
  xStar <- 0.2
  hxStar <- 5
  dhxStar <-4
  expect_equal(updateDistVals(abscissae, hx, dhx, xStar, hxStar, dhxStar)$hx,c(5,6,7,8,7,6))

})

test_that('test if xstar is in the end of abscissae',{
  abscissae <- c(3,4,5,6,7)
  hx <- c(6,7,8,7,6)
  dhx <- c(3,2,-1,-2,-3)
  xStar <- 7.2
  hxStar <- 5.9
  dhxStar <--3.2
  expect_equal(updateDistVals(abscissae, hx, dhx, xStar, hxStar, dhxStar)$hx,c(6,7,8,7,6,5.9))

})
