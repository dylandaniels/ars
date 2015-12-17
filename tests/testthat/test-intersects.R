context('test the envelopIntersecPoints function')

test_that('test the log-concavity of the function',{
  abscissae <- sort(runif(10))
  hx <- c(1:5,4:1)
  dhx <- sort(runif(10))
  expect_error(envelopeIntersectPoints(abscissae,hx,dhx),
               'The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold')
}
)

test_that('test when h(x)s are linear on abscissae',{
  set.seed(123)
  abscissae <- sort(runif(10))
  hx <- rep(1,10)
  dhx <- rep(0,10)
  expect_equal(envelopeIntersectPoints(abscissae,hx,dhx)[1:3],c(0.0455565,0.0455565,0.0455565))
})

test_that('test when some of dhx are same',{
  set.seed(123)
  abscissae <- sort(runif(10))
  hx <- c(rep(1,4),sort(runif(6),decreasing = TRUE))
  dhx <- c(rep(0,4),sort(runif(6,min=-1,max = 0),decreasing = TRUE))
  expect_equal(envelopeIntersectPoints(abscissae,hx,dhx)[1:3],
               c(0.16656700975727290,0.34827722096815705,0.43279582855757326))
})

test_that('test when all the dhx are different',{
  set.seed(123)
  abscissae <- sort(runif(10))
  hx <- c(1:5,5:1)
  dhx <- c(sort(runif(10,min=-1,max = 0),decreasing = TRUE))
  expect_equal(envelopeIntersectPoints(abscissae,hx,dhx)[1:3],
               c(434.0124712081788516,18.7986537180344584,4.9774359809600703))
})

context('test the updateIntersects fucntion')

test_that('test if xstar is inserted in the middle',{

  abscissae <- c(1,2,3,4,5)
  oldIntersects <- c(0,1.5,2.5,3.5,4.5,5.5)
  hx <- c(1,2,3,2,1)
  dhx <- c(5,4,3,2,1)
  xStar <- 2.7
  hxStar <- 2.8
  dhxStar <- 3.5
  expect_equal(updateIntersects(abscissae, oldIntersects, hx, dhx, xStar, hxStar, dhxStar)[1:3],c( 1.50,-1.30,1.30))
})

test_that('test if xstar is inserted in the beginning',{
  abscissae <- c(1,2,3,4,5)
  oldIntersects <- c(0,1.5,2.5,3.5,4.5,5.5)
  hx <- c(1,2,3,2,1)
  dhx <- c(5,4,3,2,1)
  xStar <- 0.5
  hxStar <- 0.9
  dhxStar <- 5.1
  expect_equal(updateIntersects(abscissae, oldIntersects, hx, dhx, xStar, hxStar, dhxStar)[1:3],c( -23.50,1.50,2.50))

})

test_that('test if xstar is inserted at the end',{
  abscissae <- c(1,2,3,4,5)
  oldIntersects <- c(0,1.5,2.5,3.5,4.5,5.5)
  hx <- c(1,2,3,2,1)
  dhx <- c(5,4,3,2,1)
  xStar <- 5.123
  hxStar <- 0.29
  dhxStar <- 0.38
  expect_equal(updateIntersects(abscissae, oldIntersects, hx, dhx, xStar, hxStar, dhxStar)[1:3],c(1.5,2.5,3.5))


})
