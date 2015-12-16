context('normal cases ars fuction tests')

#install.packages('NormalLaplace')
#library(NormalLaplace)

# ddgamma=function(x)
# {
#   return()
# }
# ddbeta=function(x)
# {
#   return()
# }
# ddchisqu=function(x)
# {
#
# }
# ddextremevalue=function(x)
# {
#
# }
# ddlaplace=function(x)
# {
#
# }
# ddweibull=function(x)
# {
#
# }

test_that('ars() returns the correct output when given distribution is normal(0,1) distribution', {
  set.seed(0)
  expect_equal(ars(n=1000,g=dnorm,initialPoints=c(1,3,5),leftbound=0,rightbound=6)[1:5],c(1.69148297,0.41329507,1.80531251,0.95924730,0.07400522))
})

test_that('ainline() send the correct information when given distribution is normal(0,1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=dnorm,initialPoints=c(1,3,120),leftbound=0,rightbound=121), 'derivates of initial abscissae could not be evaluted numerically.')
})


test_that('ars() returns the correct output when given distribution is uniform(0,1) distribution', {
  set.seed(0)
  expect_true(all.equal(ars(n=1000,g=dunif,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.8966972, 0.3721239, 0.9082078, 0.8983897, 0.6607978),tolerance=1e-07))
})

test_that('ainline() send the correct information when given distribution is uniform(0,1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=dunif,initialPoints=c(1,3,120),leftbound=0,rightbound=121), 'derivates of initial abscissae could not be evaluted numerically.')
})


test_that('ars() returns the correct output when given distribution is  logistic distribution(location=0,scale=1)',{
  set.seed(0)
  expect_true(all.equal(ars(n=1000,g=dlogis,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.8897672,0.3486874,0.8943206,0.8832789,0.6307263),tolerance=1e-07))
})

test_that('ainline() send the correct information when given distribution is logistic distribution(location=0,scale=1) with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=dlogis,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

#need to discuss how to import x

# FIX ALL THESE TESTS
# test_that('ars() returns the correct output when given distribution is  gamma (shape=3,rate=1) distribution', {
#   set.seed(0)
#   expect_true(all.equal(ars(n=1000,g=dgamma(x=1,shape=3),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
# })

# test_that('ars() sends the correct information when given distribution is (shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
#   set.seed(0)
#   expect_error(ars(n=1000,g=dgamma(x=1,shape=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
# })
#
# test_that('ars() returns the correct output when given distribution is  beta (shape=3,rate=1) distribution', {
#   set.seed(0)
#   expect_true(all.equal(ars(n=1000,g=dgamma(x=1,shape=4),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
# })
#
# test_that('ainline() send the correct information when given distribution is beta(shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
#   set.seed(0)
#   expect_error(ars(n=1000,g=dgamma(x=1,shape1=2,shape2=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
# })
#
# #add chi-squre
# test_that('ars() returns the correct output when given distribution is  beta (shape=3,rate=1) distribution', {
#   set.seed(0)
#   expect_true(all.equal(ars(n=1000,g=dbeta(x=1,,shape1=4,shape2=3),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
# })
#
# test_that('ainline() send the correct information when given distribution is beta(shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
#   set.seed(0)
#   expect_error(ars(n=1000,g=dgamma(x=1,shape1=2,shape2=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
# })

#add extrme value

#add laplace

#add weibull
#

#dgamma
#gamma(3)

#context('test the ars function for abnormal input functions')

test_that('test the ars function with piecewise input functions but continuous and differentiable everywhere', {
  set.seed(123)
  myfun <- Vectorize(function (x) {
    if (x >= 0 && x <= 1) {
      return(sqrt(1 - (x - 1)^2))
    } else {
      return(1)
    }
  })
  expect_equal(ars(10,myfun,initialPoints =c(1,2),leftbound=0,rightbound=5)[1:3],
               c(1.4378876006230712, 2.0448846090584993, 4.7023364214692265))
})

test_that('test the ars function with piecewise input functions but some points not differentiable', {
  set.seed(123)
  myfun<- Vectorize(function (x) {
    if (x >= 0 && x <= 1) {
      return(x^2)
    } else {
      return(1)
    }
  })
  expect_equal(ars(10,myfun,initialPoints = c(1,2),leftbound=0,rightbound=5)[1:3],
               c(1.4378876006230712, 2.0448846090584993, 4.7023364214692265))
})

test_that('test the ars functions when input functions are not continuous somewhere', {
  set.seed(123)
  myfun <- Vectorize(function(x){
    if (x>=0 && x<=3){
      return(x)
    }else if (x>3 && x<=5){
      return(5)
    }else if(x>5){
      return(-x+7)
    }
  })
  expect_equal(ars(10,myfun,initialPoints = c(0.5,1),leftbound = 0,rightbound = 6)[1:3],
               c(3.459433, 5.744093, 3.447403),
               tol=1e-07
  )
})

test_that('test the ars functions when input functions are not log-concave', {

  myfun <- Vectorize(function(x){
    return((x-5)^2)
  })
  expect_error(ars(100,myfun,initialPoints = c(4,6),leftbound = 0,rightbound = 10),
               "The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold."
  )
})
