context('normal cases Mainline fuction tests')
install.packages('NormalLaplace')
library(NormalLaplace)

ddgamma=function(x)
{
  return()
}
ddbeta=function(x)
{
  return()
}
ddchisqu=function(x)
{

}
ddextremevalue=function(x)
{

}
ddlaplace=function(x)
{

}
ddweibull=function(x)
{

}
test_that('Mainline() returns the correct output when given distribution is normal(0,1) distribution', {
 set.seed(0)
  expect_equal(Mainline(n=1000,g=dnorm,initialPoints=c(1,3,5),leftbound=0,rightbound=6)[1:5],c(1.69148297,0.41329507,1.80531251,0.95924730,0.07400522))
})
test_that('ainline() send the correct information when given distribution is normal(0,1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dnorm,initialPoints=c(1,3,120),leftbound=0,rightbound=121), 'derivates of initial abscissae could not be evaluted numerically.')
})


test_that('Mainline() returns the correct output when given distribution is uniform(0,1) distribution', {
  set.seed(0)
  expect_true(all.equal(Mainline(n=1000,g=dunif,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.8966972,0.2465487,0.8898493,0.8780676,0.5929574),tolerance=1e-07))
})
test_that('ainline() send the correct information when given distribution is uniform(0,1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dunif,initialPoints=c(1,3,120),leftbound=0,rightbound=121), 'derivates of initial abscissae could not be evaluted numerically.')
})


test_that('Mainline() returns the correct output when given distribution is  logistic distribution(location=0,scale=1)',{
  set.seed(0)
  expect_true(all.equal(Mainline(n=1000,g=dlogis,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.8897672,0.3486874,0.8943206,0.8832789,0.6307263),tolerance=1e-07))
})
test_that('ainline() send the correct information when given distribution is logistic distribution(location=0,scale=1) with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dlogis,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

#need to discuss how to import x

test_that('Mainline() returns the correct output when given distribution is  gamma (shape=3,rate=1) distribution', {
  set.seed(0)
  expect_true(all.equal(Mainline(n=1000,g=dgamma(x=1,shape=3),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
})
test_that('ainline() send the correct information when given distribution is (shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dgamma(x=1,shape=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

test_that('Mainline() returns the correct output when given distribution is  beta (shape=3,rate=1) distribution', {
  set.seed(0)
  expect_true(all.equal(Mainline(n=1000,g=dgamma(x=1,shape=4),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
})
test_that('ainline() send the correct information when given distribution is beta(shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dgamma(x=1,shape1=2,shape2=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})
#add chi-squre
test_that('Mainline() returns the correct output when given distribution is  beta (shape=3,rate=1) distribution', {
  set.seed(0)
  expect_true(all.equal(Mainline(n=1000,g=dbeta(x=1,,shape1=4,shape2=3),initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8822070,0.3266133 ,0.8771935, 0.8648125,0.6002529),tolerance=1e-07))
})
test_that('ainline() send the correct information when given distribution is beta(shape=3,rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(Mainline(n=1000,g=dgamma(x=1,shape1=2,shape2=3),initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

#add extrme value

#add laplace

#add weibull
#

dgamma
gamma(3)
