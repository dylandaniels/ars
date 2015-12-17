context('normal cases ars fuction tests')



 ddgamma=Vectorize(function(x,alpha=2,beta=3)
 {
   return(beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x))
 })
 ddbeta=Vectorize(function(x,alpha=2,beta=3)
 {
   return(x^(alpha-1)*(1-x)^(beta-1)/beta(alpha,beta))
 })
 ddchisqu=Vectorize(function(x,k=6)
 {
  return(1/(2^(k/2)*gamma(k/2))*x^(k/2-1)*exp(-x/2))
 })
 ddextremevalue=Vectorize(function(x,ep=-5,sig=5,mu=1)
 {
  return((1+(x-mu)/sig*ep)^(-1/ep))
 })
 ddlaplace=Vectorize(function(x,mu=2,b=3)
 {
   if(x>=mu)
  {return (1/(2*b)*exp((mu-x)/b))}
  if (x<mu)
{return (1/(2*b)*exp((x-mu)/b))}
 })
 ddweibull=Vectorize(function(x,lambda=2,k=3)
 {
  return((k/lambda)*(x/lambda)^(k-1)*exp(-(x/lambda)^k))
 })

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




 test_that('ars() returns the correct output when given distribution is  gamma (alpha=2,beta=3) distribution', {
   set.seed(0)
   expect_true(all.equal(ars(n=1000,g=ddgamma,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.9163975,0.3660236,0.8510459,0.8373588,0.5820760),tolerance=1e-07))
 })

 test_that('ars() sends the correct information when given distribution is (alpha=2,beta=3) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
   set.seed(0)
   expect_error(ars(n=1000,g=ddgamma,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
 })

 test_that('ars() returns the correct output when given distribution is  beta (alpha=2,beta=3) distribution', {
   set.seed(0)
   expect_true(all.equal(ars(n=1000,g=ddbeta,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(  0.38658874,0.76139262, 0.48610130, 0.09888136, 0.19441745),tolerance=1e-07))
 })

 test_that('ainline() send the correct information when given distribution is beta(alpha=2,beta=3) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
   set.seed(0)
   expect_error(ars(n=1000,g=ddbeta,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
 })


 test_that('ars() returns the correct output when given distribution is  chisqu (df=6) distribution', {
   set.seed(0)
   expect_true(all.equal(ars(n=1000,g=ddchisqu,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.6693970,0.9609523,0.9566337,0.8423037,0.3733350),tolerance=1e-07))
 })

 test_that('ainline() send the correct information when given distribution is chisqu(df=6) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
   set.seed(0)
   expect_error(ars(n=1000,g=ddchisqu,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
 })

test_that('ars() returns the correct output when given distribution is extreme value(epsion=-5,sigma=5,mu=1) distribution', {
  set.seed(0)
  expect_true(all.equal(ars(n=1000,g=ddextremevalue,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c(0.8911069,0.3572684,0.9016454 ,0.8912197,0.6449157),tolerance=1e-07))
})

test_that('ainline() send the correct information when given distribution is extreme value(epsion=-5,sigma=5,mu=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=ddextremevalue,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

test_that('ars() returns the correct output when given distribution is weibull(mu=2,b=3) distribution', {
  set.seed(0)
  expect_true(all.equal(ars(n=1000,g=ddweibull,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.6838474, 0.9629708, 0.9588711, 0.8499925, 0.3962677),tolerance=1e-07))
})

test_that('ainline() send the correct information when given distribution is weibull(mu=2,b=3) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=ddweibull,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})

test_that('ars() returns the correct output when given distribution is exponential(rate=1) distribution', {
  set.seed(0)
  expect_true(all.equal(ars(n=1000,g=dexp,initialPoints=c(0.1,0.2,0.3),leftbound=0,rightbound=1)[1:5],c( 0.8366036,0.2681764,0.8535432,0.8390765,0.5407761),tolerance=1e-07))
})

test_that('ainline() send the correct information when given distribution is exponential(rate=1) distribution with the case derivates of initial abscissae could not be evaluted numerically.',{
  set.seed(0)
  expect_error(ars(n=1000,g=dexp,initialPoints=c(1,3,12000),leftbound=0,rightbound=12100), 'derivates of initial abscissae could not be evaluted numerically.')
})





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

test_that('test the ars functions when initial points are not provided',{
  library(optimx)
  set.seed(123)
  expect_equal(ars(100,dnorm,leftbound = -10,rightbound=10)[1:3],c(-0.64761367424818173,-0.27305561111952537,2.23018308774059637))
})

test_that('test the ars functions when derivtive functions are given',{
  library(optimx)
  set.seed(123)
  derivX <- function(x){
    if(x>=0 && x<=1)return(1)
    if(x>1 && x<=2)return(-1)
  }
  myfun <- function(x){
    if (x>=0 && x<=1) return(x)
    if (x>1 && x<=2) return(2-x)
  }
  expect_equal(ars(100,myfun,dg = derivX,leftbound = 0,rightbound = 2)[1:3],c(1.11236707972281512,0.93374402364032161,1.81458526158709388))
})
