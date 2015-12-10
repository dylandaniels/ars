library(stats)
library(numDeriv)

#Some trial functions

f <- function(x, lambda = 1)
{
  if ( x >= 0 )
  {
    return( lambda*exp(-lambda*x) )
  }
  else
  {
    return( 0 )
  }
}


df <- function(x, lambda = 1)
{
  if ( x >= 0 )
  {
    return( -(lambda^2)*exp(-lambda*x) )
  }
  else
  {
    return( 0 )
  }
}

getC <- function(funcName, min = -Inf, max = Inf)
{
  c <- tryCatch(integrate(Vectorize(as.function(funcName)), min, max),
                error = function(e) {
                  print("Error in integration")
                  return (0) },
                warning = function(w) {
                  print("Warning in integration")
                  return(0)  })
  if ( class(c) == "numeric" )
  {
    return ( 0 )
  }
  else
  {
    return (c$value)
  }
}

#Given a function f, convertToLog returns
#the function object of ln f(x)
convertToLog <- function(f)
{
  lnf <- function(x)
  {
    return ( log( f(x) ) )
  }
  return ( lnf )
}


#Given a function f and df/dx, this function
#returns the function object of d(ln f)/dx
convertDerivToLog <- function(f, derivF)
{
  derivLogF <- function(x)
  {
    return( (1/f(x))*derivF(x) )
  }
  return (derivLogF)
}

#abscissae is the set of points: T_k = { x_1 , x_2, ... , x_k }
#The function sorts the points in increasing order
#distFun is g(x), (it is a function object)
#derivFun is an optional derivative function, i.e. g'(x)
#If no proof against log-concavity is found, the function returns
#a list of intersecting points z_j, j=1, ... , k-1
#Note that the returned vector may be smaller than length k-1 if there are
#abscissae points x_j where h'(x_j) == h'(x_{j+1})
envelopeIntersectPoints <- function ( abscissae, distFun, derivFun = NULL )
{
  #Sort the x_j points in increasing order:
  abscissae <- sort(abscissae)

  #Take the natural logarithm of the distribution function g(x)
  #lnDistFun is our h(x)
  logDistFun <- convertToLog(distFun)

  #Evaluate h(x) at the points in abscissae
  hx <- Vectorize(logDistFun)(abscissae)

  #If the user did not provide a derivative function
  if ( is.null(derivFun) )
  {
    #Evaluate h'(x) at the points in abscissae
    dhx <- grad(Vectorize(logDistFun), abscissae)
  }
  else
  {
    derivLogFun <- convertDerivToLog(distFun, derivFun)
    dhx <- Vectorize(derivLogFun)(abscissae)
  }

  #Check if dhx is decreasing (log-concavity requirement)
  if( sum( abs( dhx - sort(dhx, decreasing = TRUE) ) ) <= 10e-10 )
  {
    #Number of points in abscissae
    k <- length(abscissae)

    #Numerator of equation (1)
    numerator <- hx[2:k] - hx[1:k-1] - abscissae[2:k]*dhx[2:k] + abscissae[1:k-1]*dhx[1:k-1]

    #Denominator of equation (1)
    denominator <- dhx[1:k-1]-dhx[2:k]

    #For points x_j and x_{j+1}, is the derivative the same?
    zeros <- which(denominator <= 10e-10)

    #If all the differences are zero, then h(x) must be linear on the abscissae,
    #and there are no intersecting points, the function returns x_1
    if ( length(zeros) == k-1 )
    {
      intersects <- c(min(abscissae), max(abscissae))
    }
    else if (length(zeros) != 0) #If some of the differences are equal to zero
    {
      intersects <- numerator[-zeros]/denominator[-zeros]
    }
    else #If none of the differences are equal to zero.
    {
      intersects <- numerator/denominator
    }
    return (intersects)
  }
  else #If the h'(x) is not decreasing
  {
    warning ("The sampling function is not log-concave.")
    print(dhx)
    print(sort(dhx, decreasing = TRUE) )
    print(sum( abs( dhx - sort(dhx, decreasing = TRUE) ) ))
    return (NULL)
  }
}

#xStar is a point sampled from s_k(x) function
#lowerFun is the l_k(x) function
#upperFun is the u_k(x) function
#logDistFun is h(x) = ln( g(x) ) function
#Returns TRUE if xStar is ACCEPTED
#Returns FALSE if xStar is REJECTED
acceptReject <- function (xStar, lowerFun, upperFun, logDistFun)
{
  #A random uniform(0,1) value
  w <- runif(1)
  accepted <- FALSE

  if ( w <= exp( lowerFun(xStar) - upperFun(xStar) ) ) #Squeezing test
  {
    accepted <- TRUE
  }
  else if ( w <= exp ( logDistFun(xStar) - upperFun(xStar) ) ) #Rejection test
  {
    accepted <- TRUE
  }

  return (accepted)
}
