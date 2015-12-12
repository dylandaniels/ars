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


#Returns TRUE if the vec[j+1] <= vec[j] for all j=1...length(vec)-1
isDecreasing <- function (vec)
{
  k <- length(vec)

  #The differences vec[j+1] - v[j] <= 0 for all j=1...k-1
  diff <- vec[2:k] - vec[1:k-1]

  #Count the number of times the difference is <= 0
  leq <- sum ( diff <= 10e-10 )

  #If all the differences are <= 0, then the vector is decreasing
  if ( leq == k-1 )
    return(TRUE)
  else #Otherwise, the vector is not decreasing
    return (FALSE)
}


#abscissae is the set of points: T_k = { x_1 , x_2, ... , x_k } (assumption: it is sorted)
#distFun is g(x), (it is a function object)
#derivFun is an optional derivative function, i.e. g'(x)
#If no proof against log-concavity is found, the function returns
#a list of intersecting points z_j, j=1, ... , k-1
#Note that the returned vector may be smaller than length k-1 if there are
#abscissae points x_j where h'(x_j) == h'(x_{j+1})
envelopeIntersectPoints <- function ( abscissae, hx, dhx )
{
  #Check if dhx is decreasing (log-concavity requirement)
  if( isDecreasing(dhx) )
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
    #and there are no intersecting points, the function returns x_1, x_k

    #Initialize the intersects vector
    intersects <- rep(0, k-1)
    if ( length(zeros) == k-1 )
    {
      intersects <- c( rep(min(abscissae), k-2), max(abscissae) )
    }
    else if (length(zeros) != 0) #If some of the differences are equal to zero
    {
      #For the points where the denominator is not equal to zero
      intersects[-zeros] <- numerator[-zeros]/denominator[-zeros]

      for (i in 1:length(zeros))
      {
        index <- zeros[i]
        if (index != 1)
          intersects[index] <- intersects[index-1]
        else
          intersects[index] <- min(abscissae)
      }
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
    return (NULL)
  }
}

updateDistVals <- function(abscissae, hx, dhx, xStar, hxStar, dhxStar)
{
  leq <- (abscissae <= xStar)
  index <- sum(leq)

  k <- length(abscissae)
  newHx <- rep(0, k+1)
  newDhx <- rep(0, k+1)

  newHx <- c(hx[1:index], hxStar , hx[(index+1):k])
  newDhx <- c(dhx[1:index], dhxStar , dhx[(index+1):k])

  #Check if the dhx vector is still decreasing
  # TODO: correct corner case (what if abscissae added at beginning or end?)
  #if ( (newDhx[index+1] - newDhx[index] > 10e-10) || (newDhx[index+2] - newDhx[index+1] > 10e-10) )
  #  warning("In updateDistVals: Log-concavity assumption is violated, the vector of h'(x) is non-decreasing.")


  #Update the abscissae
  newAbs <- rep(0, k+1)
  newAbs <- c(abscissae[1:index], xStar, abscissae[(index+1):k])

  return( list(hx = newHx, dhx = newDhx, abscissae = newAbs) )
}

#NOT FINISHED YET
updateIntersects <- function(abscissae, intersects, xStar, hx, dhx)
{
  leq <- (abscissae <= xStar)
  index <- sum(leq) + 1

  k <- length(abscissae)

  newIntersects <- rep(0, k)
  newIntersects[1:(index-1)] <- intersects[1:index-1]

  #Intersection of the tangent lines at x_index and xStar
  xj <- abscissae[index]
  xj1

  newIntersects[index] <- hx[index+1] - hx[index] - xStar*dhx[index+1] + abscissae[index]*dhx[index] -

  #Intersection of the tangent lines at xStar and x_{index+1}
  newIntersects[index + 1] <- smthng
  newIntersects[(index+2):k] <- intersects[(index+1):k-1]
}


#xStar is a point sampled from s_k(x) function
#lowerFun is the l_k(x) function
#upperFun is the u_k(x) function
#logDistFun is h(x) = ln( g(x) ) function
#Returns TRUE if xStar is ACCEPTED
#Returns FALSE if xStar is REJECTED
acceptReject <- function (xStar, lowerFun, upperFun, logDistFun, logDistDerivFun)
{
  #A random uniform(0,1) value
  w <- runif(1)

  myList <- list(dec = FALSE, step = NULL, hx = NULL, dhx = NULL)

  upperVal <- upperFun(xStar)
  if ( w <= exp( lowerFun(xStar) - upperVal ) ) #Squeezing test
  {
    myList$dec <- TRUE
    myList$step <- 1
  }
  else if ( w <= exp ( logDistFun(xStar) - upperVal ) ) #Rejection test
  {
    myList$dec <- TRUE
    myList$step <- 2
    myList$hx <- logDistFun(xStar)
    myList$dhx <- logDistDerivFun(xStar)
  }
  else
  {
    myList$step <- 2
    myList$hx <- logDistFun(xStar)
    myList$dhx <- logDistDerivFun(xStar)
  }

  return ( myList )
}
