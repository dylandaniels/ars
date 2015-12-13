
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
    stop("The sampling function is not log-concave.")
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
  newAbs <- rep(0, k+1)

  if (index != 0 && index != k)
  {
    #Update the hx vector
    newHx <- c(hx[1:index], hxStar , hx[(index+1):k])

    #Update the dhx vector
    newDhx <- c(dhx[1:index], dhxStar , dhx[(index+1):k])

    #Update the abscissae
    newAbs <- c(abscissae[1:index], xStar, abscissae[(index+1):k])

    #Check if the dhx vector is still decreasing
    if ( (newDhx[index+1] - newDhx[index] > 10e-10) || (newDhx[index+2] - newDhx[index+1] > 10e-10) )
      stop("In updateDistVals: Log-concavity assumption is violated, the vector of h'(x) is non-decreasing.")
  }
  else if (index == 0)
  {
    #Update the hx vector
    newHx <- c(hxStar , hx[(index+1):k])

    #Update the dhx vector
    newDhx <- c(dhxStar , dhx[(index+1):k])

    #Update the abscissae
    newAbs <- c(xStar, abscissae[(index+1):k])

    #Check if the dhx vector is still decreasing
    if ( newDhx[index+2] - newDhx[index+1] > 10e-10 )
      stop("In updateDistVals: Log-concavity assumption is violated, the vector of h'(x) is non-decreasing.")
  }
  else #If the index == k
  {
    #Update the hx vector
    newHx <- c(hx[1:index], hxStar)

    #Update the dhx vector
    newDhx <- c(dhx[1:index], dhxStar)

    #Update the abscissae
    newAbs <- c(abscissae[1:index], xStar)

    #Check if the dhx vector is still decreasing
    if ( newDhx[index+1] - newDhx[index] > 10e-10 )
      stop("In updateDistVals: Log-concavity assumption is violated, the vector of h'(x) is non-decreasing.")
  }

  return( list(hx = newHx, dhx = newDhx, abscissae = newAbs) )
}

updateIntersects <- function(abscissae, oldIntersects, hx, dhx, xStar, hxStar, dhxStar)
{
  leq <- (abscissae <= xStar)
  index <- sum(leq)

  k <- length(abscissae)
  intersects <- oldIntersects[2:k]

  #cat("\n\n")
  #print(paste0("LENGTH OF ABS: ", k))
  #print(paste0("LENGTH OF OLD INTERSECT: ", length(oldIntersects)))
  #print(paste0("INTERSECTS: ", intersects))
  #cat("\n\n")

  #Drop -Inf and Inf from two ends of the intersects
  #intersects <- oldIntersects[2:k]

  newIntersects <- rep(0, k)

  #Only if index >= 2, we copy the first section of intersects
  if (index > 1)
    newIntersects[1:(index-1)] <- intersects[1:(index-1)]

  #Only if index <= k-2, we copy the last section of intersects
  if (index + 2 <= k )
    newIntersects[(index+2):k] <- intersects[(index+1):(k-1)]

  #Intersection of the tangent lines at x_index and xStar
  if (index != 0 && index != k)
  {
    #j = index, we will need to update intersects z_index and z_{index+1}
    xj <- abscissae[index]
    xj1 <- abscissae[index+1]

    #Intersect of x_j and xStar

    if ((dhx[index] - dhxStar) > 10e-10) #If the tangent values are different
      newIntersects[index] <- (hxStar - hx[index] - xStar*dhxStar + xj*dhx[index])/(dhx[index] - dhxStar)
    else if (index == 1) #If index = 1 then this is the first intersect point
      newIntersects[index] <- min(abscissae)
    else
      newIntersects[index] <- newIntersects[index - 1]

    #Intersection of the tangent lines at xStar and x_{j+1}
    if ((dhxStar - dhx[index+1]) > 10e-10) #If the tangent values are different
      newIntersects[index + 1] <- (hx[index+1] - hxStar - xj1*dhx[index+1] + xStar*dhxStar)/(dhxStar - dhx[index+1])
    else
      newIntersects[index] <- newIntersects[index-1]
  }
  else if (index == 0)
  {
    xj1 <- abscissae[index+1]

    #Intersection of the tangent lines at xStar and x_{j+1}
    if ((dhxStar - dhx[index+1]) > 10e-10) #If the tangent values are different
      newIntersects[index + 1] <- (hx[index+1] - hxStar - xj1*dhx[index+1] + xStar*dhxStar)/(dhxStar - dhx[index+1])
    else
      newIntersects[index + 1] <- min(abscissae)
  }
  else #index == k
  {
    xj <- abscissae[index]

    #Intersect of x_j and xStar
    if ((dhx[index] - dhxStar) > 10e-10) #If the tangent values are different
      newIntersects[index] <- (hxStar - hx[index] - xStar*dhxStar + xj*dhx[index])/(dhx[index] - dhxStar)
    else
      newIntersects[index] <- newIntersects[index-1]
  }

  return(newIntersects)
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
