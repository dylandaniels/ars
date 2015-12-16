library(optimx)
eps <- 1e-06

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

#If the user does not provide an initial ( or a valid) set of initial points,
#this function returns a set of initial points of size at least 2
findInitPoints <- function(h, leftbound, rightbound)
{
  #Convert h(x) to -h(x)
  convertToNeg <- function(h)
  {
    negh <- function(x)
    {
      return ( -h(x) )
    }
    return ( negh )
  }

  #Convert h(x) into -h(x) because the nlm minimizes the input function
  negh <- convertToNeg(h)

  #Some heuristic to find starting point for the optimization
  lb <- leftbound
  ub <- rightbound
  if( is.infinite(leftbound) )
    lb <- min(0.01, rightbound)
  if( is.infinite(rightbound) )
    ub <- max(0.01, leftbound)

  #If leftbound = -Inf and rightbound = Inf, starting point is 0
  startPt <- (lb + ub)/2

  #Assumption: Left Bound and Right Bound are correctly given!!
  if ( is.infinite(leftbound) && is.infinite(rightbound) )
  {
    #Minimize -h(x) using nlm function, where the max number of iterations is 50
    tryCatch({ sol <- nlm(negh, startPt)
               modeVal <- sol$estimate },
             error = function(e) { stop("Error in nlm for mode calculation.") },
             warning = function(w) { stop(print(w))})
  }
  else #If there are finite bounds at either left or right sides
  {
    tryCatch({ sol <- optimx(startPt, negh, method = "L-BFGS-B", lower = leftbound, upper = rightbound)
              modeVal <- sol[1,1] },
              error = function(e) { stop("Error in optimx for mode calculation.") },
              warning = function(w) { stop(print(w)) })
  }

  #If the mode is not within the domain, then throw an error
  if ( leftbound - eps > modeVal || modeVal > rightbound + eps)
    stop("The mode of the function is not within the domain, make sure the domain is defined correctly.")

  #The initial abscissae of size 3
  initAbs <- numeric(3)
  #The middle abscissae point is the mode
  initAbs[2] <- modeVal

  #If either left or right bound is infinite
  if ( is.infinite(leftbound) || is.infinite(rightbound))
  {
    #If left bound is -infinity
    if( is.infinite(leftbound) )
      initAbs[1] <- modeVal - 1
    else
      initAbs[1] <- (leftbound + modeVal)/10

    #If right bound is infinity
    if (is.infinite(rightbound))
      initAbs[3] <- modeVal + 1
    else
      initAbs[3] <- (rightbound + modeVal)/10
  }
  else #Both bounds are finite
  {
    initAbs[1] <- (leftbound + modeVal)/10
    initAbs[3] <- (rightbound + modeVal)/10
  }

  #We will check for uniqueness here because either the leftbound
  #or the rightbound may be equal to the mode
  #(Example: Exponential distribution, lb = 0 and mode = 0)

  #If the 1st and 2nd points are the same within some numeric tolerance
  if ( abs(initAbs[1] - initAbs[2]) < eps )
    initAbs[1] <- initAbs[2] #Set the 1st point to be the same as the 2nd
  #If the 2nd and 3rd points are the same within some numeric tolerance
  else if ( abs(initAbs[2] - initAbs[3]) < eps )
    initAbs[3] <- initAbs[2] #Set the 3rd point to be the same as the 2nd

  #Drop the repeating points
  initAbs <- unique(initAbs)

  return(initAbs)
}

#This function returns true if the values of the vector vec is decreasing,
#i.e., returns TRUE if the vec[j+1] <= vec[j] for all j=1...length(vec)-1
isDecreasing <- function (vec)
{
  k <- length(vec)

  #The differences vec[j+1] - v[j] <= 0 for all j=1...k-1
  diff <- vec[2:k] - vec[1:k-1]

  #Count the number of times the difference is <= 0
  leq <- sum ( diff < eps )

  #If all the differences are <= 0, then the vector is decreasing
  if ( leq == k-1 )
    return(TRUE)
  else #Otherwise, the vector is not decreasing
    return (FALSE)
}

#This function returns the vector of intersection points of the tangent lines of abscissae.
#Inputs:
#abscissae is the set of points: T_k = { x_1 , x_2, ... , x_k } (assumption: it is sorted)
#hx is the values of h(x) evaluated at the abscissae points
#dhx is the values of h'(x) evaluated at the abscissae points
#Output: vector of intersects of length of k-1
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
    zeros <- which(abs(denominator) < eps)


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
      intersects[zeros] <- (abscissae[zeros+1] + abscissae[zeros])/2
    }
    else #If none of the differences are equal to zero.
    {
      intersects <- numerator/denominator
    }

    return (intersects)
  }
  else #If the h'(x) is not decreasing
  {
    stop("The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold.")
    return (NULL)
  }
}

#Given the new abscissae point, xStar, this function will update the vector
#of the intersects
updateIntersects <- function(abscissae, oldIntersects, hx, dhx, xStar, hxStar, dhxStar)
{
  leq <- (abscissae <= xStar)
  index <- sum(leq)
  
  k <- length(abscissae)
  
  #Drop leftbound and rightbound from the intersects vector
  intersects <- oldIntersects[2:k]
  
  #Initialize the new set of intersects
  newIntersects <- rep(0, k)
  
  #Only if index >= 2, we copy the first section of intersects
  if (index > 1)
    newIntersects[1:(index-1)] <- intersects[1:(index-1)]
  
  #Only if index <= k-2, we copy the last section of intersects
  if (index + 2 <= k )
    newIntersects[(index+2):k] <- intersects[(index+1):(k-1)]
  
  #Intersection of the tangent lines at x_index and xStar
  
  #If the xStar is inserted to the somewhere in the middle of the abscissae
  if (index != 0 && index != k)
  {
    #j = index, we will need to update intersects z_index and z_{index+1}
    xj <- abscissae[index]
    xj1 <- abscissae[index+1]
    
    #Intersection of the tangent lines at x_j and xStar
    if (abs(dhx[index] - dhxStar) > 1e-08) #If the tangent values at x_j and xStar are different
      newIntersects[index] <- (hxStar - hx[index] - xStar*dhxStar + xj*dhx[index])/(dhx[index] - dhxStar)
    else #If the tangent values at x_j and xStar are the same
    {
      newIntersects[index] <- (xStar + xj)/2
    }
    
    #Intersection of the tangent lines at xStar and x_{j+1}
    if (abs(dhxStar - dhx[index+1]) > 1e-08) #If the tangent values at xStar and x_{j+1} are different
      newIntersects[index + 1] <- (hx[index+1] - hxStar - xj1*dhx[index+1] + xStar*dhxStar)/(dhxStar - dhx[index+1])
    else  #If the tangent values at xStar and x_{j+1} are the same
      newIntersects[index + 1] <- (xStar + xj1)/2
  }
  else if (index == 0) #If xStar is being added to the beginning of the abscissae
  {
    xj1 <- abscissae[index+1]
    
    #Intersection of the tangent lines at xStar and x_{j+1}
    if (abs(dhxStar - dhx[index+1]) > 1e-08) #If the tangent values at x_{j+1} and xStar are different
      newIntersects[index + 1] <- (hx[index+1] - hxStar - xj1*dhx[index+1] + xStar*dhxStar)/(dhxStar - dhx[index+1])
    else #If the tangent values at x_{j+1} and xStar are the same
      newIntersects[index + 1] <- (xStar + xj1)/2
  }
  else #index == k, ie xStar is being added to the end of the abscissae
  {
    xj <- abscissae[index]
    
    #Intersect of x_j and xStar
    if (abs(dhx[index] - dhxStar) > 1e-08) #If the tangent values at x_j and xStar are different
      newIntersects[index] <- (hxStar - hx[index] - xStar*dhxStar + xj*dhx[index])/(dhx[index] - dhxStar)
    else #If the tangent values at x_j and xStar are the same
    {
      newIntersects[index] <- (xStar + xj)/2
    }
  }
  
  return(newIntersects)
}

#Given the new abscissae point, xStar, this function will update the vectors
#hx, dhx and abscissae accordingly
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
    if ( (newDhx[index+1] - newDhx[index] > eps) || (newDhx[index+2] - newDhx[index+1] > eps) )
      stop("The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold.")
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
    if ( newDhx[index+2] - newDhx[index+1] > eps )
      stop("The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold.")
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
    if ( newDhx[index+1] - newDhx[index] > eps )
      stop("The sampling function failed the log-concavity test. The derivative vector is not non-increasing within the numeric threshold.")
  }

  return( list(hx = newHx, dhx = newDhx, abscissae = newAbs) )
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
