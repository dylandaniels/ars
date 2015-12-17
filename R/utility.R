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
  startPt <- 1e-04 + (lb + ub)/2

  #Assumption: Left Bound and Right Bound are correctly given!!
  if ( is.infinite(leftbound) && is.infinite(rightbound) )
  {
    #Minimize -h(x) using nlm function, where the max number of iterations is 50
    tryCatch({ sol <- nlm(negh, startPt)
    modeVal <- sol$estimate },
    error = function(e) { stop("Unable to guess initial points. Please provide them manually.") },
    warning = function(w) { stop(print(w))})
  }
  else #If there are finite bounds at either left or right sides
  {
    tryCatch({ sol <- optimx(startPt, negh, method = "L-BFGS-B", lower = leftbound, upper = rightbound)
    modeVal <- sol[1,1] },
    error = function(e) { stop("Unable to guess initial points. Please provide them manually.") },
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
      initAbs[1] <- (0.2*leftbound + 0.8*modeVal)

    #If right bound is infinity
    if (is.infinite(rightbound))
      initAbs[3] <- modeVal + 1
    else
      initAbs[3] <- (0.2*rightbound + 0.8*modeVal)
  }
  else #Both bounds are finite
  {
    initAbs[1] <- (0.2*leftbound + 0.8*modeVal)
    initAbs[3] <- (0.2*rightbound + 0.8*modeVal)
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
