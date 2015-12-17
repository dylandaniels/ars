# Given the new abscissae point, xStar, this function will update the vectors
# hx, dhx and abscissae accordingly
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
