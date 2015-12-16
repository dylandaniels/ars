eps <- 1e-04

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

