# Evaluates an individual integral between two intersection points
# lb - lower bound
# ub - upper bound
# u - envelope function
# dhx - the value of the h'(x) where x is the abscissae
evaluateIntegral <- function (lb, ub, u, dhx) {
  if (abs(dhx) < 1e-10) {
    return(exp(u(ub))*(ub - lb))
  } else {
    return((exp(u(ub)) - exp(u(lb))) / dhx)
  }
}

# Calculates the values of the integrals used for the sampling function initially
# z - intersection points of envelope function (Note that z[1] may be -Inf, and z[length(z)] may be Inf)
# u - envelope function (u(Inf) and u(-Inf) must be valid calls, both should return 0)
# dhx - vector h'(x) evaluted at each of the abscissae
calculateInitialIntegrals <- function (z, u, dhx) {
  # the evaluated integrals of the envelope function between z[j-1] and z[j] (for each j)
  integrals <- numeric(length(dhx))

  for (i in 1:length(dhx)) {
    integrals[i] <- evaluateIntegral(z[i], z[i+1], u, dhx[i])
  }

  return(integrals)
}

# Samples a new value from the envelope function
# abscissae - abscissae
# z - intersection points
# integrals - the integrated values between each pair of intersection points
# u - envelope function
sampleFromEnvelope <- function (abscissae, z, integrals, u, hx, dhx) {
  unif <- runif(1)
  numPoints <- length(abscissae)
  partSums <- cumsum(integrals)
  normalizedPartSums <-partSums/partSums[numPoints]

  # Find the value t s.t. normalizedPartSums[t - 1] < u < normalizedPartSums[t]
  t <- 1
  while (unif > normalizedPartSums[t]) {
    t <- t + 1
  }

  if (t == 1) {
    partSumsVal <- 0
  } else {
    partSumsVal <- partSums[t-1]
  }

  # Handle case where h'(x_t) = 0
  if (abs(dhx[t]) < 1e-10) {
    sampledValue <- (unif*partSums[numPoints] - partSumsVal) / exp(hx[t]) + z[t]
  } else {
    sampledValue <- ((log(dhx[t] * ((unif*partSums[numPoints]) - partSumsVal) + exp(u(z[t])))
                    - hx[t]) / dhx[t]) + abscissae[t]
  }

  return(sampledValue)
}


# Updates the integrals when a new value is evaluated.
# xStar - newly evaluated point
# oldAbscissae - abscissae from last step
# oldIntegrals - integrals from last step
# newZ - new intersection points
# u - envelope function
# dhx - h'(x) evaluated at new abscissae
updateIntegrals <- function (xStar, oldAbscissae, oldIntegrals, newZ, u, dhx) {
  # Find index i where oldAbscissae[i] <= xStar <= oldAbscissae[i+1]
  i <- sum(oldAbscissae <= xStar)

  integrals <- numeric(length(dhx))

  if (i == 0) { # on left boundary
    integrals[3:length(integrals)] <- oldIntegrals[2:length(oldIntegrals)]
    integrals[1] <- evaluateIntegral(newZ[1], newZ[2], u, dhx[1])
    integrals[2] <- evaluateIntegral(newZ[2], newZ[3], u, dhx[2])

  } else if (i == length(oldAbscissae)) { # on right boundary
    integrals[1:(length(integrals)-2)] <- oldIntegrals[1:(length(oldIntegrals)-1)]
    integrals[length(integrals)-1] <- evaluateIntegral(newZ[length(newZ)-2], newZ[length(newZ)-1], u, dhx[length(newZ)-2])
    integrals[length(integrals)] <- evaluateIntegral(newZ[length(newZ)-1], newZ[length(newZ)], u, dhx[length(newZ)-1])

  } else {  # Interior sample
    if (i > 1) {
      integrals[1:(i-1)] <- oldIntegrals[1:(i-1)]
    }
    if (i < length(oldIntegrals) - 1) {
      integrals[(i+3):length(integrals)] <- oldIntegrals[(i+2):length(oldIntegrals)]
    }
    for (j in 0:2) {
      integrals[i+j] <- evaluateIntegral(newZ[i+j], newZ[i+j+1], u, dhx[i+j])
    }
  }
  return(integrals)
}


# Accepts/rejects a candidate point xStar
# xStar is a point sampled from s_k(x) function
# lowerFun is the l_k(x) function
# upperFun is the u_k(x) function
# logDistFun is h(x) = ln( g(x) ) function
# Returns TRUE if xStar is ACCEPTED
# Returns FALSE if xStar is REJECTED
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
