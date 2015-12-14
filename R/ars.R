# lb - lower bound
# ub - upper bound
# u - envelope function
# dhx
evaluateIntegral <- function (lb, ub, u, dhx) {
  if (abs(dhx) < 1e-10) {
    return(exp(u(ub))*(ub - lb))
  } else {
    return((exp(u(ub)) - exp(u(lb))) / dhx)
  }
}

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

# Is this function needed??
# This is the function s_k(x) described in the paper
# u - envelope function
normalizedEnvelope <- function (x, u, abscissae, z, derivH) {
  numPoints <- length(abscissae)
  return(exp(u(x))/(partialSums(abscissae, z, u, derivH)[numPoints]))
}

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


  # Handle case where h'(x_t) = 0
  if (abs(dhx[t]) < 1e-10) {
    if (t == 1) {
      sampledValue <- (unif*partSums[numPoints]) / exp(hx[t]) + z[t]
    } else {
      sampledValue <- (unif*partSums[numPoints] - partSums[t - 1]) / exp(hx[t]) + z[t]
    }
  } else {
    # TODO Dylan: clean this up later
    # partSums <- c(0, partSums)
    if (t == 1) {
      sampledValue <- ((log(dhx[t] * ((unif*partSums[numPoints])) + exp(u(z[t])))
        - hx[t]) / dhx[t]) + abscissae[t]
    } else {
      print(paste0('dhx=',dhx))
      print(paste0('t=',t))
      print(paste0('partSums=', partSums))
      sampledValue <- ((log(dhx[t] * ((unif*partSums[numPoints]) - partSums[t-1]) + exp(u(z[t])))
                      - hx[t]) / dhx[t]) + abscissae[t]
    }
  }

  print(paste0('sampledValue=', sampledValue))
  return(sampledValue)
}


updateIntegrals <- function (xStar, oldAbscissae, oldIntegrals, newZ, u, dhx) {
  # Find index i where oldAbscissae[i] <= xStar <= oldAbscissae[i+1]
  #print(paste0('oldAbscissae=',oldAbscissae))
  print(paste0('newZ=',newZ))
  i <- sum(oldAbscissae <= xStar)
  #print(paste0('length oldIntegrals=', length(oldIntegrals)))
  #print(paste0('length newZ=', length(newZ)))
  print(paste0('dhx=', dhx))


  integrals <- numeric(length(dhx))

  if (i == 0) { # on left boundary
    print('left')
    integrals[3:length(integrals)] <- oldIntegrals[2:length(oldIntegrals)]
    integrals[1] <- evaluateIntegral(newZ[1], newZ[2], u, dhx[1])
    integrals[2] <- evaluateIntegral(newZ[2], newZ[3], u, dhx[2])
  } else if (i == length(oldAbscissae)) { # on right boundary
    print('right')
    integrals[1:(length(integrals)-2)] <- oldIntegrals[1:(length(oldIntegrals)-1)]
    print(integrals)
    print(evaluateIntegral(newZ[length(newZ)-2], newZ[length(newZ)-1], u, dhx[length(newZ)-2]))
    integrals[length(integrals)-1] <- evaluateIntegral(newZ[length(newZ)-2], newZ[length(newZ)-1], u, dhx[length(newZ)-2])
    integrals[length(integrals)] <- evaluateIntegral(newZ[length(newZ)-1], newZ[length(newZ)], u, dhx[length(newZ)-1])
  } else {
    print('interior')
    print(paste0('i=',i))
    # Interior sample
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
  print(paste0('oldIntegrals=', oldIntegrals))
  print(paste0('integrals=', integrals))
  return(integrals)
}


