# z - intersection points of envelope function (Note that z[1] may be -Inf, and z[length(z)] may be Inf)
# u - envelope function (u(Inf) and u(-Inf) must be valid calls, both should return 0)
# dhx - vector h'(x) evaluted at each of the abscissae
partialSums <- function (z, u, dhx) {
  lowerZ <- z[-length(z)] # z[1] through z[length(z)-1]
  upperZ <- z[-1] # z[2] through z[length(z)]

  # the evaluated integrals of the envelope function between z[j] and z[j-1] (for each j)
  integrals <- numeric(length(dhx))

  for (i in 1:length(dhx)) {
    if (abs(dhx[i]) < 1e-10) {
      integrals[i] <- exp(u(z[i+1]))*(z[i+1]-z[i])
    } else {
      integrals[i] <- (exp(u(z[i+1])) - exp(u(z[i]))) / dhx[i]
    }
  }
  return(cumsum(integrals))
}

# Is this function needed??
# This is the function s_k(x) described in the paper
# u - envelope function
normalizedEnvelope <- function (x, u, abscissae, z, derivH) {
  numPoints <- length(abscissae)
  return(exp(u(x))/(partialSums(abscissae, z, u, derivH)[numPoints]))
}

sampleFromEnvelope <- function (abscissae, z, u, hx, dhx) {
  unif <- runif(1)
  numPoints <- length(abscissae)
  partSums <- partialSums(z, u, dhx)
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
      sampledValue <- ((log(dhx[t] * ((unif*partSums[numPoints]) - partSums[t-1]) + exp(u(z[t])))
                      - hx[t]) / dhx[t]) + abscissae[t]
    }
  }
  #print(paste0('sampledValue=', sampledValue))
  return(sampledValue)
}

updatePartialSums <- function (oldPartSums, z, u, dhx) {

}


