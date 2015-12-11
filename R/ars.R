# z - intersection points of envelope function (Note that z[1] may be -Inf, and z[length(z)] may be Inf)
# u - envelope function (u(Inf) and u(-Inf) must be valid calls, both should return 0)
# dhx - vector h'(x) evaluted at each of the abscissae
partialSums <- function (z, u, dhx) {
  lowerZ <- z[-length(z)] # z[1] through z[length(z)-1]
  upperZ <- z[-1] # z[2] through z[length(z)]
  # the evaluated integrals of the envelope function between z[j] and z[j-1] (for each j)
  integrals <- (exp(u(upperZ)) - exp(u(lowerZ))) / dhx
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
  partSums <- partialSums(abscissae, z, u, dhx)
  normalizedPartSums <-partSums/partSums[numPoints]

  # Find the value t s.t. normalizedPartSums[t - 1] < u < normalizedPartSums[t]
  t <- 1
  while (unif > normalizedPartialSums[t]) {
    t <- t + 1
  }

  sampledValue <- ((log(dhx[t] * ((unif/partSums[numPoints]) - partSums[t-1]) + exp(u(z[t])))
                    - hx[t]) / dhx[t]) + abscissae[t]
  return(sampledValue)
}


