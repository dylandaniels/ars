# Thought:
# It probably makes sense to store the following variables in some sort of class/wrapper so that
# multiple functions can access them. All of the following variables remain invariant between updating
# steps:
# - abscissae
# - z
# - derivH
# - h
# - numPoints (or k)
# - partialSums

# abscissae - the x_1 through x_k values already evaluated
# z - intersection points of envelope function (Note that z[1] may be -Inf, and z[length(z)] may be Inf)
# u - envelope function (u(Inf) and u(-Inf) must be valid calls, both should return 0)
# derivH - (log f)'
partialSums <- function (abscissae, z, u, derivH) {
  lowerZ <- z[-length(z)] # z[1] through z[length(z)-1]
  upperZ <- z[-1] # z[2] through z[length(z)]
  # the evaluated integrals of the envelope function between z[j] and z[j-1] (for each j)
  integrals <- (u(upperZ) - u(lowerZ)) / derivH(abscissae)
  return(cumsum(integrals))
}


# This is the function s_k(x) described in the paper
# u - envelope function
normalizedEnvelope <- function (x, u, abscissae, z, derivH) {
  numPoints <- length(abscissae)
  return(u(x)/(partialSums(abscissae, z, u, derivH)[numPoints]))
}

sampleFromEnvelope <- function (abscissae, z, u, h, derivH) {
  unif <- runif(1)
  numPoints <- length(abscissae)
  partSums <- partialSums(abscissae, z, u, derivH)
  normalizedPartSums <-partSums/partSums[numPoints]

  # Find the value t s.t. normalizedPartSums[t - 1] < u < normalizedPartSums[t]
  t <- 1
  while (unif > normalizedPartialSums[t]) {
    t <- t + 1
  }

  sampledValue <- ((log(derivH(abscissae[t]) * ((unif/partSums[numPoints]) - partSums[t-1]) + u(z[t]))
                    - h(abscissae[t])) / derivH(abscissae[t])) + abscissae[t]
  return(sampledValue)
}


