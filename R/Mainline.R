Mainline <- function(n, g, abscissae=NULL, leftbound=-Inf, rightbound=Inf) {
  h <- function (y) {
    return(log(g(y)))
  }

  h_der <- function (y) {
    # Look into replacing this with grad() function?
    derivatives <- sapply(y, function (x) {
      env <- new.env()
      assign('x', x, envir = env)
      return(diag(attributes(numericDeriv(quote(h(x)), 'x', env))$gradient))
    })
    return(derivatives)
  }

  # TODO: put in checks for abscissae
  # (1) uniqueness
  # (2) number of abscissae is greater than 2
  # (3) first and last are at h'(x_1) < 0 and h'(x_k) > 0, respectively.

  #abscissae=initialize() # or get abscissae from user input

  hx <- h(abscissae)
  dhx <- h_der(abscissae)
  
  z <- envelopeIntersectPoints(abscissae, hx, dhx)
  z <- c(leftbound, z, rightbound)

  i <- 0
  samples <- numeric(n)

  while (i < n) {
    print('-----')
    print(z)
    print('-----')
    # TODO: actually make envelope a vector-valued function
    u <- Vectorize(function (x) {
      return(envelope(z, abscissae, x, hx, dhx))
    })
    l <- function (x) {
      return(squeezing(hx, abscissae, x))
    }
    xstar <- sampleFromEnvelope(abscissae, z, u, hx, dhx)
    print(paste0('xstar=',xstar))
    result <- acceptReject(xstar, l, u, h, h_der)
    if (result$step == 2) {
      #updateStep(z, xstar, result, abscissae, hx, dhx)
      z <- updateIntersects(abscissae, z, hx, dhx, xstar, result$hx, result$dhx)
      z <- c(leftbound, z, rightbound)
      newValues <- updateDistVals(abscissae, hx, dhx, xstar, result$hx, result$dhx)
      hx <- newValues$hx
      dhx <- newValues$dhx
      abscissae <- newValues$abscissae
    }
    if (result$dec) {
      i <- i + 1
      samples[i] <- xstar
    }
  }
  return (samples)
}

# Need to update later.
