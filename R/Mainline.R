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

  #abscissae=initialize() # or get abscissae from user input

  hx <- h(abscissae)
  dhx <- h_der(abscissae)

  i <- 0
  samples <- numeric(n)

  while (i < n) {
    z <- envelopeIntersectPoints(abscissae, hx, dhx)
    u <- function (x) {
      return(envelope(z, abscissae, x, hx, dhx))
    }
    l <- function (x) {
      return(squeezing(hx, abscissae, x))
    }
    xstar <- sampleFromEnvelope(abscissae, z, u, hx, dhx)
    result <- acceptReject(xstar, l, u, h, h_der)
    if (result$step == 2) {
      #updateStep(z, xstar, result, abscissae, hx, dhx)
      newValues <- updateDistVals(abscissae, hx, dhx, xstar, h, h_der)
      hx <- newValues$hx
      dhx <- newValues$dhx
      abscissae <- newValues$abscissae
    }
    if (result$dec) {
      i <- i + 1
      samples[i] <- xstar
    }
  }
}

# Need to update later.
