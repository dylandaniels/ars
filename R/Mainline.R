Mainline <- function(n, g, dg=NULL, initialPoints=NULL, leftbound=-Inf, rightbound=Inf, showPlot=FALSE) {

  abscissae <- initialPoints

  h <- function (y) {
    return(log(g(y)))
  }

  if (is.null(dg)) {
    h_der <- function (y) {
      # Look into replacing this with grad() function?
      derivatives <- sapply(y, function (x) {
        env <- new.env()
        assign('x', x, envir = env)
        result = tryCatch({
          diag(attributes(numericDeriv(quote(h(x)), 'x', env))$gradient)
        }, error = function(e) {
          stop('derivates of initial abscissae could not be evaluted numerically.')
        })
        return(result)
      })
      return(derivatives)
    }
  } else {
    h_der <- convertDerivToLog(g, dg)
  }

  #If the user does not provide the initial points, then run the
  #findInitPoints function
  if (is.null(initialPoints))
    abscissae <- findInitPoints(h, leftbound, rightbound)
    message(paste(c('No initial points given by user. Guessing initial points:', abscissae), collapse=" "))



  abscissae = sort(abscissae)
  hx <- h(abscissae)
  dhx <- h_der(abscissae)
  precheck(abscissae, dhx, leftbound, rightbound)

  # TODO Refactor all of the following checks into one function?

  z <- envelopeIntersectPoints(abscissae, hx, dhx)
  z <- c(leftbound, z, rightbound)

  i <- 0
  samples <- numeric(n)

  integrals <- NULL

  u <- Vectorize(function (x) {
    return(envelope(z, abscissae, x, hx, dhx))
  })

  l <- function (x) {
    return(squeezing(hx, abscissae, x))
  }

  while (i < n) {

    if (is.null(integrals)) {
      integrals <- calculateInitialIntegrals(z, u, dhx)
    }

    xstar <- sampleFromEnvelope(abscissae, z, integrals, u, hx, dhx)
    result <- acceptReject(xstar, l, u, h, h_der)
    if (result$step == 2) {
      z <- updateIntersects(abscissae, z, hx, dhx, xstar, result$hx, result$dhx)
      z <- c(leftbound, z, rightbound)

      newValues <- updateDistVals(abscissae, hx, dhx, xstar, result$hx, result$dhx)
      hx <- newValues$hx
      dhx <- newValues$dhx
      oldAbscissae <- abscissae
      abscissae <- newValues$abscissae

      u <- function (x) {
        return(envelope(z, abscissae, x, hx, dhx))
      }

      l <- function (x) {
        return(squeezing(hx, abscissae, x))
      }

      integrals <- updateIntegrals(xstar, oldAbscissae, integrals, z, u, dhx)
    }
    if (result$dec) {
      i <- i + 1
      samples[i] <- xstar
    }
  }

  if (showPlot) {
    x <- seq(0,5,0.01)
    u <- Vectorize(u)
    l <- Vectorize(l)
    plot(x, exp(u(x)), type='l')
    points(x, exp(l(x)), col='blue',type='l')
  }


  return (samples)
}

# Need to update later.
