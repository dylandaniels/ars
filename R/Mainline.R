Mainline <- function(n, g, dg=NULL, initialPoints=NULL, leftbound=-Inf, rightbound=Inf) {

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

  print(abscissae)
  abscissae = sort(abscissae)

  # TODO: put in checks for abscissae
  # (1) uniqueness
  # (2) number of abscissae is greater than 2
  # (3) first and last are at h'(x_1) < 0 and h'(x_k) > 0, respectively.

  # TODO Refactor all of the following checks into one function?
  if (length(unique(abscissae)) != length(abscissae)) {
    stop('Elements of abscissae should be unique')
  }

  if (length(abscissae) < 2) {
    stop('You must provide 2 or more abscissae.')
  }

  if (leftbound > abscissae[1] || rightbound < abscissae[length(abscissae)]) {
    stop('Abscissae should be within boundaries.')
  }

  hx <- h(abscissae)
  dhx <- h_der(abscissae)

  if ((is.infinite(leftbound) && dhx[1] <= 0) ||  (is.infinite(rightbound) && dhx[length(dhx)] >= 0)) {
    # make this more descript later.
    stop('Invalid abscissae or integral of function is divergent (cannot be normalized to a valid probability distribution)')
  }

  z <- envelopeIntersectPoints(abscissae, hx, dhx)
  z <- c(leftbound, z, rightbound)

  i <- 0
  samples <- numeric(n)

  while (i < n) {
    #print('-----')
    #print(z)
    #print('-----')
    # TODO: actually make envelope a vector-valued function
    u <- Vectorize(function (x) {
      return(envelope(z, abscissae, x, hx, dhx))
    })
    l <- function (x) {
      return(squeezing(hx, abscissae, x))
    }
    xstar <- sampleFromEnvelope(abscissae, z, u, hx, dhx)
    #print(paste0('xstar=',xstar))
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
