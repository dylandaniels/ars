#' Adaptive Rejection Sampling
#'
#' Generate samples from a distribution via the adapative rejection sampling algorithm. Adaptive rejection
#' sampling dynamically builds a ``rejection region,'' so that the user does not need to explicitly supply one.
#'
#' @param n Number of samples.
#' @param g (Unnormalized) density function to sample from. \code{g} must be (non-strictly) log-concave.
#' @param dg Derivative of the density function. If not supplied, numeric differentiation is attempted.
#' @param initialPoints A vector of the initial abscissae to generate the envelope and squeezing functions.
#'  If not supplied, an optimization routine will attempt to find initial points.
#' @param leftbound The lower bound of the domain of \code{g}.
#' @param rightbound The upper bound of the domain of \code{g}.
#'
#' @details
#' The function \code{g} must be an log-concave (unnormalized) density function. Accurate \code{leftbound}
#' \code{rightbound} values must also be supplied.
#'
#' The \code{initialPoints} are used to construct
#' the envelope and squeezing functions described in the paper (Gilks & Wild, 1992) below. If
#' not supplied, an algorithm is run which attempts to find the mode of \code{g} and generates
#' initial points near the mode. If the user wishes to supply initial points (recommended),
#' at least one initial point must be given where the derivative of \code{g} is less than zero, and
#' another where the derivative of \code{g} is greater than zero, unless the function is monotonic.
#'
#' @return Returns a vector of \code{n} samples from \code{g}.
#'
#' @references
#' Gilks, W. R., & Wild, P. (1992). Adaptive rejection sampling for Gibbs sampling.
#' \emph{Applied Statistics}, 337-348.
#'
#' @examples
#' # Sample 10 points from Normal(0,1)
#' ars(10, dnorm, initialPoints=c(-1,1))
#'
#' # Sample 15 points from Uniform[0,1]
#' ars(15, dunif, initialPoints=c(0.2, 0.3, 0.8), leftbound=0, rightbound=1)
#'
#' # Define a quadratic distribution
#' f <- function (x) {
#'  return(x^2)
#' }
#'
#' df <- function (x) {
#'  return(2 * x)
#' }
#'
#' # Sample 5 points from a quadratic distribution
#' # Note that ars() takes care of normalization
#' ars(5, f, df, initialPoints = c(1,2), leftbound=0, rightbound=10)
#' @export
ars <- function(n, g, dg=NULL, initialPoints=NULL, leftbound=-Inf, rightbound=Inf) {

  abscissae <- initialPoints

  h <- function (y) {
    return(log(g(y)))
  }

  if (is.null(dg)) {
    h_der <- function (y) {
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
  if (is.null(initialPoints)) {
    abscissae <- findInitPoints(h, leftbound, rightbound)
    message(paste(c('No initial points given by user. Guessing initial points:', abscissae), collapse=" "))
  }


  abscissae = sort(abscissae)
  hx <- h(abscissae)
  dhx <- h_der(abscissae)
  precheck(abscissae, dhx, leftbound, rightbound)

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

  return (samples)
}
