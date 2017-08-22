# Interface for adaptive sampling routines

#' Adaptive sampling routine.
#'
#' If debug is TRUE, returns a list rather than a vector of samples
#'
#' @param n number of variates to draw
#' @param log_dens log density and its derivative, a function with a single
#'                 parameter
#' @param initialPoints vector of initial points within the domain of log_dens
#' @param minRange lower boundary of support of log_dens
#' @param maxRange upper boundary of support of log_dens
#' @param debug if true, return a list of the algorithm's internal state
#' @param ... additional parameters to pass to the log density function
#' @return a vector of samples, or if debug is TRUE, a list of samples and
#'         state variables
#' @references Gilks, W. R., & Wild, P. (1992). Adaptive rejection sampling
#'             for Gibbs sampling. Applied Statistics, 337–348.
#' @seealso mklogdensf
#' @examples
#' \dontrun{
#'  h <- function(tau) {
#'     y <- log(4*sqrt(2)*sqrt(tau)*exp(-1/2*tau*(3 - 1)^2)/sqrt(pi))
#'     yprime <- 1/2*(-tau*(3 - 1)^2 + 1)/tau
#'     c(y, yprime)
#'   }
#'  xs <- adsample(n = 10, log_dens = h, initialPoints = c(0.5, 1),
#'                 minRange = 0, maxRange = Inf)
#' }
#' @export
adsample <- function(n, log_dens, initialPoints, minRange, maxRange,
                     debug=FALSE, ...) {
  stopifnot(length(initialPoints) >= 2, all(is.numeric(initialPoints)))
  stopifnot(n > 0, class(log_dens) == 'function', minRange < maxRange)
  stopifnot(all(initialPoints > minRange), all(initialPoints < maxRange))
  f <- function(x) log_dens(x, ...)
  if (debug) {
    res <- raw_ad_sample_debug(n, f, initialPoints, minRange, maxRange)
    res$f <- f
    res$n <- n
    res$minRange <- minRange
    res$maxRange <- maxRange
    class(res) <- 'adsample'
    res
  } else {
    raw_ad_sample(n, f, initialPoints, minRange, maxRange)
  }
}

#' Compute a symbolic log density function for use with \code{adsample}.
#'
#' You must specify a *single* variable to differentiate with respect to.
#' \code{mklogensf} returns a function from the abscissa into a bivariate vector
#' of the form \code{( log(density(x)), dlog(density(x))/dx )}.
#'
#' @param expr expression proportional to the density (not log density)
#' @param wrt derivative with respect to this variable
#' @param ... additional parameters to fix in density function
#' @return a function suitable for use with adsample
#' @examples
#'
#' # generated function takes any additional parameters
#' g <- mklogdensf(a + b^2, b)
#' g(3, a = 2)
#'
#' # you can also fix additional parameters when you call mklogdensf
#' g <- mklogdensf(a + b^2, b, a = 2)
#' g(3)
#' g(4)
#'
#' @export
mklogdensf <- function(expr, wrt, ...) {
  # substitute parameters given in ...
  cl <- match.call(expand.dots = FALSE)
  cl$`...` <- cl$wrt <- NULL
  cl$env <- list(...)
  cl[[1]] <- as.name('substitute')
  expr <- eval(cl)

  # change to log density
  ldens <- bquote(log(.(expr)))

  # differentiate with respect to 'wrt'
  wrts <- deparse(substitute(wrt))
  ldderiv <- deriv(ldens, wrts)

  # now manipulate log density & deriv into a function
  # from the abscissa into a bivariate vector of
  # c( log(density(x)), dlog(density(x))/dx )
  function(x, ...) {
    l <- list(...)
    eval(bquote(l[[.(wrts)]] <- .(x)))
    y <- eval(ldderiv, l)
    c(y[1], attr(y, 'gradient'))
  }
}

#' Diagnostic plots for adsample objects
#'
#' Create an adsample object by calling adsample with debug=TRUE.
#' This function plots empirical and actual density side by side.
#'
#' @param x an adsample object
#' @export
#' @rdname adsample
plot.adsample <- function(x, ...) {
  stopifnot(class(x) == 'adsample')

  par(mfrow=c(1,2))

  # LHS plot: density
  xrange <- range(x$samples)
  hist(x$samples, freq=FALSE, xlab = 'x',
       main = paste0('Samples drawn (n=',x$n,')'))
  dens <- Vectorize(function(y) exp(x$f(y)[1]))
  curve(dens, from=xrange[1], to=xrange[2], add=TRUE, col='blue')
  legend('topright', c('Density function'),
         lty=1, col=c('blue'))

  # RHS plot: log density
  log.dens <- Vectorize(function(y) x$f(y)[1])
  curve(log.dens, from=xrange[1], to=xrange[2], col='blue',
        main = paste0('Adaptive approximation (k=',x$k,')'),
        ylab = 'log(density)')
  # evaluation points
  with(x, points(`T`, H))
  # upper hull connects the points (z_i, h(x_i) + (x_i - z_i)*h'(x_i))
  # and remember there is an extra z_0 at the start of object$Z
  z <- x$Z[-1]
  y <- with(x, H + Hprime*(z - `T`))
  lines(z, y, col='darkgreen')
  # lower hull connects (T, H)
  with(x, lines(T, H, col='gray'))
  legend('topright', c('Actual density', 'Upper hull', 'Lower hull'),
         lty=1, col=c('blue', 'darkgreen', 'gray'))
}
