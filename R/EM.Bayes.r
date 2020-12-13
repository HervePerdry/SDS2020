#' Estimation of the parameters of a mixture of two univariate normal distributions
#'
#' @param x a vector of data
#' @param Beta.pars parameters of the prior for pi
#' @param lambda parameter for the prior for sigma2^2/sigma1^2
#' @param eps criterion of convergence
#' @param verbose verbosity of the function
#'
#' @details This functions gives a Maximum A Posteriori for the parameter of a 
#' Gaussian mixture \eqn{pi*N(mu1, sigma1^2) + (1-pi)*N(mu2, sigma2^2)}.
#' @details The prior on pi is a Beta distribution. There's a prior on 
#' \eqn{r = sigma1^2/sigma2^2}, which is a \eqn{gamma(lambda + 1, lambda)}.
#' 
#' @details The default values for the prior parameters correspond to flat priors.
#' 
#' 
#' @return a vector with named components
#' @export
#'
#' @examples
#' set.seed(1);
#' x <- c( rnorm(100, 0), rnorm(100, 1) )
#' EM.Bayes(x)
#' # Recommanded non flat priors:
#' EM.Bayes(x, c(2,2), 1)

EM.Bayes <- function(x, Beta.pars = c(1,1), lambda = 0, eps = 1e-8, verbose = TRUE) {
  .Call(`_SDS2020_EMbayes`, x, Beta.pars[1]-1, Beta.pars[2]-1, lambda, eps, verbose)
}
