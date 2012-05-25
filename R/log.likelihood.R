##' @title The log likelihood of a NB model
##' @param kappa shape parameter
##' @param mu mean parameter
##' @param y a n-vector of NB counts
##' @return the log likelihood of (kappa, mu)
log.likelihood.nb= function(kappa, mu, y) {
  ## p = mu/(mu+kappa);
  ## l0= sum(log(gamma(kappa+y)) - log(gamma(kappa)) - log(gamma(1+y))
  ##       + y * log(mu) + kappa * log(kappa) - (y+kappa) * log(mu+kappa));

  sum(dnbinom(y, kappa, mu=mu, log=TRUE));
}

##' The log likelihood of a NB model, non-integer counts are allowed.
##' @title The log likelihood of a NB model
##' @param kappa shape parameter
##' @param mu mean parameter
##' @param y a n-vector of NB counts
##' @return the log likelihood of (kappa, mu)
log.likelihood.nb.2= function(kappa, mu, y) {
  p = mu/(mu+kappa);
  eps = .Machine$double.eps^0.5;
  p = p + (p<eps)*eps;
  sum(lgamma(kappa+y)- lgamma(kappa) - lgamma(1+y) + y * log(p) + kappa * log(1-p));
}

test.log.likelihood= function() {
  kappa = 2;
  kappa = 1:101;
  mu = 0:100;
  y = 0:100;

  system.time({for (i in 1:100000) {
   l1 = log.likelihood.nb(kappa, mu, y);
  }});
  
  system.time({for (i in 1:100000) {
  l2 = log.likelihood.nb.2(kappa, mu, y);
  }});

  mu=0;
  y = 0;
  y = 1e-10;
  y = 1e-8;
  kappa = 1e-8;
  log.likelihood.nb(kappa, mu, y);
  log.likelihood.nb.2(kappa, mu, y);

}


