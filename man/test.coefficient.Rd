\name{test.coefficient}
\alias{test.coefficient}
\title{Large-sample Test for a Regression Coefficient in an
Negative Binomial Regression Model}
\usage{
  test.coefficient(nb, dispersion, x, beta0, tests,
    alternative = "two.sided", subset = 1:m,
    print.level = 1)
}
\arguments{
  \item{nb}{an NB data object, output from
  \code{\link{prepare.nb.data}}.}

  \item{dispersion}{dispersion estimates, output from
  \code{\link{estimate.disp}}.}

  \item{x}{an \eqn{n} by \eqn{p} design matrix describing
  the treatment structure}

  \item{beta0}{a \eqn{p}-vector specifying the null
  hypothesis. Non-NA components specify the parameters to
  test and their null values. (Currently, only
  one-dimensional test is implemented, so only one non-NA
  component is allowed).}

  \item{tests}{a character string vector specifying the
  tests to be performed, can be any subset of \code{"HOA"}
  (higher-order asymptotic test), \code{"LR"} (likelihoo
  ratio test), and \code{"Wald"} (Wald test).}

  \item{alternative}{a character string specifying the
  alternative hypothesis, must be one of \code{"two.sided"}
  (default), \code{"greater"} or \code{"less"}.}

  \item{subset}{an index vector specifying on which rows
  should be tests be performed}

  \item{print.level}{a number controlling the amount of
  messages printed: 0 for suppressing all messages, 1
  (default) for basic progress messages, and 2 to 5 for
  increasingly more detailed message.}
}
\value{
  a list containig the following components:
  \item{beta.hat}{an \eqn{m} by \eqn{p} matrix of
  regression coefficient under the full model}
  \item{mu.hat}{an \eqn{m} by \eqn{n} matrix of fitted mean
  frequencies under the full model} \item{beta.tilde}{an
  \eqn{m} by \eqn{p} matrix of regression coefficient under
  the null model} \item{mu.tilde}{an \eqn{m} by \eqn{n}
  matrix of fitted mean frequencies under the null model.}
  \item{HOA, LR, Wald}{each is a list of two
  \eqn{m}-vectors, \code{p.values} and \code{q.values},
  giving p-values and q-values of the corresponding tests
  when that test is included in \code{tests}.}
}
\description{
  \code{test.coefficient} performs large-sample tests
  (higher-order asymptotic test, likelihood ratio test,
  and/or Wald tests) for testing one of the regression
  coefficient in an NB regression model.
}
\details{
  \code{test.coefficient} performs large-sample tests for a
  one-dimensional (\eqn{q=1}) component \eqn{\psi} of the
  \eqn{p}-dimensional regression coefficient \eqn{\beta}.
  The hypothesized value \eqn{\psi_0} of \eqn{\psi} is
  specified by the non-NA component of the vector
  \code{beta0} in the input.

  The likelihood ratio statistic, \deqn{ \lambda = 2
  (l(\hat\beta) - l(\tilde\beta)),} converges in
  distribution to a chi-square distribution with \eqn{1}
  degree of freedom.  The signed square root of the
  likelihood ratio statistic \eqn{\lambda}, also called the
  directed deviance, \deqn{r = sign (\hat\psi - \psi_0)
  \sqrt \lambda} converges to a standard normal
  distribution.

  For testing a one-dimensional parameter of interest,
  Barndorff-Nielsen (1986, 1991) showed that a modified
  directed \deqn{ r^* = r - \frac{1}{r} \log(z)} is, in
  wide generality, asymptotically standard normally
  distributed to a higher order of accuracy than the
  directed deviance \eqn{r} itself, where \eqn{z} is an
  adjustment term. Tests based on high-order asymptotic
  adjustment to the likelihood ratio statistic, such as
  \eqn{r^*} or its approximation, are referred to as
  higher-order asymptotic (HOA) tests. They generally have
  better accuracy than corresponding unadjusted likelihood
  ratio tests, especially in situations where the sample
  size is small and/or when the number of nuisance
  parameters (\eqn{p-q}) is large. The implementation here
  is based on Skovgaard (2001). See Di et al. 2012 for more
  details.
}
\author{
  Yanming Di
}
\references{
  Barndorff-Nielsen, O. (1986): "Infereni on full or
  partial parameters based on the standardized signed log
  likelihood ratio," Biometrika, 73, 307-322

  Barndorff-Nielsen, O. (1991): "Modified signed log
  likelihood ratio," Biometrika, 78, 557-563.

  Skovgaard, I. (2001): "Likelihood asymptotics,"
  Scandinavian Journal of Statistics, 28, 3-32.
}

