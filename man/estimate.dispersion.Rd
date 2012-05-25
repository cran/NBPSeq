\name{estimate.dispersion}
\alias{estimate.dispersion}
\title{Estimate Negative Binomial Dispersion}
\usage{
  estimate.dispersion(nb.data, x,
    method = "log-linear-rel-mean", ...)
}
\arguments{
  \item{nb.data}{output from \code{prepare.nb.data}.}

  \item{x}{a design matrix specifiying the mean structure
  of each row.}

  \item{method}{the method for estimating the dispersion
  parameter. Currenlty, the only implemented option is
  "log-linear-rel-mean", which assumes that log dispersion
  is a log-linear function of the relative mean.}

  \item{...}{additional parameters.}
}
\value{
  a list of two components: \item{estiamtes}{dispersion
  estimates for each read count, a matrix of the same
  dimensions as the \code{counts} matrix in nb.data.}
  \item{models}{a list of dispersion models, NOT intended
  for use by end users.}
}
\description{
  Estimate NB dispersion by modeling it as a function of
  the mean frequency and library sizes.
}
\details{
  We use a negative binomial (NB) distribution to model the
  read frequency of gene \eqn{i} in sample \eqn{j}.  A
  negative binomial (NB) distribution uses a dispersion
  parameter \eqn{\phi_{ij}} to model the extra-Poisson
  variation between biological replicates. Under the NB
  model, the mean-variance relationship of a single read
  count satisfies \eqn{\sigma_{ij}^2 = \mu_{ij} + \phi_{ij}
  \mu_{ij}^2}.  Due to the typically small sample sizes of
  RNA-Seq experiments, estimating the NB dispersion
  \eqn{\phi_{ij}} for each gene \eqn{i} separately is not
  reliable.  One can pool information across genes and
  biological samples by modeling \eqn{\phi_{ij}} as a
  function of the mean frequencies and library sizes. The
  "log-linear-rel-mean" method assumes a parametric
  dispersion model \deqn{\phi_{ij} = \alpha_0 + \alpha_1
  \log(\pi_{ij}),} where \eqn{\pi_{ij} = \mu_{ij}/(N_j
  R_j)} is the relative mean frequency after normalization.
  The parameters \eqn{(\alpha_0, \alpha_1)} in this
  dispersion model are estimated by maximizing the adjusted
  profile likelihood.
}
\author{
  Yanming Di
}

