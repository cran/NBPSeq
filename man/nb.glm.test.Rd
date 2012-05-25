\name{nb.glm.test}
\alias{nb.glm.test}
\title{Fit Negative Binomial Regression Model and Test for a Regression Coefficient}
\usage{
  nb.glm.test(counts, x, beta0,
    lib.sizes = colSums(counts),
    normalization.method = NULL,
    dispersion.method = "log-linear-rel-mean",
    tests = c("HOA", "LR", "Wald"),
    alternative = "two.sided")
}
\arguments{
  \item{counts}{an m by n matrix of RNA-Seq read counts
  with rows corresponding to gene features and columns
  corresponding to independent biological samples.}

  \item{x}{an n by p design matrix specifiying the
  treatment structure.}

  \item{beta0}{a p-vector specifying the null hypothesis.
  Non-NA components specify the parameters to test and
  their null values. (Currently, only one-dimensional test
  is implemented, so only one non-NA component is
  allowed).}

  \item{lib.sizes}{a p-vector of observed library sizes,
  usually (and by default) estimated by column totals.}

  \item{normalization.method}{a character string specifying
  the method for estimating the normalization factors, can
  be \code{NULL} or \code{"AH2010"}. If \code{method=NULL},
  the normalization factors will have values of 1 (i.e., no
  normalization is applied); if \code{method="AH2010"}, the
  normalization method proposed by Anders and Huber (2010)
  will be used.}

  \item{dispersion.method}{a character string specifying
  the method for estimating the dispersion parameter.
  Currenlty, the only implemented option is
  "log-linear-rel-mean", which assumes that log dispersion
  is a log-linear function of the relative mean.}

  \item{tests}{a character string vector specifying the
  tests to be performed, can be any subset of \code{"HOA"}
  (higher-order asymptotic test), \code{"LR"} (likelihoo
  ratio test), and \code{"Wald"} (Wald test).}

  \item{alternative}{a character string specifying the
  alternative hypothesis, must be one of \code{"two.sided"}
  (default), \code{"greater"} or \code{"less"}.}
}
\value{
  A list containing the following components: \item{data}{a
  list containing the input data matrix with additional
  summary quantities, output from
  \code{\link{prepare.nb.data}}.}
  \item{dispersion}{dispersion estimates and models, output
  from \code{\link{estimate.dispersion}}.} \item{test}{test
  results, output from \code{\link{test.coefficient}}.}
}
\description{
  For each row of the input data matrix, \code{nb.glm.test}
  fits a NB log-linear regression model and perform
  large-sample tests for a one-dimensional regression
  coefficient.
}
\details{
  \code{nbp.glm.test} provides a simple, one-stop interface
  to performing a series of core tasks in regression
  analysis of RNA-Seq data: it calls
  \code{\link{estimate.norm.factors}} to estimate
  normalizaton factors; it calls
  \code{\link{prepare.nb.data}} to create a NB data
  structure; it calls \code{\link{estimate.dispersion}} to
  estimete the NB dispersion; and it calls
  \code{\link{test.coefficient}} to test the regresssion
  coefficient.

  To keep the interface simple, \code{nbp.glm.test}
  provides limited options for fine tuning
  models/parameters in each individual step. For more
  control over individual steps, advanced users can call
  \code{\link{estimate.norm.factors}},
  \code{\link{prepare.nb.data}},
  \code{\link{estimate.dispersion}}, and
  \code{\link{test.coefficient}} directly, or even
  substitue one or more of them with their own versions.
}
\examples{
## Load Arabidopsis data
 data(arab);

 ## Use a subset of arab data for demonstration purpose
 counts = arab[1:100,];

 ## Specify treatment structure
 grp.ids = as.factor(c(1, 1, 1, 2, 2, 2));
 x = model.matrix(~grp.ids);

 ## Specify the null hypothesis
 ## The null hypothesis is beta[1]=0 (beta[1] is the log fold change).
 beta0 = c(NA, 0);

 ## Fit NB regression model and perform large sample tests.
 ## The step can take long if the number of genes is large
 fit = nb.glm.test(counts, x, beta0, normalization.method="AH2010");

 ## The result contains the data, the dispersion estimates and the test results
 names(fit);

 ## Show top ten most differentially expressed genes
 subset = order(fit$test.results$HOA$p.values)[1:10];
 cbind(counts[subset,],
       HOA=fit$test.results$HOA$p.values[subset],
       LR=fit$test.results$LR$p.values[subset],
       Wald=fit$test.results$Wald$p.values[subset]
       );
}
\author{
  Yanming Di
}

