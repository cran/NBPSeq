##' For each row of the input data matrix, \code{nb.glm.test} fits a
##' NB log-linear regression model and perform large-sample tests for
##' a one-dimensional regression coefficient.
##'
##' \code{nbp.glm.test} provides a simple, one-stop interface to
##' performing a series of core tasks in regression analysis of
##' RNA-Seq data: it calls \code{\link{estimate.norm.factors}} to
##' estimate normalizaton factors; it calls
##' \code{\link{prepare.nb.data}} to create a NB data structure; it
##' calls \code{\link{estimate.dispersion}} to estimete the NB
##' dispersion; and it calls \code{\link{test.coefficient}} to test
##' the regresssion coefficient.
##'
##' To keep the interface simple, \code{nbp.glm.test} provides limited
##' options for fine tuning models/parameters in each individual
##' step. For more control over individual steps, advanced users can
##' call \code{\link{estimate.norm.factors}},
##' \code{\link{prepare.nb.data}}, \code{\link{estimate.dispersion}},
##' and \code{\link{test.coefficient}} directly, or even substitue one
##' or more of them with their own versions.
##'
##' @title Fit Negative Binomial Regression Model and Test for a Regression Coefficient
##' @param counts an m by n matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param x an n by p design matrix specifiying the treatment structure.
##' @param beta0 a p-vector specifying the null hypothesis. Non-NA
##' components specify the parameters to test and their null
##' values. (Currently, only one-dimensional test is implemented, so
##' only one non-NA component is allowed).
##' @param lib.sizes a p-vector of observed library sizes, usually (and by default)
##' estimated by column totals.
##' @param normalization.method a character string specifying the
##' method for estimating the normalization factors, can be
##' \code{NULL} or \code{"AH2010"}. If \code{method=NULL}, the
##' normalization factors will have values of 1 (i.e., no
##' normalization is applied); if \code{method="AH2010"}, the
##' normalization method proposed by Anders and Huber (2010) will be
##' used.
##' @param dispersion.method  a character string specifying  the
##' method for estimating the dispersion parameter. Currenlty, the
##' only implemented option is "log-linear-rel-mean", which assumes
##' that log dispersion is a log-linear function of the relative
##' mean. 
##' @param tests a character string vector specifying the tests to be
##' performed, can be any subset of \code{"HOA"} (higher-order
##' asymptotic test), \code{"LR"} (likelihoo ratio test), and
##' \code{"Wald"} (Wald test).
##' @param alternative a character string specifying the alternative
##' hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or
##' \code{"less"}. 
##' @return  A list containing the following components: \item{data}{a
##' list containing the input data matrix with additional summary
##' quantities, output from \code{\link{prepare.nb.data}}.}
##' \item{dispersion}{dispersion estimates and  models, output from
##' \code{\link{estimate.dispersion}}.}  \item{test}{test results,
##' output from \code{\link{test.coefficient}}.}
##' @author Yanming Di
##' @examples
##'
##'  ## Load Arabidopsis data
##'  data(arab);
##'
##'  ## Use a subset of arab data for demonstration purpose
##'  counts = arab[1:100,];
##'  
##'  ## Specify treatment structure
##'  grp.ids = as.factor(c(1, 1, 1, 2, 2, 2));
##'  x = model.matrix(~grp.ids);
##'
##'  ## Specify the null hypothesis
##'  ## The null hypothesis is beta[1]=0 (beta[1] is the log fold change).
##'  beta0 = c(NA, 0);
##'
##'  ## Fit NB regression model and perform large sample tests.
##'  ## The step can take long if the number of genes is large
##'  fit = nb.glm.test(counts, x, beta0, normalization.method="AH2010");
##'
##'  ## The result contains the data, the dispersion estimates and the test results
##'  names(fit);
##'
##'  ## Show top ten most differentially expressed genes
##'  subset = order(fit$test.results$HOA$p.values)[1:10];
##'  cbind(counts[subset,],
##'        HOA=fit$test.results$HOA$p.values[subset],
##'        LR=fit$test.results$LR$p.values[subset],
##'        Wald=fit$test.results$Wald$p.values[subset]
##'        );
##'
nb.glm.test = function(counts, x, beta0,
            lib.sizes = colSums(counts),
            normalization.method = NULL,
            dispersion.method = "log-linear-rel-mean",
            tests=c("HOA", "LR", "Wald"),
            alternative="two.sided") {

  ## Estimate normalization factors
  norm.factors = estimate.norm.factors(counts, lib.sizes, method=normalization.method);
  
  ## Create an NB object
  nb.data = prepare.nb.data(counts, lib.sizes, norm.factors);

  ## Estimate the dispersion parameters
  dispersion = estimate.dispersion(nb.data, x, method=dispersion.method);

  ## Test for the regression coefficient
  ## debug(test.coefficient);
  test.results = test.coefficient(nb.data, dispersion, x, beta0,
  ##    subset = c(3890, 3891, 3892),
    tests=tests, alternative=alternative);

  ## Return a list summarizing the NB data (nb.data), the dispersion
  ## model and estimates (dispersion) and the test results
  ## (test.results)

  list(data = nb.data, dispersion = dispersion, test.results = test.results);
}

##' \code{estimate.norm.factors} estiamtes normalization factors to
##' account for apparent reduction or increase in relative frequencies
##' of non-differentially expressing genes as a result of compensating
##' the increased or decreased relative frequencies of truly
##' differentially expressing genes.
##'
##' We take gene expression to be indicated by relative frequency of
##' RNA-Seq reads mapped to a gene, relative to library sizes (column
##' sums of the count matrix). Since the relative frequencies sum to 1
##' in each library (one column of the count matrix), the increased
##' relative frequencies of truly over expressed genes in each column
##' must be accompanied by decreased relative frequencies of other
##' genes, even when those others do not truly differently
##' express.
##'
##' The concern introduced in Robinson and Oshlack (2010) is that this
##' reduction will give a false impression of biological relevance.
##' Since the accommodation for relative frequencies summing to one is
##' shared equally by a very large number of non-differentially
##' expressing genes, we suspect that the effect is usually small, but
##' examples where it is non-ignorable have been demonstrated
##' (Robinson and Oshlack, 2010).
##'
##' A simple fix is to compute the relative frequencies relative to
##' effective library sizes---library sizes multiplied by
##' normalization factors.
##'
##' @title Estiamte Normalization Factors
##' @param counts a matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param lib.sizes  a vector of observed library sizes, usually
##' estimated by column totals.
##' @param method a character string specifying the method for
##' normalization, can be NULL or "AH2010". If method=NULL, the
##' normalization factors will have values of 1 (i.e., no
##' normalization is applied); if method="AH2010", the normalization
##' method proposed by Anders and Huber (2010) will be used.
##' @references {Anders, S. and W. Huber (2010): "Differential
##' expression analysis for sequence count data," Genome Biol., 11,
##' R106.
##'
##' Robinson, M. D. and A. Oshlack (2010): "A scaling normalization method for differential expression analysis of RNA-seq data," Genome Biol., 11, R25.}
##' @return a vector of normalization factors.
##' @author Yanming Di
estimate.norm.factors = function(counts, lib.sizes, method) {
  if (is.null(method)) {
    ## No normalization is used
    norm.factors = rep(1, dim(counts)[2]);
  } else {

    if (method=="AH2010") {
      ## The method proposed by Anders and Huber Genome Biology 2010,
      ## 11:R106.

      ## Create a reference column (geometric mean of the read
      ## frequenices in each column)
      m = exp(rowMeans(log(counts)));

      ## Estimate median fold change relative to the reference column
      size.factors = apply(counts, 2, function(y) median((y/m)[m>0]));

      ## The normalization factors are defined relative to the specified library sizes
      norm.factors =  size.factors * mean(lib.sizes) / lib.sizes;

      ## Normalize the normalizaton factors
      norm.factors = norm.factors/exp(mean(log(norm.factors)));
    }
  }

  norm.factors
}

##' Create a NB data structure to hold the RNA-Seq read counts and
##' other relevant information.
##'
##' @title Prepare the NB Data Structure for RNA-Seq Read Counts
##' @param counts an mxn matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param lib.sizes an n-vector of observed library sizes. By default, library sizes
##' are estimated to the column totals of the matrix \code{counts}.
##' @param norm.factors an n-vector of normalization factors. By default, have values 1 (no normalization is applied).
##' @param tags a matrix of tags associated with genes, one row for
##' each gene (having the same number of rows as \code{counts}. By
##' default, row names of the matrix \code{counts} are used as tags.
##' @return  A list containing the following components:
##' \item{counts}{the count matrix, same as input.}
##' \item{lib.sizes}{observed library sizes, same as input.}
##' \item{norm.factors}{normalization factors, same as input.}
##' \item{eff.lib.sizes}{effective library sizes (\code{lib.sizes} x \code{norm.factors}).}
##' \item{rel.frequencies}{relative frequencies (counts divided by the effective library sizes).}
##' \item{tags}{a matrix of gene tags, same as input.}
##' 
##' @author Yanming Di
prepare.nb.data = function(counts,
  lib.sizes=colSums(counts),
  norm.factors=rep(1, dim(counts)[2]),
  tags=matrix(row.names(counts), dim(counts)[1], 1)
  ) {

  eff.lib.sizes = lib.sizes * norm.factors;
  m = dim(counts)[1];
  n = dim(counts)[2];
  rel.freq = counts / (matrix(1, m, 1) %*% matrix(eff.lib.sizes, 1, n));

  nb.data = list(
    counts = counts,
    lib.sizes = lib.sizes,
    norm.factors = norm.factors,
    eff.lib.sizes = eff.lib.sizes,
    rel.frequencies=rel.freq,
    tags = tags);
}

##' Estimate NB dispersion by modeling it as a function of the mean
##' frequency and library sizes.
##'
##' We use a negative binomial (NB) distribution to model the read
##' frequency of gene \eqn{i} in sample \eqn{j}.  A negative binomial
##' (NB) distribution uses a dispersion parameter \eqn{\phi_{ij}} to
##' model the extra-Poisson variation between biological replicates.
##' Under the NB model, the mean-variance relationship of a single
##' read count satisfies \eqn{\sigma_{ij}^2 = \mu_{ij} + \phi_{ij}
##' \mu_{ij}^2}.  Due to the typically small sample sizes of RNA-Seq
##' experiments, estimating the NB dispersion \eqn{\phi_{ij}} for each
##' gene \eqn{i} separately is not reliable.  One can pool information
##' across genes and biological samples by modeling \eqn{\phi_{ij}} as
##' a function of the mean frequencies and library sizes. The
##' "log-linear-rel-mean" method assumes a parametric dispersion model
##' \deqn{\phi_{ij} = \alpha_0 + \alpha_1 \log(\pi_{ij}),} where
##' \eqn{\pi_{ij} = \mu_{ij}/(N_j R_j)} is the relative mean frequency
##' after normalization. The parameters \eqn{(\alpha_0, \alpha_1)} in
##' this dispersion model are estimated by maximizing the adjusted
##' profile likelihood.
##'
##' @title Estimate Negative Binomial Dispersion
##' @param nb.data output from \code{prepare.nb.data}.
##' @param x a design matrix specifiying the mean structure of each row.
##' @param method the method for estimating the dispersion
##' parameter. Currenlty, the only implemented option is
##' "log-linear-rel-mean", which assumes that log dispersion is a
##' log-linear function of the relative mean. 
##' @param ... additional parameters.
##' @return a list of two components:
##' \item{estiamtes}{dispersion
##' estimates for each read count, a matrix of the same dimensions as
##' the \code{counts} matrix in nb.data.}
##' \item{models}{a list of
##' dispersion models, NOT intended for use by end users.}
##' @author Yanming Di
estimate.dispersion = function(nb.data, x, method="log-linear-rel-mean", ...) {
 
 res = estimate.disp.mapl.nbp(nb.data$counts, nb.data$eff.lib.sizes, x,
                         ...);

  dispersion = list(estimates = res$phi, models = list(res));
}

##' \code{test.coefficient} performs large-sample tests (higher-order
##' asymptotic test, likelihood ratio test, and/or Wald tests) for
##' testing one of the regression coefficient in an NB regression
##' model. 
##'
##' \code{test.coefficient} performs large-sample tests for a
##' one-dimensional (\eqn{q=1}) component \eqn{\psi} of the
##' \eqn{p}-dimensional regression coefficient \eqn{\beta}. The
##' hypothesized value \eqn{\psi_0} of \eqn{\psi} is specified by the
##' non-NA component of the vector \code{beta0} in the input.
##'
##' The likelihood ratio statistic, \deqn{ \lambda = 2 (l(\hat\beta) -
##' l(\tilde\beta)),} converges in distribution to a chi-square
##' distribution with \eqn{1} degree of freedom.  The signed square
##' root of the likelihood ratio statistic \eqn{\lambda}, also called
##' the directed deviance, \deqn{r = sign (\hat\psi - \psi_0) \sqrt
##' \lambda} converges to a standard normal distribution.
##'
##' For testing a one-dimensional parameter of interest,
##' Barndorff-Nielsen (1986, 1991) showed that a  modified directed
##' \deqn{ r^* = r - \frac{1}{r} \log(z)} is, in wide generality,
##' asymptotically standard normally distributed to a higher order of
##' accuracy than the directed deviance \eqn{r} itself, where \eqn{z}
##' is an adjustment term. Tests based on high-order asymptotic
##' adjustment to the likelihood ratio statistic, such as \eqn{r^*} or
##' its approximation, are referred to as higher-order asymptotic
##' (HOA) tests. They generally have better accuracy than
##' corresponding unadjusted likelihood ratio tests, especially in
##' situations where the sample size is small and/or when the number
##' of nuisance parameters (\eqn{p-q}) is large. The implementation
##' here is based on Skovgaard (2001). See Di et al. 2012 for more
##' details.
##'
##' @title Large-sample Test for a Regression Coefficient in an
##' Negative Binomial Regression Model
##' @param nb an NB data object, output from \code{\link{prepare.nb.data}}.
##' @param dispersion dispersion estimates, output from \code{\link{estimate.disp}}.
##' @param x an \eqn{n} by \eqn{p} design matrix describing the treatment structure
##' @param beta0 a \eqn{p}-vector specifying the null hypothesis. Non-NA
##' components specify the parameters to test and their null
##' values. (Currently, only one-dimensional test is implemented, so
##' only one non-NA component is allowed).
##' @param tests a character string vector specifying the tests to be
##' performed, can be any subset of \code{"HOA"} (higher-order
##' asymptotic test), \code{"LR"} (likelihoo ratio test), and
##' \code{"Wald"} (Wald test).
##' @param alternative a character string specifying the alternative
##' hypothesis, must be one of \code{"two.sided"} (default),
##' \code{"greater"} or \code{"less"}. 
##' @param subset an index vector specifying on which rows should be
##' tests be performed
##' @param print.level a number controlling the amount of messages
##' printed: 0 for suppressing all messages, 1 (default) for basic
##' progress messages, and 2 to 5 for increasingly more detailed
##' message.
##' @references Barndorff-Nielsen, O. (1986): "Infereni on full or
##' partial parameters based on the standardized signed log likelihood
##' ratio," Biometrika, 73, 307-322
##'
##' Barndorff-Nielsen, O. (1991): "Modified signed log likelihood
##' ratio," Biometrika, 78, 557-563.
##' 
##' Skovgaard, I. (2001): "Likelihood asymptotics," Scandinavian
##' Journal of Statistics, 28, 3-32.
##'
##' @return a list containig the following components:
##' \item{beta.hat}{an \eqn{m} by \eqn{p} matrix of regression
##' coefficient under the full model} \item{mu.hat}{an \eqn{m} by
##' \eqn{n} matrix of fitted mean frequencies under the full model}
##' \item{beta.tilde}{an \eqn{m} by \eqn{p} matrix of regression
##' coefficient under the null model} \item{mu.tilde}{an \eqn{m} by
##' \eqn{n} matrix of fitted mean frequencies under the null model.}
##' \item{HOA, LR, Wald}{each is a list of two \eqn{m}-vectors,
##' \code{p.values} and \code{q.values}, giving p-values and q-values
##' of the corresponding tests when that test is included in
##' \code{tests}.}
##'
##' @author Yanming Di
test.coefficient = function(nb, dispersion, x, beta0,
  tests,
  alternative="two.sided",
  subset = 1:m,
  print.level=1) {

  if (print.level>0)
    message("HOA test for regression coefficients.");

  ## Determine the alternative
  alternatives = c("two.sided", "less", "greater");
  alt = pmatch(alternative, c("two.sided", "less", "greater"));
  if (is.na(alt)) {
    warning('alternative should be "two.sided", "less" or "greater"');
    return(NULL);
  }

  counts = nb$counts;
  lib.sizes = nb$eff.lib.sizes;
  phi = dispersion$estimates;
  m = dim(counts)[1];
  n = dim(counts)[2];
  p = length(beta0);
    
  obj = list(
    beta.hat = matrix(NA, m, p),
    mu.hat = matrix(NA, m, n),
    beta.tilde = matrix(NA, m, p),
    mu.tilde = matrix(NA, m, n));

  if ("HOA" %in% tests) {
    HOA = list(p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if ("LR" %in% tests) {
    LR = list(p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if ("Wald" %in% tests) {
    Wald = list(p.values=rep(NA, m), q.values=rep(NA, m));
  }

  if (print.level>1) pb = txtProgressBar(style=3);
  for (i in subset) {
    if (print.level>1) setTxtProgressBar(pb, i/m);

    res.hoa = try(
      hoa.1d(counts[i,], lib.sizes, x, phi[i,], beta0,
             alternative=alternative,
             print.level=print.level-1),
      silent=TRUE); 

    if (!("try-error" %in% class(res.hoa))) {
      ## print(i);
      obj$beta.hat[i,] = res.hoa$beta.hat;
      obj$mu.hat[i,] = res.hoa$mu.hat;
      obj$beta.tilde[i,] = res.hoa$beta.tilde;
      obj$mu.tilde[i,] = res.hoa$mu.tilde;

      if ("HOA" %in% tests) {
        HOA$p.values[i] = res.hoa$pstar;
      }
      
      if ("LR" %in% tests) {
        LR$p.values[i] = res.hoa$p;
      }

      if ("Wald" %in% tests) {
        Wald$p.values[i] = res.hoa$p.wald;
      }
    }
  }

  if (print.level>1) close(pb);

  ## Compute q-values
  compute.q.values = function(p.values) {
    q.values = rep(NA, length(p.values));
    id = !is.na(p.values);
    q.values[id] = qvalue(p.values[id])$qvalues;
  }

  if ("HOA" %in% tests) {
    HOA$q.values = compute.q.values(HOA$p.values);
    obj$HOA = HOA;
  }
  
  if ("LR" %in% tests) {
    LR$q.values = compute.q.values(LR$p.values);
    obj$LR = LR;
  }

  if ("Wald" %in% tests) {
    Wald$q.values = compute.q.values(Wald$p.values);
    obj$Wald = Wald;
  }

  obj;
}

