\alias{nbp.test}
\name{nbp.test}
\title{
  NBP Test for Differential Gene Expression from RNA-Seq Counts
}
\description{
  \code{nbp.test} fits a NBP model to the RNA-Seq counts and performs Robinson and Smyth's
exact NB test on each gene to assess differential gene expression between two groups.
}
\usage{
  nbp.test(counts, grp.ids, grp1, grp2, norm.factors, method.disp = "NBP", print.level, ...)
}
\arguments{
  \item{counts}{
    an \eqn{n} by \eqn{r} matrix of RNA-Seq read counts with rows
  corresponding to genes (exons, gene isoforms, etc) and columns
  corresponding to libraries (independent biological samples).
  }
  \item{grp.ids}{
    an \eqn{r} vector of treatment group identifiers (e.g. integers).
  }
  \item{grp1, grp2}{
    identifiers of the two groups to be compared.
  }
  \item{norm.factors}{
    normalization factors.
  }
  \item{method.disp}{
    method for estimating the dispersion parameter(s). Current options
    are "NBP" (default) or "NB2". 
  }
  \item{print.level}{
    controls the amount of messages printed: 0 for suppressing all
    messages, 1 (default) for basic progress messages, and 2 to 5
    for increasingly more detailed messages.
    }
%%  \item{method.mean}{
%%    MLE or MOM, method for estimating the pooled mean from samples in
%%    two groups.
%%  }
  \item{\dots}{
  optional parameters to be passed to \code{\link{estimate.disp}}, the function that estimates the
  dispersion parameters.
  }
}
\details{

  \code{nbp.test} calls \code{\link{prepare.nbp}} to create the NBP data
  structure, perform optional normalization and adjust library sizes,
  calls \code{\link{estimate.disp}} to estimate the NBP dispersion parameters and
  \code{\link{exact.nb.test}} to perform the exact NB test for differential
  gene expression on each gene. The results are summarized using p-values and q-values
  (FDR).

  \subsection{Overview}{
  For assessing evidence for differential gene expression from RNA-Seq
  read counts, it is critical to adequately model the count variability
  between independent biological replicates.  Negative binomial (NB)
  distribution offers a more realistic model for RNA-Seq count
  variability than Poisson distribution and still permits an exact
  (non-asymptotic) test for comparing two groups.

  For each individual gene, a NB distribution uses a dispersion
  parameter \eqn{\phi_i} to model the extra-Poisson variation between
  biological replicates.  Across all genes, the NBP parameterization of
  the NB distribution (the NBP model) uses two parameters \eqn{(\phi,
  \alpha)} to model extra-Poisson variation over the entire range of
  expression levels. The NBP model allows the NB dispersion parameter to
  be an arbitrary power function of the mean (\eqn{\phi_i =
  \phi\mu_i^{2-\alpha}}). The NBP model includes the Poisson model as a
  limiting case (as \eqn{\phi} tends to \eqn{0}) and the NB2 model as a
  special case (when \eqn{\alpha=2}). Under the NB2 model, the
  dispersion parameter is a constant and does not vary with the mean
  expression levels. NBP model is more flexible and is the recommended
  default option.}
  
  \subsection{Count Normalization}{
  We take gene expression to be indicated by relative frequency of
  RNA-Seq reads mapped to a gene, relative to library sizes (column sums
  of the count matrix). Since the relative frequencies sum to 1 in each
  library (one column of the count matrix), the increased relative
  frequencies of truly over expressed genes in each column must be
  accompanied by decreased relative frequencies of other genes, even
  when those others do not truly differently express. Robinson and
  Oshlack (2010) presented examples where this problem is
  noticeable.

  A simple fix is to compute the relative frequencies relative to
  effective library sizes---library sizes multiplied by normalization
  factors. 
  By default, \code{nbp.test} assumes the normalization factors are 1 (i.e. no
normalization is needed). Users can specify normalization factors through the argument \code{norm.factors}.
  Many authors (Robinson and Oshlack (2010), Anders and Huber
  (2010)) propose to estimate the normalization factors based on the
  assumption that most genes are NOT differentially expressed.
}

 
  \subsection{Library Size Adjustment}{
  The exact test requires that the effective library sizes (column sums
  of the count matrix multiplied by normalization factors) are
  approximately equal. By default, \code{nbp.test} will thin
  (downsample) the counts to make the effective library sizes
  equal. Thinning may lose statistical efficiency, but is unlikely to
  introduce bias.
}
  
}

\value{
  \code{nbp.test} returns a list with the following components:
  \item{counts}{the count matrix, same as input.}
  \item{lib.sizes}{column sums of the count matrix.}
  \item{grp.ids}{a vector of identifiers of treatment groups, same as input.}
  \item{grp1, grp2}{identifiers of the two groups to be compared, same as input.}
  \item{eff.lib.sizes}{effective library sizes, lib.sizes multiplied by the normalization
  factors.}
  \item{pseudo.counts}{count matrix after thinning.}
  \item{pseduo.lib.sizes}{effective library sizes of pseudo counts,
  i.e., column sums of the pseudo count matrix multiplied by the normalization.}
  \item{phi, alpha}{dispersion parameters.}
  \item{pie}{estimated mean relative frequencies of RNA-Seq reads mapped
  to each gene.}
  \item{pooled.pie}{estimated pooled mean of relative frequencies
  in the two groups being compared.}
  \item{expression.levels}{a matrix of estimated gene expression
  levels as indicated by mean relative frequencies of RNA-Seq reads. It has three
  columns \code{grp1}, \code{grp2}, \code{pooled} corresponding to
  the two treatment groups and the pooled mean.}
  \item{log.fc}{base 2 log fold change in mean relative frequency between two groups.}
  \item{p.values}{p-values of the exact NB test applied to each gene (row).}
  \item{q.values}{q-values (estimated FDR) corresponding to the p-values.}
}

\references{
  Di, Y, D. W. Schafer, J. S. Cumbie, and J. H. Chang: "The NBP Negative Binomial Model
  for Assessing Differential Gene Expression from RNA-Seq", SAGMB,
  accepted.
  
  Robinson, M. D. and G. K. Smyth (2007): "Moderated statistical tests
  for assessing differences in tag abundance," Bioinformatics, 23, 2881-2887.

  Robinson, M. D. and G. K. Smyth (2008): "Small-sample estimation of
  negative binomial dispersion, with applications to SAGE data," Biostatistics,
  9, 321-332.

  Anders, S. and W. Huber (2010): "Differential expression analysis for
  sequence count data," Genome Biol., 11, R106.

  Robinson, M. D. and A. Oshlack (2010): "A scaling normalization method
  for differential expression analysis of RNA-seq data," Genome Biol., 11,
  R25.

}
\author{
Yanming Di <diy@stat.oregonstate.edu>, Daniel W Schafer <schafer@stat.oregonstate.edu>
}
\note{
  Due to thinning (random downsampling of counts), two identical calls
  to \code{nbp.test} may yield slightly different results. A random
  number seed can be used to make the results reproducible.  We are
  working on tests based on higher-order asymptotics (HOA) that do not
  require thinning.

  Advanced users can call \code{prepare.nbp}, \code{estimate.disp},
  \code{exact.nb.test} directly to have more control over modeling
  and testing.

%   For example, \code{nbp.test} assumes that different
%  treatment groups have a common set of dispersion parameters
%  \eqn{(\phi,\alpha)}. Advanced users can use \code{estimate.disp} to
%  fit separate NBP model to each treatment group.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, \code{\link{exact.nb.test}}. 
}
\examples{
  ## Load Arabidopsis data
  data(arab);
  
  ## Specify treatment groups and ids of the two groups to be compared
  grp.ids = c(1, 1, 1, 2, 2, 2);

  grp1 = 1;
  grp2 = 2;

  ## Set random number seed to make results reproducible
  set.seed(999);

  ## Fit the NBP model and perform exact NB test for differential gene expression
  res = nbp.test(arab, grp.ids, grp1, grp2, print.level = 5);

  ## Show top ten most differentially expressed genes
  subset = order(res$p.values)[1:10];
  print.nbp(res, subset);

  ## Count the number of differentially expressed genes (e.g. qvalue < 0.05)
  alpha = 0.05;
  sig.res = res$q.values < alpha;
  table(sig.res);

  ## Show boxplots, MA-plot, mean-variance plot and mean-dispersion plot
  plot.nbp(res);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RNA-Seq}
\keyword{negtive binomial}
\keyword{overdispersion}
