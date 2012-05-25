\name{estimate.norm.factors}
\alias{estimate.norm.factors}
\title{Estiamte Normalization Factors}
\usage{
  estimate.norm.factors(counts, lib.sizes, method)
}
\arguments{
  \item{counts}{a matrix of RNA-Seq read counts with rows
  corresponding to gene features and columns corresponding
  to independent biological samples.}

  \item{lib.sizes}{a vector of observed library sizes,
  usually estimated by column totals.}

  \item{method}{a character string specifying the method
  for normalization, can be NULL or "AH2010". If
  method=NULL, the normalization factors will have values
  of 1 (i.e., no normalization is applied); if
  method="AH2010", the normalization method proposed by
  Anders and Huber (2010) will be used.}
}
\value{
  a vector of normalization factors.
}
\description{
  \code{estimate.norm.factors} estiamtes normalization
  factors to account for apparent reduction or increase in
  relative frequencies of non-differentially expressing
  genes as a result of compensating the increased or
  decreased relative frequencies of truly differentially
  expressing genes.
}
\details{
  We take gene expression to be indicated by relative
  frequency of RNA-Seq reads mapped to a gene, relative to
  library sizes (column sums of the count matrix). Since
  the relative frequencies sum to 1 in each library (one
  column of the count matrix), the increased relative
  frequencies of truly over expressed genes in each column
  must be accompanied by decreased relative frequencies of
  other genes, even when those others do not truly
  differently express.

  The concern introduced in Robinson and Oshlack (2010) is
  that this reduction will give a false impression of
  biological relevance. Since the accommodation for
  relative frequencies summing to one is shared equally by
  a very large number of non-differentially expressing
  genes, we suspect that the effect is usually small, but
  examples where it is non-ignorable have been demonstrated
  (Robinson and Oshlack, 2010).

  A simple fix is to compute the relative frequencies
  relative to effective library sizes---library sizes
  multiplied by normalization factors.
}
\author{
  Yanming Di
}
\references{
  {Anders, S. and W. Huber (2010): "Differential expression
  analysis for sequence count data," Genome Biol., 11,
  R106.

  Robinson, M. D. and A. Oshlack (2010): "A scaling
  normalization method for differential expression analysis
  of RNA-seq data," Genome Biol., 11, R25.}
}

