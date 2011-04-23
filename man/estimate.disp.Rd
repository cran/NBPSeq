\name{estimate.disp}
\alias{estimate.disp}
\title{
Estimate Negative Binomial Dispersion Parameters
}
\description{
Fit a NBP (or NB2) model to the RNA-Seq count data. Estimate negative binomial
dispersion parameter(s).
}
\usage{
estimate.disp(obj, method = "NBP", print.level = 1, ...)
}
\arguments{
    \item{obj}{output from \code{\link{prepare.nbp}}. 
    }
    \item{method}{
	"NBP" (default) or "NB2", the model for the count variance.
    }
    \item{print.level}{
	controls the amount of messages printed: 0 for
	    suppressing all messages, 1 for basic progress messages, larger
	    values for more detailed messages.
    }
    \item{\dots}{
	additional parameters controlling the estimation of the parameters.
    }
}

\value{
    The list \code{obj} from the input with the following added components:
	\item{phi, alpha}{dispersion parameters.}
    \item{pie}{a matrix the same size as \code{counts}, estimated mean relative frequencies.}
}

\details{
  For each individual gene \eqn{i}, a negative binomial (NB) distribution uses a
  dispersion parameter \eqn{\phi_i} to model the extra-Poisson variation
  between biological replicates: the NB model imposes a mean-variance
  relationship \eqn{\sigma_i^2 = \mu_i + \phi_i \mu_i^2}.
  Across all genes, the NBP
  parameterization of the NB distribution (the NBP model) uses two
  parameters \eqn{(\phi, \alpha)} to model extra-Poisson variation over
  the entire range of expression levels. The NBP model allows the NB
  dispersion parameter to be an arbitrary power function of the mean
  (\eqn{\phi_i = \phi\mu_i^{2-\alpha}}). The NBP model includes the
  Poisson model as a limiting case (as \eqn{\phi} tends to \eqn{0}) and
  the NB2 model as a special case (when \eqn{\alpha=2}). Under the NB2
  model, the dispersion parameter is a constant and does not vary with
  the mean expression levels. NBP model is more flexible and is the
  recommended default option.
  
  The dispersion parameters are estimated from the pseudo counts (counts
  adjusted to have same effective library sizes).  The parameters are
  estimated by maximizing the log conditional likelihood of \eqn{(\phi,
  \alpha)} given the row sums. The log conditional likelihood is
  computed for each gene in each treatment group and then summed over
  genes and treatment groups.
}

\note{
  Users should call \code{\link{prepare.nbp}} before calling this
  function. The function \code{\link{prepare.nbp}} will normalize
  the counts and adjust the counts so that the effective library sizes are
  approximately the same (computing the conditional likelihood requires
  the library sizes to be the same).
}

\references{
  Di, Y, D. W. Schafer, J. S. Cumbie, and J. H. Chang: "The NBP Negative Binomial Model
  for Assessing Differential Gene Expression from RNA-Seq", SAGMB,
  accepted.
}
\author{
Yanming Di <diy@stat.oregonstate.edu>, Daniel W Schafer <schafer@stat.oregonstate.edu>
}

\seealso{
  See Also \code{\link{nbp.test}}.
}
\examples{

  ## Load Arabidopsis data
  data(arab);

  ## Specify treatment groups 
  grp.ids = c(1, 1, 1, 2, 2, 2);

  ## Prepare an NBP object, adjust the library sizes by thinning the counts.
  set.seed(999);
  obj = prepare.nbp(arab, grp.ids, print.level=5);

  ## Estimate the NBP dispersion parameters
  obj = estimate.disp(obj, print.level=5);
  
  ## Print the NBP object
  print.nbp(obj);

}
