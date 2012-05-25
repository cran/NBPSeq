## Thin the counts
thin.counts = function(y, current.lib.sizes = colSums(y),
  target.lib.sizes = min(current.lib.sizes)) {
  ## Arguments:
  ##
  ##   y: a n x r matrix of counts
  ##
  ##   current.lib.sizes: an r vector of currrent library sizes
  ##
  ##   norm.factors: an r vector of normalization factors
  ##
  ##   target.lib.sizes: target library sizes 

  n = dim(y)[1];
  ncols = dim(y)[2];

  if (any (target.lib.sizes > current.lib.sizes)) {
    stop("Error: taget.lib.sizes > current.lib.sizes!");
  }

  if (length(target.lib.sizes) == 1) {
    target.lib.sizes = rep(target.lib.sizes, ncols);
  }

  y.thinned = y;

  for (k in 1:ncols) {
    keep.p =  target.lib.sizes[k]/current.lib.sizes[k]; # proportion of reads to keep
    y.thinned[,k] = rbinom(n, y[,k], keep.p);
  }

  ## target.lib.sizes = rep(mean(colSums(y.thinned)), ncols);

  list(counts =y.thinned, lib.sizes = target.lib.sizes);
}


## Create an NBP object, perform coutn normalization, and adjust for
## unequal library sizes.
prepare.nbp = function(counts, grp.ids, lib.sizes = colSums(counts),
  norm.factors = NULL,
  thinning = TRUE,
  print.level=1) {

  ## Create NBP object
  if (print.level > 0) cat("Create NBP data structure.\n");
  obj = list(counts = counts, lib.sizes = lib.sizes, grp.ids = grp.ids);

  ## Normalization

  ## Currently, no normalization is implemented, but user-specified
  ## normlization factors are supported.

  if (is.null(norm.factors)) {
    if (print.level>0) cat("No normalization is performed.\n");
    norm.factors = rep(1, length(lib.sizes));
  } else if (is.numeric(norm.factors)) {
    if (print.level>0) cat("Use specified normalization factors.\n");
    obj$norm.factors = norm.factors;
  }

  if (print.level > 1) {
    cat("  Normalization factors:"); 
    cat(norm.factors, sep=","); 
    cat("\n");
  }
  obj$eff.lib.sizes = lib.sizes * norm.factors;
    
  ## Library size adjustment
  if (thinning) {
    if (print.level > 0) cat("Thin the counts to make library sizes approximately equal.\n");
    res = thin.counts(counts, obj$eff.lib.sizes);
    obj$pseudo.counts = res$counts;
    obj$pseudo.lib.sizes = res$lib.sizes;
  } else {
    obj$pseudo.counts = counts;
    obj$pseudo.lib.sizes = rep(mean(obj$eff.lib.sizes), length(obj$eff.lib.sizes));
  }

  class(obj)="nbp";
  obj
}

## Estimating NBP dispersion parameters
estimate.disp = function(obj, method = "NBP", print.level=1, ...) {
  ## Estimate the dispersion parameters using counts from ALL
  ## treatment groups. In this step, we do not assume that genes have
  ## the same mean relative counts under under different treatments.

  ## Arguments:
  ##
  ##   obj: output from <prepare.nbp.obj>.
  ##
  ##   method: "NBP" (default) or "NB2", model for count variance (dispersion).
  ##
  ##   print.level: controls messages printed.

  if (method =="NBP") {
    nbp.est = nbp.mcle(obj$pseudo.counts, obj$pseudo.lib.sizes, obj$grp.ids,
      print.level=print.level, ...);
  } else if (method=="NB2") {
    nbp.est = nb2.mcle(obj$pseudo.counts, obj$pseudo.lib.sizes, obj$grp.ids,
      print.level=print.level, ...);
  } else {
    stop("Unknown method for estimating the dispersion parameters! Try NBP or NB2.");
  }

  res = obj;
  res$phi = nbp.est$phi;
  res$alpha = nbp.est$alpha;
  res$pie = nbp.est$pie;

  res;
}

## Perform exact test based on the NBP model.
nbp.test = function(counts, grp.ids, grp1, grp2,
  norm.factors = NULL,
##  thinning = TRUE,
  method.disp = "NBP",
  print.level = 1,
  ...) {

  ## Arguments
  ##
  ##   counts: an n x r matrix of counts with rows corresponding to
  ##   genes (exons, etc) and columns corresponding to libraries
  ##   (samples).
  ##
  ##   grp.ids: an r vector denoting the treatment groups. For each
  ##   gene, the mean counts within each group are assumed to be the
  ##   same.
  ##
  ##   grp1, grp2: the two groups to be compared
  ##
  ##   norm.factors: normalization factors
  ##
  ##   method.disp: NBP (default) or NB2, method for estimating the
  ##   dispersion parameters
  ## 
  ##   thinning: whethere the counts should be thinned so that the
  ##   library sizes (adjusted for normalization factors) are roughy
  ##   equal
  ##
  ##   ...: optional parameters to be passed to estimate.disp

  thinning = TRUE;

  ## Create an NBP object
  obj = prepare.nbp(counts, grp.ids, norm.factors = norm.factors, thinning=thinning,
    print.level=print.level);

  ## Estimating count variance
  obj = estimate.disp(obj, method = method.disp, print.level=print.level, ...);

  ## Testing for differential gene expression
  obj = exact.nb.test(obj, grp1, grp2, print.level=print.level);

  ## class(obj) = "nbp.obj";
  obj
}
