## Functions for accessing elements in an "nbjObj"

## Arguments:
##
##   obj: a list, output from nbp-mcle or nbp-test
##

get.mean.hat = function(obj, grp.id) {
  r = length(obj$grp.ids);
  grp = (1:r)[obj$grp.ids == grp.id][1];
  mu = obj$pie[,grp] * obj$lib.sizes[grp];
}

get.var.hat = function(obj, grp.id) {
  r = length(obj$grp.ids);
  grp = (1:r)[obj$grp.ids == grp.id][1];
  mu = obj$pie[,grp] * obj$lib.sizes[grp];
  phi = obj$phi[grp];
  v = mu + phi * mu^obj$alpha;
}

get.phi = function(obj, grp.id) {
  r = length(obj$grp.ids);
  grp = (1:r)[obj$grp.ids == grp.id][1];
  phi = obj$phi[grp];
}

get.alpha = function(obj, grp.id) {
  alpha = obj$alpha;
}


## Retrieve nbp parameters for one of the treatment groups from an nbp object
get.nbp.pars = function(obj, grp.id) {
## Arguments:
##
##   obj: output form nbp.mcle
##
##   grp.id: the id of a treatment grp
##
## Value:
##
##   n: number of genes
##
##   r: number of replicates
##
##   lib.sizes:
##
##   phi, alpha, pie:w
##   
  r = length(obj$grp.ids);
  grp = (obj$grp.ids == grp.id);
  idx = (1:r)[grp][1];

  ## nbp.pars = new("nbpPars",
  nbp.pars = list(
       n = dim(obj$pseudo.counts)[1],
       r = sum(grp),
       grp.ids = obj$grp.ids[grp],
       lib.sizes = obj$lib.sizes[idx],
       pie = obj$pie[,idx],
       phi = obj$phi[idx],
       alpha = obj$alpha);
}

## Print nbp object
print.nbp.counts = function(obj, subset=1:10) {
## Arguments:
##
##   obj: an nbp obj
##
##   subset: 
##
## Value:
##
##   No return value.
  n = dim(obj$counts)[1];
  r = dim(obj$counts)[2];
  cat(sprintf("Number of rows: %d\n", n));
  cat(sprintf("Number of columns: %d\n", r));

  cat("Groups:", sprintf("%s", obj$grp.ids), "\n");

  cat("Counts:\n");
  print(obj$counts[subset,]);
  cat("...\n");
  cat("Library sizes:", obj$lib.sizes, "\n");
  cat("\n");
  
  cat("Pseudo counts:\n");
  print(obj$pseudo.counts[subset,]);
  cat("Pseudo library sizes:", obj$pseudo.lib.sizes, "\n");
  cat("\n");

  invisible();

}

## Print nbp.pars
print.nbp.pars = function(obj, subset=1:10) {
## Arguments:
##
##   obj: output form nbp.mcle
##
##   grp.id: the id of a treatment grp
##
## Value:
##
##   n: number of genes
##
##   r: number of replicates
##
##   lib.sizes:
##
##   phi, alpha, pie:w

  cat(sprintf("Number of rows: %d\n", obj$n));
  cat(sprintf("Number of columns: %d\n", obj$r));
  cat("Groups:", obj$grp.ids, "\n");
  cat("Library sizes:", obj$lib.sizes, "\n");
  cat("Adjusted library sizes:", obj$pseudo.lib.sizes, "\n");
  cat("Dispersion parameters:\n");
  cat("    phi:", obj$phi,"\n");
  cat("  alpha:", obj$alpha, "\n");
  
  cat("Mean relative frequencies:\n");
  print(obj$pie[subset,], digits=8);
  cat("...\n");
  cat("\n");
  
  invisible();
}


## Print NBP test results
print.nbp.test = function(obj, subset=1:10) {
  cat("Expression levels and measures of differential expression:\n");
  print(cbind(obj$expression.levels[subset,],
              log.fc = obj$log.fc[subset],
              p.value = obj$p.values[subset],
              q.value = obj$q.values[subset]));
  cat("\n");

  invisible();
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Print nbp object
##' @param x 
##' @param subset 
##' @param ... 
##' @return  NULL
##' @author Yanming Di
print.nbp = function(x, subset=1:10, ...) {
  obj=x;
  print.nbp.counts(obj, subset);
  print.nbp.pars(obj, subset);
  print.nbp.test(obj, subset);
  invisible();
}
