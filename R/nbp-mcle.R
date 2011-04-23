## Functions for fitting the NBP model for count variance.

## rm(list=ls());
## Computing the MLE estimator of the vector pi.
optimize.loglik.pie = function(pie=pie, phi=phi, alpha=alpha, y=y, m=m) {
  ## Arguments:
  ##
  ##   pie: an n vector, the current estimates of pi[i]'s.
  ##
  ##   phi, alpha: the dirpersion parameters.
  ##
  ##   y: the n x r matrix of counts.
  ##
  ##   m: an r vector of library sizes 
  ##
  ## Value:
  ##
  ##   pie.new: an n vector, the MLE estimator of pi.
  ##
  ## Details:
  ##
  ##   This function searches for the MLE of pi[i], the mean relative
  ##   gene counts, for each gene i.  When maximizing the likelihood,
  ##   the parameter used is actually log(pi[i]). This guarantees that
  ##   the estimated pi[i] is always positive. For each i, R function
  ##   optimize() is used to search the interval (log(pi[i].old) - 5,
  ##   log(max(y[i,])/m)) for a maximum point of the likelihood
  ##   function of log(pi[i]). The log likelihood of logpie can be
  ##   multimodal, we require that pie.hat < max(y)/m.
  ##
  ##   If the row sum is 0 for row i, the corresponding pi[i] will not
  ##   be updated.
  
  ## The log likelihood function of log(pie)
  l.logpie = function(logpie=logpie,  y=y) {
    ## Arguments:
    ##
    ##   logpie: a scalar
    ##
    ##   y: an r vector (one row in the original y matrix);
    ##
    ## Value:
    ##
    ##   l: a scalar, the log likelihood of logpie
    ##
    ## Details:
    ##
    ##   This is not a general purpose function for computing log
    ##   likelihood. In particular, it will not check values of logpie
    ##   and y to see whether they will lead to NaN values.
    mu =  exp(logpie) *  m;
    theta =  mu^(2 - alpha) / phi;
    l = sum(dnbinom(y, theta, mu=mu, log=TRUE));
  }

  n = dim(y)[1];
  r = dim(y)[2];
  t = y %*% matrix(1, r, 1);

  pie.new = pie;
  logpie = log(pie);

  ## The log likelihood of logpie can be multimodal, we require that pie.hat < max(y)/m.
  logpie.upper = log(apply(y / (matrix(1, n, 1) %*% m), 1, max));

  for (i in 1:n) {
    if (t[i] > 0) {
      ## Maximize the likelihood function for log(pi[i]);
      obj =  optimize(l.logpie, interval=c(logpie[i]-5, logpie.upper[i]), y=y[i,], maximum=TRUE);
      pie.new[i] = exp(obj$maximum);
    }
    ## cat(sprintf("%d, ", y[i,]));
    ## cat(sprintf("(mu.old, mu.new) = (%f, %f)\n", pie[i] * m, pie.new[i] * m));
  }
  pie.new;
}

## Compute the conditional log likelihood of (log(phi), alpha)
cloglik.logphi = function(logphi=logphi, alpha=alpha, pie=pie, y=y, m=m,
  n.grps=n.grps, grps=grps, grp.sizes=grp.sizes, mu.lower=mu.lower, mu.upper=mu.upper) {
  ## Arguments:
  ## 
  ##   logphi, alpha: the dispersion parameters
  ##
  ##   pie: an n x r matrix mean relative counts
  ##
  ##   y: the n x r matrix of gene counts
  ##
  ##   m: an r vector of library sizes
  ##
  ##   n.grps: number of treatment groups
  ##
  ##   grps: an n.grps x r group membership matrix
  ##
  ##   grp.sizes: an n.grps vector of group sizes
  ##
  ##   mu.lower, mu.upper: only genes with mean counts between
  ##   mu.lower and mu.upper will be used to compute the conditional
  ##   likilihood of the dispersion parameters
  ##
  ## Value:
  ##
  ##   the conditional log likelihood of (log(phi), alpha)
  ##
  ## Note:
  ##
  ##   We require (but do not check) that the library sizes within
  ##   each treatment group are the same

  l = double(n.grps);

  ## Compute the log conditional likelihood group by group
  for (i in 1:n.grps) {
    r = grp.sizes[i]; 

    ## We can only compute the conditional likelihood for groups with replicates
    if ( r > 1) {
      ## Compute the n vector of mean counts. The columns of pie are
      ## identical within a group, so we only need to take the first
      ## column.
      mu = pie[, grps[i,]][,1] * m[grps[i,]][1];

      theta =  mu^(2 - alpha) * exp(-logphi); # an n vector

      ## thresholding
      s = (mu > max(mu.lower, 1)) & (mu < mu.upper); # an n vector

      theta = theta[s];
      ## p = mu / (mu + theta);

      y0 = y[,grps[i,]][s,];
      
      l[i] = sum(lgamma(theta + y0)); # theta is n x 1, y0 is n x r 

      t = rowSums(y0);

      l[i] = l[i] + sum(lgamma(r * theta) - r * lgamma(theta) - lgamma(r * theta + t));
    }
  }

  sum(l);
}
## debug(cloglik.logphi);

## Maximize the conditional log likelihood wrt logphi = log(phi) for a
## fixed alpha.
optimize.cloglik.logphi = function(alpha, pie, y, m, n.grps, grps, grp.sizes,
  mu.lower, mu.upper, logphi.bounds) {
  ## Arguments:
  ##
  ##   alpha: the dispersion parameter
  ##
  ##   pie: an n x r matrix of the current estimates of the relative counts
  ##
  ##   y: the n x r matrix of counts
  ##
  ##   m: an r vector of library size
  ##
  ##   n.grps: number of treatment groups
  ##
  ##   grps: an n x r group membership matrix
  ##
  ##   grp.sizes: an n.grps vector of group sizes
  ##
  ##   mu.lower, mu.upper: only genes with mean counts between
  ##   mu.lower and mu.upper will be used to compute the conditional
  ##   likilihood of the dispersion parameters
  ##
  ## Value:
  ##
  ##   logphi: the log(phi) value that maximizes the conditional log
  ##   likelihood
  ##
  ##   l: the conditional log likelihood computed at the maximum point 
  ##
  ## Details:
  ##
  ##   This function uses the R function optimize() to search for a
  ##   maximum point of the conditional log likelihood for a fixed
  ##   alpha in the interval specified by logphi.bounds.
  
  obj = optimize(cloglik.logphi, interval=logphi.bounds,
    alpha=alpha, pie=pie, y=y, m=m,
    n.grps = n.grps, grps = grps, grp.sizes = grp.sizes,
    mu.lower = mu.lower, mu.upper=mu.upper,
    maximum="TRUE");

  ## print(obj);
  list(logphi=obj$maximum, l=obj$objective);
}

## Find values of (phi, alpha) that maximize the conditional log
## likelihood.
optimize.pcl = function(pie, y, m, n.grps, grps, grp.sizes, mu.lower, mu.upper) {
  ##  Arguments:
  ##
  ##    pie: an n x r matrix of the current estimate of relative counts
  ##
  ##    y: the n x r matrix of rnaseq counts
  ##
  ##    m: an r vector of library sizes
  ##
  ##    n.grps: number of treatment groups
  ##
  ##    grps: an n x r group membership matrix
  ##
  ##    grp.sizes: an n.grps vector of group sizes
  ##
  ##    mu.lower, mu.upper: only genes with mean counts between
  ##    mu.lower and mu.upper will be used to compute the conditional
  ##    likilihood of the dispersion parameters estimate the
  ##    dispersion parameters
  ##
  ##
  ##  Values:
  ##
  ##    tau: tau =(phi, alpha), the parameter values that maximize the
  ##    conditional log likelihood function.
  ##
  ##    l: the conditional log likelihood computed at the estimated tau.
  ##
  ##  Details:
  ##
  ##    This function maximizes the conditional likelihood by
  ##    maximizing the profile conditional likelihood
  ##
  ##       pcl(lapha) = max_logphi (l(logphi, alpha))
  ##
  ##    using the function opitmize().

  ## Bounds for alpha and logphi, these defaults should be good for
  ## most situations in practice.
  alpha.bounds = c(0, 3);
  logphi.bounds = c(-5, 5);

  ## The profile conditional likelihood of alpha
  pcl = function(alpha) {
    optimize.cloglik.logphi(alpha, pie, y, m, n.grps, grps, grp.sizes, mu.lower, mu.upper, logphi.bounds)$l;
  }

  ## Maximize the pcl of alpha
  obj = optimize(pcl, interval=alpha.bounds, maximum="TRUE");

  ## print(obj);

  ## Find the phi that maximizes the conditional log-likelihood when
  ## alpha is fixed at the optimal value.
  alpha.hat = obj$maximum;
  phi.hat = exp(optimize.cloglik.logphi(alpha.hat, pie, y, m, n.grps, grps, grp.sizes, mu.lower, mu.upper, logphi.bounds)$logphi);

  list(tau=c(phi.hat, alpha.hat), l=obj$objective);
}

## Estimate the parameters in the NBP model.
nbp.mcle = function(counts, lib.sizes, grp.ids,
  mu.lower=0, mu.upper=Inf, tol=0.005, maxit=0,
  print.level=1) {
  ## Arguments
  ##
  ##   counts: an n x r matrix of counts with rows corresponding to
  ##   genes (exons, etc) and columns corresponding to libraries
  ##   (samples).
  ##
  ##   grp.ids: an r vector denoting the treatment groups. For each
  ##   gene. The mean counts within each group are assumed to be the
  ##   same.
  ##
  ##   mu.lower, mu.upper: only genes with mean counts between
  ##   mu.lower and mu.upper will be used to compute the conditional
  ##   likilihood of the dispersion parameters estimate the dispersion
  ##   parameters
  ##
  ##   tol: the convergence tolerance
  ##
  ##   maxit: the maximum number of iterations
  ##
  ##   trace: trace=TRUE turns on printing of updated estimates at
  ##   each iteration
  ##
  ## Value:
  ##
  ##   counts, lib.sizes, grp.ids: same as input
  ##
  ##   phi, alpha: the estimated dispersion parameters
  ##
  ##   pie: a n x r matrix of estimated relative mean counts

  ## Details:
  ##
  ##   n is the number of genes, r is the number of replicates.  The
  ##   entry in row i and column j is the number of reads for gene i
  ##   in replicate j.
  ##
  ##   This function assumes that the dispersion parameters (phi,
  ##   alpha) for all treatment groups are the same.
  ## 
  ##   This function estimates a pi (theoretical relative frequency)
  ##   for each gene by maximum likelihood estimation and two
  ##   dispersion parameters of the negative binomial-P distribution,
  ##   phi and alpha, which are thought to be common to all genes.  In
  ##   this parameterization, the mean count is mu = pi*m where m is
  ##   the number of reads (library length) for a sample. The variance
  ##   is mu + phi*mu^alpha.

  ## CHECK y

  if (print.level > 0) cat("Estimate count variance assuming an NBP model. \n");

  if (print.level > 1)
    cat("  The NBP use two dispersion parameters (phi, alpha) to model extra-Poisson count variation.\n");
    
  if (is.matrix(counts)==FALSE) stop ("counts not a matrix");

  if (diff(range(lib.sizes)) / mean(lib.sizes) > 0.05)
    warning("The lib.sizes are not the same!");

  n = dim(counts)[1];		# number of genes 
  r = dim(counts)[2];		# number of libraries

  ## Identify the treatment groups
  gids = unique(grp.ids);  # unique group ids
  n.grps = length(gids);  # number of treatment groups
  grps = matrix(FALSE, n.grps, r);
  grp.sizes = integer(n.grps);
  for (i in  1:n.grps) {
    grps[i,] = (grp.ids == gids[i]);
    grp.sizes[i] = sum(grps[i,]);
  }

  if (max(grp.sizes) < 2) {
    stop ("No replicates in any treatment group. Dispersion parameters cannot be estimated.");
  }

  ## INITIAL ESTIMATES OF UNKNOWN PARAMETERS
  iter = 0;

  pie.current = matrix(0, n, r);

  ## Estimate pie by mean relative counts
  rel.counts = counts / (matrix(1, n, 1) %*% lib.sizes);

  for (i in 1:n.grps) {
    pie.current[, grps[i,]] = rel.counts[,grps[i,]] %*% matrix(1/grp.sizes[i], grp.sizes[i], 1);
    ## pie.current[, grps[i,]] = rowMeans(rel.counts[,grps[i,]]);
  }

  ## DEBUG
  ## print(counts[1:10,]);
  ## print(pie.current[1:10,] * (matrix(1,10,1) %*% lib.sizes));

  eps = .Machine$double.eps^(0.5); # A very small number

  ## For the initial estimate of (phi, alpha), all rows with non-zero
  ## totals will be used.
  tau = optimize.pcl(pie.current, counts, lib.sizes, n.grps, grps, grp.sizes, max(mu.lower, eps), mu.upper)$tau;

  phi.current = tau[1];
  alpha.current = tau[2];

  if (print.level>1) {
    cat(sprintf("  Iter  0: (phi, alpha) = (%6.3f, %6.3f).\n", phi.current, alpha.current));
  }

  diff = tol + 1;

  ## s0 = (t>0);
  ## START COMBINED CONDITIONAL LIKELIHOOD AND MAXIMUM LIKELIHOOD ITERATIONS
  while ((diff > tol) & (iter < maxit)) {
    iter = iter + 1;
    phi.old = phi.current;
    alpha.old = alpha.current;
    pie.old = pie.current;

    ## Update pi;
    ## if (trace==T) print("start maximum likelihood iterations")
    for (i in 1:n.grps) {
      if (grp.sizes[i] > 1) {
        pie.current[, grps[i,]] = optimize.loglik.pie(pie.current[, grps[i,]][,1], phi.current, alpha.current,
                     counts[,grps[i,]], lib.sizes[grps[i,]]);
      }
    }

    ## Update tau = (phi, alpha)
    ## if(trace==T) print("start conditional likelihood iterations")	

    ## For MLE estimator of pi, it is recommended use only rows with m
    ## * pi.hat > 1.
    tau = optimize.pcl(pie.current, counts, lib.sizes, n.grps, grps, grp.sizes, max(mu.lower, 1), mu.upper)$tau;
    phi.current = tau[1];
    alpha.current = tau[2];

    ## Compute the max change in this iteration. For pie and phi, we
    ## use the relative change. For alpha, we use the absolute change.
    phi.diff = abs(phi.current-phi.old)/phi.old;
    alpha.diff = abs(alpha.current-alpha.old);
    s0 =  (pie.old > 0);
    pie.diff = max(abs(pie.current[s0]-pie.old[s0])/pie.old[s0]);

    diff = max(phi.diff, alpha.diff, pie.diff);

    if (print.level>1) {
      cat(sprintf("  Iter %2d: (phi, alpha) = (%6.3f, %6.3f), delta(phi, alpha, max(mu_i)) = (%f, %f, %f)\n",
                  iter, phi.current, alpha.current,
                  phi.diff, alpha.diff, pie.diff));
    }
  } # End iterations

  list(counts=counts,
       lib.sizes=lib.sizes,
       grp.ids = grp.ids,
       pie=pie.current,
       phi=phi.current, 
       alpha=alpha.current);
} # End nbp.mcle


## Estimate the parameters in the NBP model, assuming alpha = 2.
nb2.mcle = function(counts, lib.sizes, grp.ids,
  mu.lower=0, mu.upper=Inf, alpha = 2, tol=0.005, maxit=0, 
  print.level = 1) {
  ## Arguments
  ##
  ##   counts: an n x r matrix of counts with rows corresponding to
  ##   genes (exons, etc) and columns corresponding to libraries
  ##   (samples).
  ##
  ##   grp.ids: an r vector denoting the treatment groups. For each
  ##   gene. The mean counts within each group are assumed to be the
  ##   same.
  ##
  ##   mu.lower, mu.upper: only genes with mean counts between
  ##   mu.lower and mu.upper will be used to compute the conditional
  ##   likilihood of the dispersion parameters estimate the dispersion
  ##   parameters
  ##
  ##   tol: the convergence tolerance
  ##
  ##   maxit: the maximum number of iterations
  ##
  ##   trace: trace=TRUE turns on printing of updated estimates at
  ##   each iteration
  ##
  ## Value:
  ##
  ##   pseudo.counts: counts after adjusting for unequal library sizes
  ##
  ##   lib.sizes: a r vector, the library sizes for the treatment
  ##   groups of the pseudo.counts. After adjustment, the library
  ##   sizes within each treatgroup are the same.
  ##
  ##   phi: the estimated dispersion parameter
  ##
  ##   alpha: the dispersion parameter alpha will be fixed at given
  ##   value (default is 2).
  ##
  ##   pie: a n x r matrix of estimated relative mean counts
  ##
  ## Details:
  ##
  ##   n is the number of genes, r is the number of replicates.  The
  ##   entry in row i and column j is the number of reads for gene i
  ##   in replicate j.
  ##
  ##   This function assumes that the dispersion parameters (phi,
  ##   alpha) for all treatment groups are the same.
  ## 
  ##   This function estimates a pi (theoretical relative frequency)
  ##   for each gene by maximum likelihood estimation and the
  ##   dispersion parameter of the negative binomial-P distribution,
  ##   phi (the other dispersion parameter alpha is fixed at the given
  ##   value, default is 2), which are thought to be common to all
  ##   genes.  In this parameterization, the mean count is mu = pi*m
  ##   where m is the number of reads (library length) for a
  ##   sample. The variance is mu + phi*mu^alpha.

  if (print.level > 0) cat("Estimate count variance assuming an NB2 model.\n");

  if (print.level > 1)
    cat("  The NB2 use one dispersion parameter phi to model extra-Poisson variation.\n");
    

  ## CHECK y
  if (is.matrix(counts)==FALSE) stop ("counts not a matrix");

  if (diff(range(lib.sizes)) / mean(lib.sizes) > 0.05)
    warning("The lib.sizes are not the same!");

  n = dim(counts)[1];		# number of genes 
  r = dim(counts)[2];		# number of libraries

  ## Identify the treatment groups
  gids = unique(grp.ids);  # unique group ids
  n.grps = length(gids);  # number of treatment groups
  grps = matrix(FALSE, n.grps, r);
  grp.sizes = integer(n.grps);
  for (i in  1:n.grps) {
    grps[i,] = (grp.ids == gids[i]);
    grp.sizes[i] = sum(grps[i,]);
  }

  if (max(grp.sizes) < 2) {
    stop ("No replicates in any treatment group. Dispersion parameters cannot be estimated.");
  }

  ## DEBUG
  ## print(grps); print(grp.sizes);

   
  ## t = apply(pseudo.counts,1,sum); 	# row sums

  ## INITIAL ESTIMATES OF UNKNOWN PARAMETERS
  iter = 0;

  pie.current = matrix(0, n, r);

  ## Estimate pie by mean relative counts
  rel.counts = counts / (matrix(1, n, 1) %*% lib.sizes);

  for (i in 1:n.grps) {
    pie.current[, grps[i,]] = rel.counts[,grps[i,]] %*% matrix(1/grp.sizes[i], grp.sizes[i], 1);
  }

  ## DEBUG
  ## print(counts[1:10,]);
  ## print(pie.current[1:10,] * (matrix(1,10,1) %*% lib.sizes));

  eps = .Machine$double.eps^(0.5); # A very small number
  logphi.bounds = c(-5, 5);

  ## For the initial estimate of phi, all rows with non-zero
  ## totals will be used.
  obj = optimize(cloglik.logphi, alpha = alpha, interval=logphi.bounds,
    pie=pie.current, y = counts, m = lib.sizes,
    n.grps = n.grps, grps = grps, grp.sizes = grp.sizes,
    mu.lower = max(mu.lower,eps), mu.upper=mu.upper,
    maximum="TRUE");

  phi.current = exp(obj$maximum);

  if (print.level>1) {
    cat(sprintf("  Iter  0: (phi, alpha) = (%6.3f, %6.3f).\n", phi.current, alpha));
  }

  diff = tol + 1;

  ## s0 = (t>0);
  ## START COMBINED CONDITIONAL LIKELIHOOD AND MAXIMUM LIKELIHOOD ITERATIONS
  while ((diff > tol) & (iter < maxit)) {
    iter = iter + 1;
    phi.old = phi.current;
    pie.old = pie.current;

    ## Update pi;
    ## if (trace==T) print("start maximum likelihood iterations")
    for (i in 1:n.grps) {
      if (grp.sizes[i] > 1) {
        pie.current[, grps[i,]] = optimize.loglik.pie(pie.current[, grps[i,]][,1], phi.current, alpha,
                     counts[,grps[i,]], lib.sizes[grps[i,]]);
      }
    }

    ## Update phi;
    ## if(trace==T) print("start conditional likelihood iterations")	

    ## For MLE estimator of pi, it is recommended use only rows with m
    ## * pi.hat > 1.
    obj = optimize(cloglik.logphi, alpha = alpha, interval=logphi.bounds,
      pie=pie.current, y = counts, m = lib.sizes,
      n.grps = n.grps, grps = grps, grp.sizes = grp.sizes,
      mu.lower = max(mu.lower,1), mu.upper=mu.upper,
      maximum="TRUE");

    phi.current = exp(obj$maximum);

    ## Compute the max change in this iteration. For pie and phi, we
    ## use the relative change. For alpha, we use the absolute change.
    phi.diff = abs(phi.current-phi.old)/phi.old;
    s0 =  (pie.old > 0);
    pie.diff = max(abs(pie.current[s0]-pie.old[s0])/pie.old[s0]);

    diff = max(phi.diff, pie.diff);

    if (print.level>1) {
      cat(sprintf("  Iter %2d: (phi, alpha) = (%6.3f, %6.3f), delta(phi, max(mu_i)) = (%f, %f)\n",
                  iter, phi.current, alpha, phi.diff, pie.diff));
    }
  } # End iterations

  list(counts=counts,
       grp.ids = grp.ids,
       lib.sizes=lib.sizes,
       pie=pie.current,
       phi=phi.current,
       alpha=alpha);
} # End nb2.mcle
