## Some diagnostic plots for NBP models: mean-variance plots on log-log scale, MA-plot

##  Plot log fold change vs log mean relative counts (MA plots).
##  Highlight differentially expressed genes.  Differentially
##  expressed genes are those with smallest p-values.  The user can
##  specify either the p-value cutoff or the number of genes to be
##  claimed as significant.
ma.plot = function(test.out,
  top = NULL,
  q.cutoff = NULL,
  p.cutoff = NULL,
  col.sig = "magenta",
  main = "MA Plot",
  ...
  ) {
  ## Arguments:
  ##
  ##   test.out: output of nbp.test or nb2.test
  ##
  ##   top: the number of genes to be claimed as differentially expressed
  ##
  ##   p.cutoff: p-value cutoff
  ##
  ##   q.cutoff: q-value cutoff
  ##
  ##   title: title of the MA plot
  ##
  ## Value
  ##
  ##    id: indices of top genes.

  if (all(is.null(c(top, q.cutoff, p.cutoff)))) {
    stop("Need value for top, p.cutoff or q.cutoff.");
  }

  if (is.null(top)) {
    if (!is.null(q.cutoff)) {
      top = sum(test.out$q.values <= q.cutoff, na.rm=TRUE);
    } else {
      top = sum(test.out$p.values <= p.cutoff, na.rm=TRUE);
    }
  }

  log2.pie = log2(test.out$pooled.pie);
  log2.fc = log2(exp(test.out$log.fc));

  ## p-value plot
  smart.plot(log2.pie, log2.fc,
       xlab = "log2 mean relative count",
       ylab = "log2 fold change",
       main = main,
       pch = "+", ...
             );

  id = NULL;
  if (top > 0) {
    ## List the most significant test results
    id = order(test.out$p.values)[1:top];

    ## Highlight differentially expressed genes
    points(log2.pie[id], log2.fc[id], col=col.sig, pch="+");
  }

  invisible(id);
}


## Plot mean variance 
mv.plot = function(counts,
  xlab = "mean", ylab = "variance",
  main = "variance vs mean",
  log ="xy",
  ...
  ) {
  mu = rowMeans(counts);
  v = apply(counts, 1, var);

  id = (mu>0) & (v>0);
  mu = mu[id];
  v = v[id];

  smart.plot(mu, v, xlab=xlab, ylab=ylab, main=main, log=log, ...);

  invisible();
}

## Highlight a subset of points on the mean-variance plot
mv.points = function(counts, subset, ...) {
  mu = rowMeans(counts);
  v = apply(counts, 1, var);
  points(mu[subset], v[subset], ...);
  invisible();
}

## Overlay an estimated mean-variance line 
mv.line = function(mu, v, ...) {
  id = order(mu);
  mu = mu[id];
  v = v[id];

  if (length(mu) > 1000) {
    id = c(seq(1, length(mu)-1, length=1000), length(mu));
  }
  lines(mu[id], v[id],  ...);

  invisible();
}

## Fit polynomial regression to (log(mean), log(var))
mv.line.fitted = function(obj, ...) {

  ## obj: a fitted curve, with components mu and v 
  mv.line(obj$mu, obj$v, ...);
  invisible();
}

mv.line.nbp = function(nbp.obj, grp.id, ...) {
  ## nbp.obj: output from nbp-mcle.
  mv.line(get.mean.hat(nbp.obj, grp.id), get.var.hat(nbp.obj, grp.id), ...);
  invisible();
}

## Plot variance-to-mean ratio
vmr.plot = function(counts, alpha = 1, 
         xlab = "mean",
         ylab = "variance",
         main = "variance/mean vs mean",
         log = "xy",
  ...) {

  mu = apply(counts, 1, mean);
  v = apply(counts, 1, var);

  id = (mu>0) & (v>0);
  mu = mu[id];
  v = v[id];
  vmr = v/mu^alpha;

  smart.plot(mu, vmr, xlab=xlab, ylab=ylab, main=main, log=log, ...);

  invisible();
}

## Overlay an esimtated VMR line 
vmr.line = function(mu, v, alpha=1, ...) {

  id = order(mu);
  mu = mu[id];
  v = v[id];
  vmr = v/mu^alpha;

  if (length(mu) > 1000) {
    id = c(seq(1, length(mu)-1, length=1000), length(mu));
  }
  lines(mu[id], vmr[id], ...);

  invisible();
}

vmr.line.fitted = function(obj, alpha=1, ...) {
  ## obj: a fitted curve, with components mu and v 
  vmr.line(obj$mu, obj$v, alpha=alpha, ...);
  invisible();
}

vmr.line.nbp = function(nbp.obj, grp.id, alpha=1, ...) {
  vmr.line(get.mean.hat(nbp.obj,1), get.var.hat(nbp.obj, 1), alpha = alpha, ...);
  invisible();
}

## Plot phi.hat assuming NB2
phi.plot = function(counts, alpha = 2, 
         xlab = "mean",
         ylab = "phi.hat",
         main = "phi.hat vs mean",
         log="xy",
  ...) {

  v = apply(counts, 1, var);
  mu = apply(counts, 1, mean);
  phi = (v - mu) / mu^alpha;

  id = (phi > 0) & (mu>0);
  mu  = mu[id];
  phi = phi[id];

  smart.plot(mu, phi, xlab=xlab, ylab=ylab, main=main, log=log, ...);

  invisible(phi);
}

phi.line = function(mu, v, alpha=2, ...) {
  id = order(mu);
  mu = mu[id];
  v = v[id];

  phi = (v - mu) / mu^alpha;

  id = (phi > 0) & (mu>0);
  mu  = mu[id];
  phi = phi[id];

  if (length(mu) > 1000) {
    id = c(seq(1, length(mu)-1, length=1000), length(mu));
  }
  
  lines(mu[id], phi[id], ...);
  invisible();
}

phi.line.fitted = function(obj, alpha=2, ...) {
  ## obj: a fitted curve, with components mu and v 
  phi.line(obj$mu, obj$v, alpha=alpha, ...);
  invisible();
}

phi.line.nbp = function(nbp.obj, grp.id, alpha=2, ...) {
  phi.line(get.mean.hat(nbp.obj,1), get.var.hat(nbp.obj, 1), alpha = alpha, ...);
  invisible();
}

## 
plot.counts = function(counts, fitted.lines,
  plots = c("mv", "vmr", "phi"),
  cols = 1:length(fitted.lines), lwds = 2, ltys = 1,
  alpha = 2, 
  ...) {
  ##
  ## counts:
  ##
  ## fitted.lines:  a list of fitted lines

  n.lines = length(fitted.lines);

  cols = rep(cols, n.lines)[1:n.lines];
  lwds = rep(lwds, n.lines)[1:n.lines];
  ltys = rep(ltys, n.lines)[1:n.lines];

  if ("mv" %in% plots) {
    ## mean-variance plot
    mv.plot(counts, ...);

    if (n.lines > 0) {
      for (i in 1:n.lines) {
        mv.line(fitted.lines[[i]]$mu, fitted.lines[[i]]$v,
                col=cols[i], lty=ltys[i], lwd=lwds[i]);
      }
    }
  }

  if ("vmr" %in% plots) {
    ## mean-variance plot
    vmr.plot(counts, ...);

    if (n.lines > 0) {
      for (i in 1:n.lines) {
        vmr.line(fitted.lines[[i]]$mu, fitted.lines[[i]]$v,
                 col=cols[i], lty=ltys[i], lwd=lwds[i]);
      }
    }
  }

  if ("phi" %in% plots) {
    ## mean-variance plot
    phi.plot(counts, alpha=alpha, ...);

    if (n.lines>0) {
      for (i in 1:n.lines) {
        phi.line(fitted.lines[[i]]$mu, fitted.lines[[i]]$v,
                 col=cols[i], lty=ltys[i], lwd=lwds[i], alpha=alpha);
      }
    }
  }

  invisible();
}


plot.nbp = function(x, ...){
  ## Arguments:
  ##
  ##   obj: an NBP object
  obj = x;

  grp1 = obj$grp1;

  par.old= par(mfrow=c(2,2));

  ## Boxplot
  boxplot(log(obj$pseudo.counts+1), main="Boxplots of log(counts+1)");

  ## MA plot
  ma.plot(obj, q.cutoff=0.05, main="MA plot");

  ## MV plot
  mv.plot(obj$pseudo.counts, clip = 32,
          main="Mean-Variance Plot",
          xlab="Average Gene Count",
          ylab="Estimated Variance of Gene Counts"
          );
  mv.line.nbp(obj, grp1, col="magenta");
  ## mv.line.nbp(obj, grp2);

  ## Dispersion plot
  phi.plot(obj$pseudo.counts, clip = 32,
           main="Mean-Dispersion Plot",
           xlab="Average Gene Count",
           ylab="Estimated NB Dispersion"
           );
  phi.line.nbp(obj, grp1, col="magenta");
  ## phi.line.nbp(obj, grp1);

  par(par.old);

  invisible();
}

