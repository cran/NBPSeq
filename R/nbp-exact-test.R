## Changlog:
##
##   0. Added nb2.test()
##
##   1. Added phi, alpha, pie, q.values to the output of nbp-test and
##   nb2-test.
## 

## Compute the tail conditional probablity of S1 | S1 + S2, assuming
## S1, S2 are both NB random varaibles with overdispersion parameters
## (theta1, p1) and (theta2, p2) respectively.
compute.tail.prob = function(s1, s2, theta1, theta2, p1, p2) {
  ## Arguments:
  ##   
  ##   s1, s2: two NBP random variables distribuated as NB(p, theta1)
  ##   and NB(p, theta2)
  ##
  ##   theta1, theta2: the theta parameters
  ##
  ##   p1, p2: the p parameters
  ##
  ## Value:
  ##
  ##   pval: the p-value (tail probablity)
  ##
  ## Details:
  ##
  ##   If S1, S2 are independent NB r.v.s with a commmon p parameter,
  ##   then S = S1 + S2 is also an NB random variable with theta
  ##   parameter thetas = theta1 + theta2.
  ##
  ##   This routine computes the probability of observing (S1, S2)
  ##   that are more extreme than (s1, s2) conditional on the value of
  ##   S1 + S2 = s1 + s2.


  ## Probability of the observed (s1, s2);
  pr.obs = dnbinom(s1, theta1, 1-p1) * dnbinom(s2, theta2, 1-p2);
  ## Probability of all pairs of (S1, S2) such that S1 + S2 = s.
  s = s1 + s2;
  pr = dnbinom(0:s, theta1, 1-p1) * dnbinom(s:0, theta2, 1-p2);

  ## Indices of the pairs that are more extreme than the observed pair. 
  ## WARNING: this line can be slow
  ## this line can be inaccurate? 
  id.extreme = (pr <= pr.obs);

  ## thetas = theta1 + theta2;
  ## ps = dnbinom(s, thetas, 1-p);

  pval = sum(pr[id.extreme])/sum(pr);
};


## Test <compute.tail.prob>
test.compute.tail.prob = function() {
  
  ## NB model parameters
  mu = 1:1000;
  ## mu = rep(1000, 1000);
  phi = 1.5; alpha = 1.5;
  ## phi = 2.0; alpha = 0.3;
  theta = mu^(2 - alpha)/phi;
  p = mu/(mu+theta);

  ## Simulate a sample from a true NB model
  n = length(mu);
  r = 6;
  y = rnbinom(n * r, theta, mu = mu);
  dim(y) = c(n,r);

  plot(mu, rowMeans(y));

  ## Compare the total count in the first column with that in columns
  ## 2:3 for each row.
  s1 = rowSums(y[,1:3]);
  s2 = rowSums(y[,4:6]);

  p.values = numeric(n);
  debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    p.values[i] = compute.tail.prob(s1[i], s2[i], theta[i]*3, theta[i]*3, p[i], p[i]);
  }
  hist(p.values);

  ## Compare the total count in the first column with that in columns
  ## 2:3 for each row.
  ##s1 = rowSums(y[,1:3]);
  s1 = y[,1];
  s2 = rowSums(y[,2:3]);

  p.values = numeric(n);
  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    p.values[i] = compute.tail.prob(s1[i], s2[i], theta[i], theta[i]*2, p[i], p[i]);
  }

  ## The histogram of the p.values should be roughly uniform.
  hist(p.values, prob=TRUE);

  invisible();
}

## Perform Robinson and Smyth exact NB test based on the fitted NBP model.

exact.nb.test = function(obj, grp1, grp2, print.level=1) {
  ## Arguments
  ##
  ##   obj: output of estimate.disp, a list with the following components:
  ##  
  ##     pseudo.counts: an n x r matrix of "normalized" counts with
  ##     rows corresponding to genes (exons, etc) and columns
  ##     corresponding to libraries (samples).
  ##
  ##     pseudo.lib.sizes: the library sizes.
  ##
  ##     grp.ids: an r vector denoting the treatment groups. For each
  ##     gene, the mean counts within each group are assumed to be the
  ##     same.
  ##
  ##     phi, alpha: dispersion parameters
  ##
  ##   grp1, grp2: the two groups to be compared
  ##
  ##   count.threshold: a test will be performed only when the total
  ##   row counts is greater than or equal to this value.
  ##
  ##   print.level: controls messages printed.
  ##
  ## Value:
  ##
  ##   The object <obj> from the input with the following additional components:
  ##
  ##   expression.levels: a data frame of grp1 , grp2 and pooled
  ##   expression levels
  ##
  ##   log.fc: base 2 log fold change in mean relative frequencies of reads
  ## 
  ##   p.values:
  ##
  ##   q.values:

  ## Details:
  ##
  ##   Currently, we assume that the dispersion parameters for the two
  ##   groups are the same.
  ##

  count.threshold=1;

  if (print.level>0) cat("Perform exact NB test for differential gene expression. \n");

  y = obj$pseudo.counts; # The pseudo data.
  m = obj$pseudo.lib.sizes; # The library sizes, assumed to be the same for libraries within each group
  phi = obj$phi;
  alpha = obj$alpha;
  grp.ids = obj$grp.ids;
  
  ## Estimate the relative mean counts assuming common means between group 1 and 2.
  if (print.level>1) cat("  Estimating pooled mean counts ... \n");

  ## MOM
  n = dim(y)[1];
  rel.counts = y / (matrix(1, n, 1) %*% m);
  grp12.ids = grp.ids %in% c(grp1, grp2);
  grp12.size = sum(grp12.ids);
  pie = rel.counts[,grp12.ids] %*% matrix(1/grp12.size, grp12.size, 1);

  ## Extract columns corresponding to the two treatment groups
  grp1.ids = grp.ids %in% grp1;
  grp2.ids = grp.ids %in% grp2;
  
  y1 = y[, grp1.ids];
  y2 = y[, grp2.ids];

  ## Compute the total counts for each gene within each treatment group
  s1 = apply(y1, 1, sum);
  s2 = apply(y2, 1, sum);
  ## s = apply(y, 1, sum);

  r1 = dim(y1)[2];
  r2 = dim(y2)[2];

  ## Compute the theta, p corresponding to the sums
  ## mu, theta, p, are n-vectors
  m1  = m[grp1.ids][1];
  m2  = m[grp2.ids][1];
  
  phi1 = phi2 = phi;
  ## phi1  = phi[grp1.ids][1];
  ## phi2  = phi[grp2.ids][1];

  mu1 = pie * m1;
  theta11 =  mu1^(2 - alpha) / phi1;
  p1 = mu1 / (mu1 + theta11);
  theta1 = theta11 * r1; 

  mu2 = pie * m2;
  theta21 =  mu2^(2 - alpha) / phi2;
  p2 = mu2 / (mu2 + theta21);
  theta2 = theta21 * r2;

  if (print.level>1) cat("  Performing the NB test for each gene ...\n");
  ## Performe the exact tests
  p.values = numeric(n);
  p.values[] = NA;

  threshold = max(count.threshold, 1);

  ## debug(compute.tail.prob);
  ## undebug(compute.tail.prob);
  for (i in 1:n) {
    if (s1[i] + s2[i] > threshold) {
      ## Perform the NBP test
      p.values[i] = compute.tail.prob(s1[i], s2[i], theta1[i], theta2[i], p1[i], p2[i]);
    }
  }

  log.fc = log2((s2/m2/r2)/(s1/m1/r1));

  expression.levels = cbind(s1/m1/r1, s2/m2/r2, pie);
  colnames(expression.levels) = c("grp1", "grp2", "pooled");

  ## Compute q-values
  q.values = rep(NA, length(p.values));
  id = !is.na(p.values);
  q.values[id] = qvalue(p.values[id])$qvalues;

  obj$grp1 = grp1;
  obj$grp2 = grp2;
  obj$pooled.pie = pie;
  obj$expression.levels = expression.levels;
  obj$log.fc = log.fc;
  obj$p.values = p.values;
  obj$q.values = q.values;

  obj
}

