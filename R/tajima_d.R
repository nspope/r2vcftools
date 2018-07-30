
.taj_d <- function(afs) # tajima's d
{
  n <- length(afs) + 1
  Sn <- sum(afs)
  theta_pi <- 2/(n*(n-1)) * sum( (1:(n-1)) * (n - 1:(n-1)) * afs)
  theta_w <- sum(afs)/sum(1/(1:(n-1)))
  a1 <- sum(1/(1:(n-1)))
  a2 <- sum(1/(1:(n-1))^2)
  b1 <- (n + 1)/(3*(n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9*n*(n-1))
  c1 <- (b1 - 1/a1)
  c2 <- (b2 + a2/a1^2 - (n+2)/(a1*n))
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  c(theta_pi - theta_w, e1*Sn + e2*Sn*(Sn - 1))# num, denom^2
}

.get_stats <- function(geno)
{
  # stats needed for VCFtools method of calculating D
  cnt <- apply(geno, 2, function(x) sum(x[x>0]))
  nonmiss <- apply(geno, 2, function(x) 2*sum(x!=-1))
  p <- ifelse(nonmiss > 0, cnt/nonmiss, 0)
  He <- sum(p*(1-p))
  N <- nrow(geno)*2
  Sn <- sum(p > 0)
  return (list(He=He, N=N, Sn=Sn))
}

.taj_d_vcftools <- function(stats) # tajima's D, ala vcftools
{
  He <- stats$He
  N <- stats$N
  Sn <- stats$Sn
  a1 <- sum(1/(1:(N-1)))
  a2 <- sum(1/(1:(N-1))^2)
  b1 <- (N + 1)/(3*(N - 1))
  b2 <- 2 * (N^2 + N + 3)/(9*N*(N-1))
  c1 <- (b1 - 1/a1)
  c2 <- (b2 + a2/a1^2 - (N+2)/(a1*N))
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  theta_pi <- 2 * He * N/(N-1)
  theta_w <- Sn/a1
  (theta_pi - theta_w)/sqrt(e1*Sn + e2*Sn*(Sn - 1)) 
}

.sim_stats <- function(geno, maf=0, ploidy=2)
{
  # simulate given null model, stuff needed for VCFtools method of calculating D
  N <- nrow(geno)*ploidy

  # get probability of allele count given maf filter, sample size
  pr <- c(0, 1/(1:(N-1)), 0)
  maf <- ceiling(N*maf)
  if (maf > 0)
    pr <- c(rep(0, maf-1), 1/(maf:(N-maf)), rep(0, maf-1)) 
  else
    pr <- c(1/(1:(N-1)))

  # simulate data respecting missing data pattern
  nonmiss <- apply(geno, 2, function(x) 2*sum(x!=-1))
  cnt <- rep(0, ncol(geno))
  for (i in unique(nonmiss))
    cnt[nonmiss == i] <- sample(1:(i-1), sum(nonmiss == i), prob=pr[1:(i-1)], replace=TRUE)

  p <- ifelse(nonmiss > 0, cnt/nonmiss, 0)
  He <- sum(p*(1-p))
  Sn <- sum(p > 0)
  return (list(He=He, N=N, Sn=Sn))
}

.sim_afs <- function(N, Sn, maf=0, ploidy=2) 
{
  # simulate afs under neutral model with missing data and a MAF cutoff -- *at least* MAF chromosomes must have minor allele
  pr <- c(0, 1/(1:(N-1)), 0)
  maf <- ceiling(N*maf)
  if (maf > 0)
    pr <- c(rep(0, maf-1), 1/(maf:(N-maf)), rep(0, maf-1)) 
  else
    pr <- c(1/(1:(N-1)))
  rmultinom(1, size=Sn, prob=pr)
}

.TajimaD <- function(geno, ploidy=2)
{
  # my way of calculating D, just an aggregate estimator across missing data patterns
  mis_pat <- apply(geno,2,function(x) ploidy*sum(x!=-1))
  sums <- apply(geno,2,function(x) sum(x[x!=-1]))

  afs <- list()
  for (i in unique(mis_pat))
  {
    af <- table(sums[mis_pat==i])
    afs[[as.character(i)]] <- rep(0, i-1)
    afs[[as.character(i)]][as.numeric(names(af))] <- af[]
  }

  tmp <- sapply(afs, .taj_d)
  tmp <- apply(tmp, 1, mean)
  tmp[1]/sqrt(tmp[2])
}

#' Tajima's D
#'
#' @description
#' Calculate a genome-wide estimate of Tajima's D, and optionally use simulation from the neutral model to correct for bias due to a minor-allele-frequency filter and perform a significance test.
#'
#' @param object an object of class "vcfLink"
#' @param maf the minimum allele frequency allowed in the null model simulations, as a proportion; SNPs with frequency < MAF are excluded from the simulations. The data is also thinned using this criterion prior to calculating D.
#' @param nboot the number of bootstrap replicates used to generate a confidence interval around D; if 0, no bootstrap is performed
#' @param nsim the number of simulations from the null model, used to correct bias due to the MAF and to test the significance of D; if 0, no correction is performed
#' @param ploidy the ploidy of the organism
#' @param use_vcftools_D if TRUE, use the method for calculating D implemented in vcftools, which is biased when there is missing data. If FALSE (default), use an unbiased (and faster) method
#' @rdname TajimaD
#' @export
TajimaD <- function(object, maf=0, nboot=0, nsim=1000, ploidy=2, use_vcftools_D = FALSE)
{
  if (maf>0)
    object <- Filter(object, filterOptions(maf=maf), verbose=FALSE)

  geno <- GenotypeMatrix(object, verbose=FALSE)

  min_count <- min(apply(geno, 2, function(x) sum(x[x > -1])))
  min_maf <- min(apply(geno, 2, function(x) ceiling(maf*ploidy*sum(x > -1))))

  if (min_count < min_maf) stop("'maf' does not agree with data: minimum count in data was ", min_count, " while minimum count allowed by maf filter was", min_maf)
  if (min_count > min_maf) warning("'maf' may not agree with data: minimum count in data was ", min_count, " while minimum count allowed by maf filter was", min_maf)

  if (!use_vcftools_D)
  {

  mis_pat <- table(apply(geno,2,function(x) ploidy*sum(x!=-1)))
  sims <- rep(NA, nsim)
  if (nsim > 0)
    for (i in 1:nsim)
    {
      if (i %% 100 == 0) cat("Simulation: iter", i, "\n")
      sim_num <- rep(0, length(mis_pat)) # numerator of D
      sim_den <- rep(0, length(mis_pat)) # squared denominator of D
      for (j in 1:length(mis_pat))
      {
        tmp <- .taj_d(.sim_afs(as.numeric(names(mis_pat)[j]), mis_pat[j], maf=maf))
        sim_num[j] <- tmp[1]
        sim_den[j] <- tmp[2]
      }
      
      # for each missing data pattern, we've calculated
      # \theta_\pi and \theta_w and the variance thereof.
      # To get an aggregate estimator, we take the mean.
      sims[i] <- mean(sim_num)/sqrt(mean(sim_den))
    }

  boot <- rep(NA, nboot)
  if (nboot > 0)
    for (i in 1:nboot)
    {
      if (i %% 100 == 0) cat("Bootstrap: iter", i, "\n")
      boot[i] <- .TajimaD(geno[,sample(1:ncol(geno), ncol(geno), replace=T)])
    }
  obs <- .TajimaD(geno)

  }
  else
  {

  sims <- rep(NA, nsim)
  if (nsim > 0)
    for (i in 1:nsim)
    {
      if (i %% 100 == 0) cat("Simulation: iter", i, "\n")
      sims[i] <- .taj_d_vcftools(.sim_stats(geno, maf))
    }

  boot <- rep(NA, nboot)
  if (nboot > 0)
    for (i in 1:nboot)
    {
      if (i %% 100 == 0) cat("Bootstrap: iter", i, "\n")
      boot[i] <- .taj_d_vcftools(.get_stats(geno[,sample(1:ncol(geno), ncol(geno), replace=T)]))
    }

  obs <- .taj_d_vcftools(.get_stats(geno))

  }

  bias <- if(nsim > 0) mean(sims) else 0

  sims_corr <- (sims - bias)
  boot_corr <- (boot - bias)
  obs_corr <- (obs - bias)

  p <- sum(abs(sims_corr) >= abs(obs_corr - bias))/nsim

  list (results = list(
        "Tajima D, raw" = obs,
        "Tajima D, bias-corrected" = obs_corr,
        "Tajima D, bias-corrected 95CI" = quantile(boot_corr, c(0.025, 0.975)),
        "Pr(abs(D) > abs(D_observed)|Null)" = ifelse(p, p, paste("<", 1/nsim)),
        "Bias" = bias),
        simulations = list(
        "Null" = sims,
        "Null, bias-corrected" = sims_corr,
        "Bootstrap" = boot,
        "Bootstrap, bias-corrected" = boot_corr))
}

#TajimaD_test_vcftools <- function(object, maf=0, nboot=1000, nsim=1000, null_quantiles=seq(0,1,0.1))
#{
#  geno <- GenotypeMatrix(object)
#  
#  sims <- rep(NA, nsim)
#  if (nsim > 0)
#    for (i in 1:nsim)
#    {
#      if (i %% 100 == 0) cat("Simulation: iter", i, "\n")
#      sims[i] <- taj_d_vcftools(sim_stats(geno, maf))
#    }
#
#  boot <- rep(NA, nboot)
#  if (nboot > 0)
#    for (i in 1:nboot)
#    {
#      if (i %% 100 == 0) cat("Bootstrap: iter", i, "\n")
#      boot[i] <- taj_d_vcftools(get_stats(geno[,sample(1:ncol(geno), ncol(geno), replace=T)]))
#    }
#
#  obs <- taj_d_vcftools(get_stats(geno))
#  bias <- mean(sims)
#  sims <- (sims - bias)
#  boot <- (boot - bias)
#  obs_corr <- (obs - bias)
#
#  p <- sum(abs(sims) >= abs(obs_corr))/nsim
#
#  list ("Tajima D, raw" = obs,
#        "Bias" = bias,
#        "Tajima D, bias-corrected" = obs_corr,
#        "Tajima D, bias-corrected 95CI" = quantile(boot, c(0.025, 0.975)),
#        "Null model quantiles" = quantile(sims, null_quantiles),
#        "Pr(abs(D) > abs(D_observed)|Null)" = ifelse(p, p, paste("<", 1/nsim)))
#}
#
