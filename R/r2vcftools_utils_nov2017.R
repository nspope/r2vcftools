#####################################################################################################
######################################### r2vcftools utilities ######################################
#####################################################################################################

# Rodolfo Jaffe, Eder Lanes, et al

# To run this script you need to install VCFtools
# GitHub: https://github.com/vcftools/vcftools
# Manual: https://vcftools.github.io/man_latest.html

## Install the r2vcftools library from GitHub
#install.packages("devtools")
#devtools::install_github("Bioconductor-mirror/LEA", force=T)
#devtools::install_github("nspope/r2vcftools", force=T)

### Summary stats
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

####CI
CI <- function(df, var) {
  error <- qt(0.975,df=length(df[, var])-1)*sd(df[, var])/sqrt(length(df[, var]))
  M <-  mean(df[, var])
  left <- M-error
  right <- M+error
  res <- data.frame(species=sp, mean=M, low=left, high=right)
  print(res)
}

### Genetic diversity measures
GenDiv <- function(vcf){
  
  HE <- Query(vcf, type="het")
  
    HE$HO <- (HE$N_SITES-HE$O.HOM.)/(HE$N_SITES) ## Observed heterozygosity (HO)
     error_ho <- qt(0.975,df=length(HE$HO)-1)*sd(HE$HO)/sqrt(length(HE$HO))
     ho <- mean(HE$HO)
     left_ho <- ho-error_ho
     right_ho <- ho+error_ho
    HO.df <- data.frame(HO=ho, lowHO=left_ho, upHO=right_ho)
    
    HE$HE <- (HE$N_SITES-HE$E.HOM.)/(HE$N_SITES) ## Expected heterozygosity (HE)
     error_he <- qt(0.975,df=length(HE$HE)-1)*sd(HE$HE)/sqrt(length(HE$HE))
      he <- mean(HE$HE)
      left_he <- he-error_he
      right_he <- he+error_he
    HE.df <- data.frame(HE=he, lowHE=left_he, upHE=right_he)
  
    error_f <- qt(0.975,df=length(HE$F)-1)*sd(HE$F)/sqrt(length(HE$F)) ## Inbreeding (F)
    f <- mean(HE$F)
    left_f <- f-error_f
    right_f <- f+error_f
  F.df <- data.frame(F=f, lowF=left_f, upF=right_f)
  
PI <- Query(vcf, type="site-pi") ## Nucleotide diversity (PI)
  error_pi <- qt(0.975,df=length(PI$PI)-1)*sd(PI$PI)/sqrt(length(PI$PI))
  pi <- mean(PI$PI)
  left_pi <- pi-error_pi
  right_pi <- pi+error_pi
PI.df <- data.frame(PI=pi, lowPI=left_pi, upPI=right_pi)

RES <- cbind(HO.df,  HE.df,  F.df,  PI.df)
 return(RES)
}

### Select best run (lowest cross-entropy) from snmf projects
Best.run <- function(nrep, optimalK, p1, p2, p3, p4){
  ce1 = cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, ", run = ", Best_run)) 
}

### Fst function (LEA)
fst = function(project, run = 1, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t}
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}

### Genomic Inflation Factor (lambda)
GIF <- function(project, run, K, fst.values){
  n = dim(Q(project, K, run))[1]
  fst.values[fst.values<0] = 0.000001
  z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
  lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
  print(lambda)
  return(lambda)
}

###Adjusted p-values
Apvals <- function(project, run, K, fst.values, lambda) {
  n = dim(Q(project, K, run))[1]
  z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
  AP = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
  return(AP)
}
                      
## FDR control: Benjamini-Hochberg at level q
candidates <- function(alpha=0.05, adj.p.values){ ## alpha is significance level
  L = length(adj.p.values) ## Number of loci
  w = which(sort(adj.p.values) < alpha * (1:L) / L)
  candidates = order(adj.p.values)[w]
  print(length(candidates))
  return(candidates)
}

###Manhatan plot
ManPlot <- function(adj.p.values, candidates){
  plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
  points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")
}

