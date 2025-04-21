##Estimating Ne using neutral SNPs

library(data.table)
library(R.utils)
library(boot)
library(tidyverse)

# Read in allele frequency data

# MAFS
abol_mafs_neutral <- fread("abol_neutral_notrans.mafs.gz", header=TRUE)
cbol_mafs_neutral <- fread("cbol_neutral_notrans.mafs.gz", header=TRUE)
amta_mafs_neutral <- fread("amta_neutral_notrans.mafs.gz", header=TRUE)
cmta_mafs_neutral <- fread("cmta_neutral_notrans.mafs.gz", header=TRUE)

#Merge population comparisons

setnames(abol_mafs_neutral, c("knownEM", 'nInd'), c("freq1", 'nInd1'))
setnames(cbol_mafs_neutral, c("knownEM", 'nInd'), c("freq1", 'nInd1'))
setnames(amta_mafs_neutral, c("knownEM", 'nInd'), c("freq1", 'nInd1'))
setnames(cmta_mafs_neutral, c("knownEM", 'nInd'), c("freq1", 'nInd1'))


abol_cbol_mafs_neutral <- abol_mafs_neutral
abol_cbol_mafs_neutral$freq2 <- cbol_mafs_neutral$freq1
abol_cbol_mafs_neutral$NInd2 <- cbol_mafs_neutral$nInd1

amta_cmta_mafs_neutral <- amta_mafs_neutral
amta_cmta_mafs_neutral$freq2 <- cmta_mafs_neutral$freq1
amta_cmta_mafs_neutral$NInd2 <- cmta_mafs_neutral$nInd1

#Subset to maf>0.001

abol_cbol_mafs_001_neutral <- subset(abol_cbol_mafs_neutral, freq1 > 0.001 & freq2 > 0.001)
amta_cmta_mafs_001_neutral <- subset(amta_cmta_mafs_neutral, freq1 > 0.001 & freq2 > 0.001)

#Find where all populations overlap with maf>0.001

all_mafs_neutral <- abol_cbol_mafs_neutral
all_mafs_neutral$freq3 <- amta_mafs_neutral$freq1
all_mafs_neutral$nInd3 <- amta_mafs_neutral$nInd1
all_mafs_neutral$freq4 <- cmta_mafs_neutral$freq1
all_mafs_neutral$nInd4 <- cmta_mafs_neutral$nInd1

all_mafs_001_neutral <- subset(all_mafs_neutral, freq1 > 0.001 & freq2 > 0.001 & freq3 > 0.001 & freq4 > 0.001)
#703,178 SNPs meet this criteria

abol_cbol_allmafs_001_neutral <- select(all_mafs_001_neutral, -11:-14)
amta_cmta_allmafs_001_neutral <- select(all_mafs_001_neutral, -7:-10)
setnames(amta_cmta_allmafs_001_neutral, c("freq3", 'nInd3', "freq4", 'nInd4'), c("freq1", 'nInd1', "freq2", 'nInd2'))
abol_cbol_allmafs_001_neutral <- select(abol_cbol_allmafs_001_neutral, -3:-6)
amta_cmta_allmafs_001_neutral <- select(amta_cmta_allmafs_001_neutral, -3:-6)

setnames(abol_cbol_allmafs_001_neutral, c("NInd2"), c("nInd2"))

#Run calculations

### Jorde & Ryman/NeEstimator approach
# Jorde & Ryman 2007

# Ne in # diploid individuals
# based on NeEstimator manual v2.1

jrNe2 <- function(maf1, maf2, n1, n2, gen){
  Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles

  z <- (maf1+maf2)/2 # for the first allele
  z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
  Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
  Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual

  sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals

  S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
  S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
  Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
  return(gen/(2*Fsprime)) # calculation of Ne in # diploid individuals
}

abol_cbol_allmafs_001_neutral[, jrNe2(freq1, freq2, nInd1, nInd2, 114)] #2311.604
amta_cmta_allmafs_001_neutral[, jrNe2(freq1, freq2, nInd1, nInd2, 113)] #2230.299

#Bootstrap over loci to get CIs

## Jorde & Ryman Ne estimator, for boot() to use

jrNe2boot <- function(data, gen, indices){
  maf1 <- data$freq1[indices]
  maf2 <- data$freq2[indices]
  n1 <- data$nInd1[indices]
  n2 <- data$nInd2[indices]

  Fsnum <- (maf1-maf2)^2 + (1-maf1 - (1-maf2))^2 # the numerator, summing across the two alleles

  z <- (maf1+maf2)/2 # for the first allele
  z2 <- ((1-maf1)+(1-maf2))/2 # z for the 2nd allele
  Fsdenom <- z*(1-z) + z2*(1-z2) # the denominator of Fs, summing across the 2 alleles
  Fs <- sum(Fsnum)/sum(Fsdenom) # from NeEstimator calculations manual

  sl <- 2/(1/n1 + 1/n2) # harmonic mean sample size for each locus, in # individuals

  S <- length(maf1)*2/sum(2/sl) # harmonic mean sample size in # individuals, across loci and across both times. 2 alleles. Eq. 4.10 in NeEstimator v2.1 manual
  S2 <- length(maf2)*2/sum(2/n2) # harmonic mean sample size of 2nd sample in # individuals, across loci. all 2 alleles. See NeEstimator v2.1 below Eq. 4.13
  Fsprime <- (Fs*(1 - 1/(4*S)) - 1/S)/((1 + Fs/4)*(1 - 1/(2*S2))) # Eq. 4.13 in NeEstimator v2.1
  Ne <- gen/(2*Fsprime)
  if(Ne < 0) Ne <- Inf
  return(Ne) # calculation of Ne in # diploid individuals
}

boot_abol_cbol_neutral <- boot(data = abol_cbol_allmafs_001_neutral, statistic = jrNe2boot, R = 1000, gen = 114)
#original   bias    std. error
#t1* 2311.604 0.465257     9.32831
boot.ci(boot_abol_cbol_neutral, type='perc') #(2294, 2331)

boot_amta_cmta_neutral <- boot(data = amta_cmta_allmafs_001_neutral, statistic = jrNe2boot, R = 1000, gen = 113)
#original     bias    std. error
#t1* 2230.299 -0.6369695    11.75173
boot.ci(boot_amta_cmta_neutral, type='perc') #(2207, 2253)


