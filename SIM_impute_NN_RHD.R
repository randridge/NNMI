####################
# NN MI Simulation
# Impute with NN or RHD + ABB (original)
# Author: Rebecca Andridge
# Last Modified: 5/3/2021
####################
rm(list=ls())

require(doParallel)
require(foreach)
require(doRNG)
require(tidyverse)
require(survey)
require(rsimsum)

# For testing
# N <- 1000
# SAMPSIZE <- 100
# RHO <- 0.9
# MISSMECH <- "mar50neg"  # mcar25, mcar50, mar25pos, mar25neg, mar50pos mar50neg
# K <- 1   # size of donor pool (0 = use all respondents)
# D <- 20  # number of imputations
# TYPE <- "pps"  # "srs" or "pps"

impute <- function(N, SAMPSIZE, RHO, MISSMECH, K, D, TYPE)
{
  # Load data
  rho <- RHO*10
  if (TYPE=="srs") {
    load(paste(paste("./data_samp/samples_srs_norm_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
    load(paste(paste("./data_samp/samples_srs_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
    load(paste(paste("./data_samp/samples_srs_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
    # Sample weights (N/n since SRS)
    samples_norm$w <- samples_unif$w <- samples_logn$w <- N/SAMPSIZE
  } else if (TYPE=="pps") {
    load(paste(paste("./data_samp/samples_pps_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
    load(paste(paste("./data_samp/samples_pps_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
    # Sample weights
    samples_unif$w <- 1/samples_unif$pik
    samples_logn$w <- 1/samples_logn$pik
  }
  # Induce missingness for specified missingness mechanism
  if (TYPE=="srs")
  {
    miss <- samples_norm[[MISSMECH]]
    samples_norm$y[miss] <- NA
  }
  miss <- samples_unif[[MISSMECH]]
  samples_unif$y[miss] <- NA
  miss <- samples_logn[[MISSMECH]]
  samples_logn$y[miss] <- NA
  a <- Sys.time()
  # Parallel processing
  NREPS <- max(samples_unif$sample)
  cl <- makeCluster(8, setup_strategy="sequential")
  registerDoParallel(cl)
  # Looping through simulated data sets
  post <- foreach(i=1:NREPS, .combine=rbind, .packages=c("tidyverse","mice","survey")) %dorng% {
    # Functions
    source("./functions/nn_impute.R")
    source("./functions/miStats_srs.R")
    source("./functions/miStats_pps.R")
    # Get samples
    if (TYPE=="srs") samp_norm <- as_tibble(subset(samples_norm, sample==i))
    samp_unif <- as_tibble(subset(samples_unif, sample==i))
    samp_logn <- as_tibble(subset(samples_logn, sample==i))
    # Multiply impute using nearest neighbor/variant
    if (TYPE=="srs") yimp_norm <- nn_impute(samp_norm$y, samp_norm$x, K, D)
    yimp_unif <- nn_impute(samp_unif$y, samp_unif$x, K, D)
    yimp_logn <- nn_impute(samp_logn$y, samp_logn$x, K, D)
    # Apply Rubin's rules
    if (TYPE=="srs") {
      postmi_norm <- miStats_srs(yimp_norm, 1-SAMPSIZE/N)
      postmi_unif <- miStats_srs(yimp_unif, 1-SAMPSIZE/N)
      postmi_logn <- miStats_srs(yimp_logn, 1-SAMPSIZE/N)
      res <- cbind(varmethod=0,distnum=1:3, rbind(postmi_norm, postmi_unif, postmi_logn))
    } else if (TYPE=="pps") {
      postmi_unif <- miStats_pps(yimp_unif, samp_unif$pik)
      postmi_logn <- miStats_pps(yimp_logn, samp_logn$pik)
      res <- cbind(distnum=c(2,2,3,3), rbind(postmi_unif, postmi_logn))
    }
    return(res)
  }
  stopCluster(cl)
  b <- Sys.time()
  
  # Get true means from population
  popN <- N
  true <- read_csv("./data_pop/pop_sumstats.csv",col_types = cols()) %>% filter(N==popN & rho==RHO)
  post$truemean[post$distnum==1] <- true$meany_norm
  post$truemean[post$distnum==2] <- true$meany_unif
  post$truemean[post$distnum==3] <- true$meany_logn
  
  # Additional info from replicate
  post <- post %>% mutate(N=N, n=SAMPSIZE, rho=RHO, missmech=MISSMECH, k=K, D=D, sampling=TYPE)
  return(post)
}

############
### TYPE = PPS
TYPE <- "pps"
############
set.seed(97423)
SEEDS <- round(runif(1080, 1, 99999))
j <- 1
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (n in c(100,200,400))
  {
    # 10 correlations
    for (r in 0:9)
    {
      # 6 missingness mechanisms
      for (miss in c("mcar25","mcar50","mar25pos","mar25neg","mar50pos","mar50neg"))
      {
        # 3 imputation methods
        for (k in c(1,5,0))
        {
          print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10, "k =",k, "missing mech =",miss))
          set.seed(SEEDS[j])
          postmi <- impute(N, n, r/10, miss,  k, 20, TYPE)
          # save info from replicates
          # write_csv(postmi, file=paste("./results/nn_rhd/results_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k",k,".csv",sep=""))
          # get summary info across replicates and save
          postmi$se <- sqrt(postmi$var)
          stats <- postmi %>%
            simsum(estvarname="mean", se="se", true="truemean", by=c("varmethod","distnum","N","n","rho","missmech","k","D","sampling")) %>%
            summary() %>%
            get_data()
          write_csv(stats, file=paste("./results/nn_rhd/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k",k,".csv",sep=""))
          j <- j+1
        }
      }
    }
  }
}



############
### TYPE = SRS
TYPE <- "srs"
############
set.seed(1002231)
SEEDS <- round(runif(1080, 1, 99999))
j <- 1
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (n in c(100,200,400))
  {
    # 10 correlations
    for (r in 0:9)
    {
      # 6 missingness mechanisms
      for (miss in c("mcar25","mcar50","mar25pos","mar25neg","mar50pos","mar50neg"))
      {
        # 3 imputation methods
        for (k in c(1,5,0))
        {
          print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10, "k =",k, "missing mech =",miss))
          set.seed(SEEDS[j])
          postmi <- impute(N, n, r/10, miss,  k, 20, TYPE)
          # save info from replicates
          # write_csv(postmi, file=paste("./results/nn_rhd/results_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k",k,".csv",sep=""))
          # get summary info across replicates and save
          postmi$se <- sqrt(postmi$var)
          stats <- postmi %>%
            simsum(estvarname="mean", se="se", true="truemean", by=c("varmethod","distnum","N","n","rho","missmech","k","D","sampling")) %>%
            summary() %>%
            get_data()
          write_csv(stats, file=paste("./results/nn_rhd/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k",k,".csv",sep=""))
          j <- j+1
        }
      }
    }
  }
}

