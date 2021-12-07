####################
# NN MI Simulation
# Impute with weighted RHD + ABB
# Author: Rebecca Andridge
# Last Modified: 5/14/2021
####################
rm(list=ls())

require(doParallel)
require(foreach)
require(doRNG)
require(tidyverse)
require(survey)
require(rsimsum)

# # For testing
# N <- 1000
# SAMPSIZE <- 100
# RHO <- 0.5
# MISSMECH <- "mcar50"
# D <- 20  # number of imputations

#### ONLY USED FOR PPS SAMPLING (not different from RHD when sampling is SRS)

impute <- function(N, SAMPSIZE, RHO, MISSMECH, D)
{
  # Load data
  rho <- RHO*10
  load(paste(paste("./data_samp/samples_pps_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
  load(paste(paste("./data_samp/samples_pps_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",rho,".RData",sep="")))
  # Sample weights
  samples_unif$w <- 1/samples_unif$pik
  samples_logn$w <- 1/samples_logn$pik
  # Induce missingness for specified missingness mechanism
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
    source("./functions/wrhd_impute.R")
    source("./functions/miStats_pps.R")
    # Get samples
    samp_unif <- as_tibble(subset(samples_unif, sample==i))
    samp_logn <- as_tibble(subset(samples_logn, sample==i))
    # Multiply impute using weighted RHD
    yimp_unif <- wrhd_impute(samp_unif$y, samp_unif$w, D)
    yimp_logn <- wrhd_impute(samp_logn$y, samp_logn$w, D)
    # Apply Rubin's rules
    postmi_unif <- miStats_pps(yimp_unif, samp_unif$pik)
    postmi_logn <- miStats_pps(yimp_logn, samp_logn$pik)
    # Return statistics as matrix
    return(cbind(distnum=c(2,2,3,3), rbind(postmi_unif, postmi_logn)))
  }
  stopCluster(cl)
  b <- Sys.time()
  
  # Get true means from population
  popN <- N
  true <- read_csv("./data_pop/pop_sumstats.csv",col_types = cols()) %>% filter(N==popN & rho==RHO)
  post$truemean[post$distnum==2] <- true$meany_unif
  post$truemean[post$distnum==3] <- true$meany_logn
  
  # Additional info from replicate
  post <- post %>% mutate(N=N, n=SAMPSIZE, rho=RHO, missmech=MISSMECH, k=10, D=D, sampling="pps") #k=10 --> wRHD
  rm(samples_logn, samples_unif)
  return(post)
}

############
### TYPE = PPS
############
set.seed(87243)
SEEDS <- round(runif(1000, 1, 99999))
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
        print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10, "k = 10 missing mech =",miss))
        set.seed(SEEDS[j])
        a <- Sys.time()
        postmi <- impute(N, n, r/10, miss, 20)
        b <- Sys.time()
        print(b-a)
        # save info from replicates
        write_csv(postmi, file=paste("./results/wrhd/results_pps_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k10.csv",sep=""))
        # get summary info across replicates and save
        postmi$se <- sqrt(postmi$var)
        stats <- postmi %>%
          simsum(estvarname="mean", se="se", true="truemean", by=c("varmethod","distnum","N","n","rho","missmech","k","D","sampling")) %>% 
          summary() %>%
          get_data()
        write_csv(stats, file=paste("./results/wrhd/sumstats_pps_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_k10.csv",sep=""))
        j <- j+1
      }
    }
  }
}

