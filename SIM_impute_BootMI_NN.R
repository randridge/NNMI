####################
# NN MI Simulation
# Impute with BootMI (nImp=2, nBoot=500)
# Author: Rebecca Andridge
# Last Modified: 5/17/2021
####################
rm(list=ls())

require(doParallel)
require(foreach)
require(doRNG)
require(bootImpute)
require(tidyverse)
require(rsimsum)

# # # For testing
# N <- 1000
# SAMPSIZE <- 100
# RHO <- 0.5
# MISSMECH <- "mcar50"  # mcar25, mcar50, mar25pos, mar25neg, mar50pos mar50neg
# TYPE <- "pps"  # "srs" or "pps"
impute <- function(N, SAMPSIZE, RHO, MISSMECH, TYPE)
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
  # Parallel processing
  NREPS <- max(samples_unif$sample)
  cl <- makeCluster(8, setup_strategy="sequential")
  registerDoParallel(cl)
  # Looping through simulated data sets
  post <- foreach(i=1:NREPS, .combine=rbind, .packages=c("bootImpute","tidyverse")) %dorng% {
    # Functions
    source("./functions/nn_impute_df.R")
    source("./functions/bootImputeAnalyse_fpc.R")
    # Get samples
    if (TYPE=="srs") samp_norm <- as_tibble(subset(samples_norm, sample==i))
    samp_unif <- as_tibble(subset(samples_unif, sample==i))
    samp_logn <- as_tibble(subset(samples_logn, sample==i))
    # Multiply impute using bootstrapped MI with D=2 and NBOOT=500
    nBoot <- 500 # 500 bootstraps
    nImp <- 2    # 2 imputations per bootstrap sample
    imps_norm <- imps_unif <- imps_logn <- vector("list", nBoot*nImp)
    count <- 1
    for (b in 1:nBoot)
    {
      # check to ensure bootstrapped sample included at least some respondents, and if not then resample
      noResp <- TRUE
      while (noResp)
      {
        bsIndices <- sample(1:SAMPSIZE, replace=TRUE)
        if (TYPE=="srs") {
          noResp <- (sum(!is.na(samp_norm$y[bsIndices]))==0) | (sum(!is.na(samp_unif$y[bsIndices]))==0) | (sum(!is.na(samp_logn$y[bsIndices]))==0)
        } else {
          noResp <- (sum(!is.na(samp_unif$y[bsIndices]))==0) | (sum(!is.na(samp_logn$y[bsIndices]))==0)
        }
      }
      for (m in 1:nImp)
      {
        if (TYPE=="srs") imps_norm[[count]] <- nn_impute_df(samp_norm[bsIndices,c("x","y","w")], K=1)
        imps_unif[[count]] <- nn_impute_df(samp_unif[bsIndices,c("x","y","w")], K=1)
        imps_logn[[count]] <- nn_impute_df(samp_logn[bsIndices,c("x","y","w")], K=1)
        count <- count + 1
      }
    }
    # Add attributes needed for analyze function
    attributes(imps_norm) <- attributes(imps_unif) <- attributes(imps_logn) <- list(nBoot=nBoot, nImp=nImp)
    # Combining rules per Von Hippel (2018)
    # Function to estimate the mean of Y from the dataset
    meanfunc <- function(DATA) sum(DATA$y*DATA$w)/sum(DATA$w)
    if (TYPE=="srs") postmi_norm <- bootImputeAnalyse(imps_norm, meanfunc, quiet=TRUE)
    postmi_unif <- bootImputeAnalyse(imps_unif, meanfunc, quiet=TRUE)
    postmi_logn <- bootImputeAnalyse(imps_logn, meanfunc, quiet=TRUE)
    # Combining rules per Von Hippel (2018) - with incorporation of FPC based on # respondents
    if (TYPE=="srs") fpc_norm <- 1 - sum(!is.na(samp_norm$y))/N  # FPC based on number of respondents
    fpc_unif <- 1 - sum(!is.na(samp_unif$y))/N
    fpc_logn <- 1 - sum(!is.na(samp_logn$y))/N
    if (TYPE=="srs") postmi_norm_fpc <- bootImputeAnalyse_fpc(fpc_norm, imps_norm, meanfunc, quiet=TRUE)
    postmi_unif_fpc <- bootImputeAnalyse_fpc(fpc_unif, imps_unif, meanfunc, quiet=TRUE)
    postmi_logn_fpc <- bootImputeAnalyse_fpc(fpc_unif, imps_logn, meanfunc, quiet=TRUE)
    # Return MI statistics as matrix
    if (TYPE=="srs"){
      res_nofpc <- cbind(distnum=1:3, k=c(99,99,99),
                         rbind(unlist(postmi_norm), unlist(postmi_unif), unlist(postmi_logn)),   # ests, var, lb, ub, df
                         var_ml=rep(NA,3), var_bd=rep(NA,3))
      res_fpc   <- cbind(distnum=c(1,2,3), k=c(88,88,88),
                         rbind(unlist(postmi_norm_fpc), unlist(postmi_unif_fpc), unlist(postmi_logn_fpc)))  # var, lb, ub, df, var_ml, var_bd
    } else {
      res_nofpc <- cbind(distnum=2:3, k=c(99,99),
                         rbind(unlist(postmi_unif), unlist(postmi_logn)),   # ests, var, lb, ub, df
                         var_ml=rep(NA,2), var_bd=rep(NA,2))
      res_fpc   <- cbind(distnum=c(2,3), k=c(88,88),
                         rbind(unlist(postmi_unif_fpc), unlist(postmi_logn_fpc)))  # var, lb, ub, df, var_ml, var_bd
      
    }
    return(rbind(res_nofpc, res_fpc))
  }
  stopCluster(cl)
  # result is matrix
  post <- as.data.frame(post)
  row.names(post) <- NULL
  names(post) <- c("distnum","k","mean","var","lb","ub","df","var_ml","var_bd")
  # Get true means from population
  load(paste("./data_pop/pop_N",format(N,scientific=F),"_rho",rho,".RData",sep=""))
  post$truemean[post$distnum==1] <- pop$stats$meany_norm
  post$truemean[post$distnum==2] <- pop$stats$meany_unif
  post$truemean[post$distnum==3] <- pop$stats$meany_logn
  # Additional info from replicate
  post <- post %>% mutate(N=N, n=SAMPSIZE, rho=RHO, missmech=MISSMECH, sampling=TYPE)
  rm(pop)
  return(post)
}

############
### TYPE = PPS
TYPE <- "pps"
############
set.seed(8520019)
SEEDS <- runif(1000, 1, 99999)
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
        print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10, "k = BootMI missing mech =",miss))
        set.seed(SEEDS[j])
        a <- Sys.time()
        postmi <- impute(N, n, r/10, miss, TYPE)
        b <- Sys.time()
        print(b-a)
        # save info from replicates
        write_csv(postmi, file=paste("./results/bootmi/results_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_bootmi.csv",sep=""))
        # get summary info across replicates and save
        postmi$varmethod <- 9
        postmi$se <- sqrt(postmi$var)
        stats <- postmi %>%
          simsum(estvarname="mean", se="se", true="truemean", by=c("varmethod","distnum","N","n","rho","missmech","k","sampling")) %>% 
          summary() %>%
          get_data()
        write_csv(stats, file=paste("./results/bootmi/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_bootmi.csv",sep=""))
        j <- j+1
      }
    }
  }
}


############
### TYPE = SRS
TYPE <- "srs"
############
set.seed(178432)
SEEDS <- runif(1000, 1, 99999)
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
        postmi <- impute(N, n, r/10, miss, TYPE)
        # save info from replicates
        # write_csv(postmi, file=paste("./results/bootmi/results_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_bootmi.csv",sep=""))
        # get summary info across replicates and save
        postmi$varmethod <- 9
        postmi$se <- sqrt(postmi$var)
        stats <- postmi %>%
          simsum(estvarname="mean", se="se", true="truemean", by=c("distnum","N","n","rho","missmech","k","sampling")) %>% 
          summary() %>%
          get_data()
        write_csv(stats, file=paste("./results/bootmi/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,"_",miss,"_bootmi.csv",sep=""))
        j <- j+1
      }
    }
  }
}
