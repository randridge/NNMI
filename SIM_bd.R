####################
# NN MI Simulation
# Before deletion design-based estimates
# Author: Rebecca Andridge
# Last Modified: 5/1/2021
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
# SAMPSIZE <- 400
# RHO <- 0.5
# TYPE <- "pps"  # "srs" or "pps"

bd <- function(N, SAMPSIZE, RHO, TYPE)
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
  NREPS <- max(samples_unif$sample)
  a <- Sys.time()
  # Parallel processing
  cl <- makeCluster(8, setup_strategy="sequential")
  registerDoParallel(cl)
  # Looping through simulated data sets
  post <- foreach(i=1:NREPS, .combine=rbind, .packages=c("tidyverse","survey")) %dorng% {
    # Get samples
    if (TYPE=="srs") samp_norm <- as_tibble(subset(samples_norm, sample==i))
    samp_unif <- as_tibble(subset(samples_unif, sample==i))
    samp_logn <- as_tibble(subset(samples_logn, sample==i))
    # Sample design
    if (TYPE=="srs"){
      des_norm <- svydesign(id=~1, weights=~w, fpc=~N, data=samp_norm)
      des_unif <- svydesign(id=~1, weights=~w, fpc=~N, data=samp_unif)
      des_logn <- svydesign(id=~1, weights=~w, fpc=~N, data=samp_logn)
    } else if (TYPE=="pps") {
      des_unif_wr <- svydesign(id=~1, probs=~pik, data=samp_unif) # with-replacement approximation
      des_logn_wr <- svydesign(id=~1, probs=~pik, data=samp_logn)
      des_unif_br <- svydesign(id=~1, fpc=~pik, data=samp_unif, pps="brewer") # Brewer PPSWOR approximation
      des_logn_br <- svydesign(id=~1, fpc=~pik, data=samp_logn, pps="brewer")
      des_unif_ov <- svydesign(id=~1, fpc=~pik, data=samp_unif, pps="overton") # Overton PPSWOR approximation
      des_logn_ov <- svydesign(id=~1, fpc=~pik, data=samp_logn, pps="overton")
      des_unif_hr <- svydesign(id=~1, fpc=~pik, data=samp_unif, pps=HR()) # Hartley-Rao PPSWOR approximation
      des_logn_hr <- svydesign(id=~1, fpc=~pik, data=samp_logn, pps=HR())
    }
    # Sample estimates
    # Return statistics as matrix
    if (TYPE=="srs")
    {
      norm <- svymean(~y, design=des_norm)
      unif <- svymean(~y, design=des_unif)
      logn <- svymean(~y, design=des_logn)
      norm_ci <- confint(norm)
      unif_ci <- confint(unif)
      logn_ci <- confint(logn)
      stats <- as.data.frame(cbind(distnum=1:3,
                                   mean=c(coef(norm),coef(unif),coef(logn)),
                                   var=c(vcov(norm),vcov(unif),vcov(logn)),
                                   lb=c(norm_ci[1],unif_ci[1],logn_ci[1]),
                                   ub=c(norm_ci[2],unif_ci[2],logn_ci[2])))
    } else if (TYPE=="pps") {
      unif_wr <- svymean(~y, design=des_unif_wr) # with-replacement approximation
      logn_wr <- svymean(~y, design=des_logn_wr)
      unif_br <- svymean(~y, design=des_unif_br) # Brewer PPSWOR approximation
      logn_br <- svymean(~y, design=des_logn_br)
      unif_ov <- svymean(~y, design=des_unif_ov) # Overton PPSWOR approximation
      logn_ov <- svymean(~y, design=des_logn_ov)
      unif_hr <- svymean(~y, design=des_unif_hr) # Hartley-Rao PPSWOR approximation
      logn_hr <- svymean(~y, design=des_logn_hr)
      unif_wr_ci <- confint(unif_wr)
      logn_wr_ci <- confint(logn_wr)
      unif_br_ci <- confint(unif_br)
      logn_br_ci <- confint(logn_br)
      unif_ov_ci <- confint(unif_ov)
      logn_ov_ci <- confint(logn_ov)
      unif_hr_ci <- confint(unif_hr)
      logn_hr_ci <- confint(logn_hr)
      stats_unif <- as.data.frame(cbind(distnum=2, varmethod=1:4,
                                        mean=c(coef(unif_wr),coef(unif_br),coef(unif_ov),coef(unif_hr)),
                                        var=c(vcov(unif_wr),vcov(unif_br),vcov(unif_ov),vcov(unif_hr)),
                                        lb=c(unif_wr_ci[1],unif_br_ci[1],unif_ov_ci[1],logn_hr_ci[1]),
                                        ub=c(unif_wr_ci[2],unif_br_ci[2],unif_ov_ci[2],logn_hr_ci[2])))
      stats_logn <- as.data.frame(cbind(distnum=3, varmethod=1:4,
                                        mean=c(coef(logn_wr),coef(logn_br),coef(logn_ov),coef(logn_hr)),
                                        var=c(vcov(logn_wr),vcov(logn_br),vcov(logn_ov),vcov(logn_hr)),
                                        lb=c(logn_wr_ci[1],logn_br_ci[1],logn_ov_ci[1],logn_hr_ci[1]),
                                        ub=c(logn_wr_ci[2],logn_br_ci[2],logn_ov_ci[2],logn_hr_ci[2])))
      stats <- rbind(stats_unif, stats_logn)
    }
    row.names(stats) <- NULL
    return(stats)
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
  post <- post %>% mutate(N=N, n=SAMPSIZE, rho=RHO, k=99999, sampling=TYPE) #k=99999 --> BD
  
  return(post)
}

############
### TYPE = PPS
TYPE <- "pps"
############
# 2 population sizes
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (n in c(100,200,400))
  {
    # 10 correlations
    for (r in 0:9)
    {
      print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10))
      result <- bd(N, n, r/10, TYPE)
      # get summary info across replicates and save
      result$se <- sqrt(result$var)
      stats <- result %>%
        simsum(estvarname="mean", se="se", true="truemean", by=c("varmethod","distnum","N","n","rho","k","sampling")) %>% 
        summary() %>%
        get_data()
      write_csv(stats, file=paste("./results/bd/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,".csv",sep=""))
    }
  }
}

############
### TYPE = PPS
TYPE <- "srs"
############
# 2 population sizes
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (n in c(100,200,400))
  {
    # 10 correlations
    for (r in 0:9)
    {
      print(paste("N =",format(N,scientific=N), "n =",n, "rho =",r/10))
      result <- bd(N, n, r/10, TYPE)
      # get summary info across replicates and save
      result$se <- sqrt(result$var)
      stats <- result %>%
        simsum(estvarname="mean", se="se", true="truemean", by=c("distnum","N","n","rho","k","sampling")) %>% 
        summary() %>%
        get_data()
      write_csv(stats, file=paste("./results/bd/sumstats_",TYPE,"_N",format(N,scientific=F),"_n",n,"_rho",r,".csv",sep=""))
    }
  }
}
