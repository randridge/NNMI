####################
# NN MI Simulation
# Generate MAR missing data in samples
# Author: Rebecca Andridge
# Last Modified: 5/3/2021
####################
rm(list=ls())

require(tidyverse)

########################################
# Induce MAR in SRS samples
marSRS <- function()
{
  # Load samples
  load(paste("./data_samp_bd/samples_srs_norm_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  load(paste("./data_samp_bd/samples_srs_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  load(paste("./data_samp_bd/samples_srs_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  
  NREPS <- max(samples_unif$sample)
  
  # Missingness indicators - MAR
  # Normal - 25% missing
  p_pos <- 1 - 1/(1+exp(-1.14 + log(1.5)*samples_norm$x))   # larger units missing
  p_neg <- 1 - 1/(1+exp(-1.14 - log(1.5)*samples_norm$x))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_norm <- samples_norm %>% mutate(mar25pos=(m_pos==1),mar25neg=(m_neg==1))
  # Normal - 50% missing
  p_pos <- 1 - 1/(1+exp(0 + log(1.5)*samples_norm$x))   # larger units missing
  p_neg <- 1 - 1/(1+exp(0 - log(1.5)*samples_norm$x))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_norm <- samples_norm %>% mutate(mar50pos=(m_pos==1),mar50neg=(m_neg==1))
  # Uniform - 25% missing
  p_pos <- 1 - 1/(1+exp(-1.14 + log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # larger units missing
  p_neg <- 1 - 1/(1+exp(-1.14 - log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_unif <- samples_unif %>% mutate(mar25pos=(m_pos==1),mar25neg=(m_neg==1))
  # Uniform - 50% missing
  p_pos <- 1 - 1/(1+exp(0 + log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # larger units missing
  p_neg <- 1 - 1/(1+exp(0 - log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_unif <- samples_unif %>% mutate(mar50pos=(m_pos==1),mar50neg=(m_neg==1))
  # Lognormal - 25% missing
  p_pos <- 1 - 1/(1+exp(-1.14 + log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # larger units missing
  p_neg <- 1 - 1/(1+exp(-1.14 - log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_logn <- samples_logn %>% mutate(mar25pos=(m_pos==1),mar25neg=(m_neg==1))
  # Lognormal - 50% missing
  p_pos <- 1 - 1/(1+exp(0 + log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # larger units missing
  p_neg <- 1 - 1/(1+exp(0 - log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_logn <- samples_logn %>% mutate(mar50pos=(m_pos==1),mar50neg=(m_neg==1))

  # Save
  rm(p_pos, p_neg, m_pos, m_neg)
  save(samples_norm, file=paste("./data_samp/samples_srs_norm_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_unif, file=paste("./data_samp/samples_srs_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_logn, file=paste("./data_samp/samples_srs_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  
  # Percent missing
  pmis <- cbind(distnum=1:3, rbind(samples_norm %>% select(starts_with("m")) %>% colMeans,
                                   samples_unif %>% select(starts_with("m")) %>% colMeans,
                                   samples_logn %>% select(starts_with("m")) %>% colMeans))
  pmis <- as_tibble(cbind(N=N, n=SAMPSIZE, rho=j/10, pmis))
  return(pmis)
}

pmiss <- vector("list", length=2*3*10)
i <- 0
# 2 population sizes
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (SAMPSIZE in c(100,200,400))
  {
    # 10 correlations
    for (j in 0:9)
    {
      i <- i + 1
      print(paste("N=",format(N,scientific=F)," n=",SAMPSIZE," rho=",j/10,sep=""))
      pmiss[[i]] <- marSRS()
    }
  }
}
srs_pmiss <- bind_rows(pmiss)
hist(srs_pmiss$mar25pos)
hist(srs_pmiss$mar25neg)
hist(srs_pmiss$mar50pos)
hist(srs_pmiss$mar50neg)
write_csv(srs_pmiss,"./data_samp/srs_pmiss.csv")

########################################
# Induce MAR in PPS samples
marPPS <- function()
{
  # Load samples
  load(paste("./data_samp_bd/samples_pps_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  load(paste("./data_samp_bd/samples_pps_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  
  NREPS <- max(samples_unif$sample)
  
  # Missingness indicators - MAR
  # Uniform - 25% missing
  p_pos <- 1 - 1/(1+exp(-1.32 + log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # larger units missing
  p_neg <- 1 - 1/(1+exp(-0.94 - log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_unif <- samples_unif %>% mutate(mar25pos=(m_pos==1),mar25neg=(m_neg==1))
  # Uniform - 50% missing
  p_pos <- 1 - 1/(1+exp(-0.19 + log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # larger units missing
  p_neg <- 1 - 1/(1+exp(0.19 - log(1.5)/sqrt((.9-.1)^2/12)*(samples_unif$x-0.5)))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_unif <- samples_unif %>% mutate(mar50pos=(m_pos==1),mar50neg=(m_neg==1))
  # Lognormal
  if (N==1000 & SAMPSIZE==100)
  {
    ALPHA <- c(-1.59, -0.81, -0.36, 0.36)
  } else if (N==1000 & SAMPSIZE==200)
  {
    ALPHA <- c(-1.49, -0.85, -0.31, 0.31)
  } else if (N==1000 & SAMPSIZE==400)
  {
    ALPHA <- c(-1.35, -0.94, -0.20, 0.20)
  } else 
  {
    ALPHA <- c(-1.69, -0.76, -0.43, 0.43)
  }
  # Lognormal - 25% missing
  p_pos <- 1 - 1/(1+exp(ALPHA[1] + log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # larger units missing
  p_neg <- 1 - 1/(1+exp(ALPHA[2] - log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_logn <- samples_logn %>% mutate(mar25pos=(m_pos==1),mar25neg=(m_neg==1))
  # Lognormal - 50% missing
  p_pos <- 1 - 1/(1+exp(ALPHA[3] + log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # larger units missing
  p_neg <- 1 - 1/(1+exp(ALPHA[4] - log(1.5)/sqrt(exp(1)*(exp(1)-1))*(samples_logn$x-exp(0.5))))   # smaller units missing
  m_pos <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_pos)
  m_neg <- rbinom(n=SAMPSIZE*NREPS, size=1, prob=p_neg)
  samples_logn <- samples_logn %>% mutate(mar50pos=(m_pos==1),mar50neg=(m_neg==1))
  
  # Save
  rm(p_pos, p_neg, m_pos, m_neg)
  save(samples_unif, file=paste("./data_samp/samples_pps_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_logn, file=paste("./data_samp/samples_pps_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  
  # Percent missing
  pmis <- cbind(distnum=2:3, rbind(samples_unif %>% select(starts_with("m")) %>% colMeans,
                                   samples_logn %>% select(starts_with("m")) %>% colMeans))
  pmis <- as_tibble(cbind(N=N, n=SAMPSIZE, rho=j/10, pmis))
  return(pmis)
}

pmiss <- vector("list", length=2*3*10)
i <- 0
# 2 population sizes
for (N in c(1000,100000))
{
  # 3 sample sizes
  for (SAMPSIZE in c(100,200,400))
  {
    # 10 correlations
    for (j in 0:9)
    {
      i <- i + 1
      print(paste("N=",format(N,scientific=F)," n=",SAMPSIZE," rho=",j/10,sep=""))
      pmiss[[i]] <- marPPS()
    }
  }
}
pps_pmiss <- bind_rows(pmiss)
hist(pps_pmiss$mar25pos)
hist(pps_pmiss$mar25neg)
hist(pps_pmiss$mar50pos)
hist(pps_pmiss$mar50neg)
write_csv(pps_pmiss,"./data_samp/pps_pmiss.csv")
