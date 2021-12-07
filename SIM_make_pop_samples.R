####################
# NN MI Simulation
# Generate finite population and draw samples
# Author: Rebecca Andridge
# Last Modified: 5/6/2021
####################
rm(list=ls())

require(sampling)
require(tidyverse)

########################################
# MAKE POPULATIONS
########################################
# Function to generate finite population
makepop <- function(RHO, N, SEED)
{
  set.seed(SEED)
  ## Normal -- E(X)=0, V(X)=1; E(Y)=1, V(Y)=1
  x_norm <- rnorm(N)                                 # X ~ N(0,1)
  y_norm <- 1 + RHO*x_norm + rnorm(N)*sqrt(1-RHO^2)  # Y|X ~ N(1+RHO*X, 1-RHO^2) --> Corr(X,Y)=RHO
  ## Uniform -- E(X)=0.5, V(X)=(.9-.1)^2/12; E(Y)=0.5, V(Y)=1/12
  x_unif <- 0.1+(0.9-0.1)*pnorm(x_norm)              # X ~ Uniform(0.1,0.9)
  y_unif <- pnorm(y_norm-1)                          # Y ~ Uniform(0,1)
  ## Lognormal -- E(X)=exp(0.5), V(X)=exp(1)*(exp(1)-1); E(Y)=exp(0.5), V(Y)=exp(1)*(exp(1)-1)
  x_logn <- exp(x_norm)
  y_logn <- exp(y_norm-1)
  # Make dataframe
  df <- data.frame(N=N, rho=RHO, id=1:N, 
                   x_norm=x_norm, y_norm=y_norm,
                   x_unif=x_unif, y_unif=y_unif,
                   x_logn=x_logn, y_logn=y_logn)
  # Population summary stats
  stats <- list(N=N, rho=RHO,
                meany_norm=mean(y_norm), meany_unif=mean(y_unif), meany_logn=mean(y_logn), 
                meanx_norm=mean(x_norm), meanx_unif=mean(x_unif), meanx_logn=mean(x_logn), 
                pcor_norm=cor(x_norm,y_norm), pcor_unif=cor(x_unif,y_unif), pcor_logn=cor(x_logn,y_logn), 
                scor=cor(x_norm,y_norm,method="spearman"))
  return(list(dat=df, stats=stats))
}

# Make populations
# using same seed for all means the X (size variable) is the same for all RHO values
sumstats <- vector()
for (i in 0:9)
{
  print(i)
  # Small population (N=1,000)
  pop <- makepop(i/10, 1000, 531217)
  sumstats <- rbind(sumstats,unlist(pop$stats))
  save(pop, file=paste("./data_pop/pop_N1000_rho",i,".RData",sep=""))
  # Large population (N=100,000)
  pop <- makepop(i/10, 100000, 531217)
  sumstats <- rbind(sumstats,unlist(pop$stats))
  save(pop, file=paste("./data_pop/pop_N100000_rho",i,".RData",sep=""))
}
# save population summary stats
write_csv(as.data.frame(sumstats), file="./data_pop/pop_sumstats.csv")
# # check XY correlations in finite populations
# long <- as_tibble(sumstats) %>% pivot_longer(starts_with("pcor"), names_to="dist", values_to="pcor")
# ggplot(long, aes(x=rho, y=pcor, group=dist, shape=dist, col=dist)) + theme_bw() + geom_line() + geom_point() + facet_grid(.~N)
# rm(list=ls())
# as_tibble(sumstats) %>% filter(N==1000)
# as_tibble(sumstats) %>% filter(N==100000)

########################################
# Draw SRS samples
NREPS = 5000
########################################
samp_srs <- function(N, SAMPSIZE, RHO, NREPS, SEED)
{
  # Load finite population dataset
  j <- RHO*10
  load(paste("./data_pop/pop_N",format(N,scientific=F),"_rho",j,".RData",sep=""))
  
  # Take repeated SRS from finite population
  # List to hold sampled IDs
  mat_ids <- vector("list", length=NREPS)
  set.seed(SEED)
  for (i in 1:NREPS){
    mat_ids[[i]] <- sample(1:N, SAMPSIZE, replace=FALSE)
  }
  # create dataframes stacking sampled IDs
  samples <- as_tibble(cbind(sample=sort(rep(1:NREPS,SAMPSIZE)), n=SAMPSIZE, pop$dat[unlist(mat_ids),]))
  
  # Missingness indicators - MCAR
  rand <- runif(SAMPSIZE*NREPS)
  samples <- samples %>% mutate(mcar25=(rand<0.25), mcar50=(rand<0.50))
  
  # split into separate DF for each distribution
  samples_norm <- samples %>% select(-(contains("unif")|contains("logn"))) %>% rename(x=x_norm, y=y_norm)
  samples_unif <- samples %>% select(-(contains("norm")|contains("logn"))) %>% rename(x=x_unif, y=y_unif)
  samples_logn <- samples %>% select(-(contains("norm")|contains("unif"))) %>% rename(x=x_logn, y=y_logn)
 
  # Save
  rm(rand, samples, pop, mat_ids)
  save(samples_norm, file=paste("./data_samp_bd/samples_srs_norm_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_unif, file=paste("./data_samp_bd/samples_srs_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_logn, file=paste("./data_samp_bd/samples_srs_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
}
# Draw samples
set.seed(971726)
SEEDS <- round(runif(60, 1, 999999))
j <- 1
for (i in 0:9)
{
  print(i)
  samp_srs(  1000, 100, i/10, NREPS, SEEDS[j]); j <- j + 1
  samp_srs(  1000, 200, i/10, NREPS, SEEDS[j]); j <- j + 1
  samp_srs(  1000, 400, i/10, NREPS, SEEDS[j]); j <- j + 1
  samp_srs(100000, 100, i/10, NREPS, SEEDS[j]); j <- j + 1
  samp_srs(100000, 200, i/10, NREPS, SEEDS[j]); j <- j + 1
  samp_srs(100000, 400, i/10, NREPS, SEEDS[j]); j <- j + 1
}
rm(list=ls())


########################################
# Draw PPS samples
NREPS = 5000
########################################
samp_pps <- function(N, SAMPSIZE, RHO, NREPS, SEED)
{
  # Load finite population dataset
  j <- RHO*10
  load(paste("./data_pop/pop_N",format(N,scientific=F),"_rho",j,".RData",sep=""))
  
  # Take repeated PPS from finite population (Uniform and Lognormal only)
  # PPS sampling proportional to X
  # Inclusion probabilities
  pik_unif <- inclusionprobabilities(pop$dat$x_unif, SAMPSIZE)
  pik_logn <- inclusionprobabilities(pop$dat$x_logn, SAMPSIZE)
  # table(pik_logn==1)  # N=1,000: 120 are certainty for n=400, 18 for n=200, 3 for n=100
  
  # lists to hold sampled IDs
  mat_ids_unif <- mat_ids_logn <- vector("list", length=NREPS)
  # sample via systematic PPS
  set.seed(SEED)
  for (i in 1:NREPS){
    # if (i %% 500 == 0) print(paste("i=",i,sep=""))
    s_unif <- UPrandomsystematic(pik_unif)
    s_logn <- UPrandomsystematic(pik_logn)
    mat_ids_unif[[i]] <- c(1:N)[s_unif==1]
    mat_ids_logn[[i]] <- c(1:N)[s_logn==1]
  }
  # create dataframe stacking sampled IDs
  samples_unif <- as_tibble(cbind(sample=sort(rep(1:NREPS,SAMPSIZE)), N=N, n=SAMPSIZE, rho=RHO, id=unlist(mat_ids_unif)))
  samples_logn <- as_tibble(cbind(sample=sort(rep(1:NREPS,SAMPSIZE)), N=N, n=SAMPSIZE, rho=RHO, id=unlist(mat_ids_logn)))
  samples_unif <- samples_unif %>% mutate(pik=pik_unif[samples_unif$id], x=pop$dat$x_unif[samples_unif$id], y=pop$dat$y_unif[samples_unif$id])
  samples_logn <- samples_logn %>% mutate(pik=pik_logn[samples_logn$id], x=pop$dat$x_logn[samples_logn$id], y=pop$dat$y_logn[samples_logn$id])

  # Missingness indicators - MCAR
  rand <- runif(SAMPSIZE*NREPS)
  samples_unif <- samples_unif %>% mutate(mcar25=(rand<0.25), mcar50=(rand<0.50))
  samples_logn <- samples_logn %>% mutate(mcar25=(rand<0.25), mcar50=(rand<0.50))
  
  # Save
  rm(rand, mat_ids_unif, mat_ids_logn, pik_unif, pik_logn, s_unif, s_logn, pop)
  save(samples_unif, file=paste("./data_samp_bd/samples_pps_unif_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
  save(samples_logn, file=paste("./data_samp_bd/samples_pps_logn_N",format(N,scientific=F),"_n",SAMPSIZE,"_rho",j,".RData",sep=""))
}

# using same seed for all means would mean the same IDs chosen for all RHOs
set.seed(8810028)
SEEDS <- round(runif(60, 1, 999999))
j <- 1
for (i in 0:9)
{
  print(paste("RHO =",i))
  samp_pps(  1000, 100, i/10, NREPS, SEEDS[j]); print("1000/100")
  j <- j + 1
  samp_pps(  1000, 200, i/10, NREPS, SEEDS[j]); print("1000/200")
  j <- j + 1
  samp_pps(  1000, 400, i/10, NREPS, SEEDS[j]); print("1000/400")
  j <- j + 1
  samp_pps(100000, 100, i/10, NREPS, SEEDS[j]); print("100000/100")
  j <- j + 1
  samp_pps(100000, 200, i/10, NREPS, SEEDS[j]); print("100000/200")
  j <- j + 1
  samp_pps(100000, 400, i/10, NREPS, SEEDS[j]); print("100000/400")
  j <- j + 1
}
