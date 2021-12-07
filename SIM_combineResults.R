####################
# NN MI Simulation
# Author: Rebecca Andridge
# Last Modified: 8/5/2021
####################

library(tidyverse)

# load results - NN1, NN5, RHD --> k=0 (RHD), 1,5 (NN)
files <- dir("./results/nn_rhd", pattern = "sumstats") # get file names
nn_rhd <- files %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(file.path("./results/nn_rhd", .), col_types=cols())) %>% 
  reduce(rbind)

# load results - wRHD --> k=10
files <- dir("./results/wrhd", pattern = "sumstats") # get file names
wrhd <- files %>%
  map(~ read_csv(file.path("./results/wrhd", .), col_types=cols())) %>% 
  reduce(rbind)

# load results - BootMI --> k=88 (w/FPC), 99 (no FPC)
files_pps <- dir("./results/bootmi", pattern = "sumstats_pps") # get file names
bmi_pps <- files_pps %>%
  map(~ read_csv(file.path("./results/bootmi", .), col_types=cols())) %>% 
  reduce(rbind)
files_srs <- dir("./results/bootmi", pattern = "sumstats_srs") # get file names
bmi_srs <- files_srs %>%
  map(~ read_csv(file.path("./results/bootmi", .), col_types=cols())) %>% 
  reduce(rbind)
bmi_srs$varmethod <- 9
bmi <- bind_rows(bmi_pps, bmi_srs)
rm(bmi_pps, bmi_srs)

ftable(nn_rhd$sampling, nn_rhd$N, nn_rhd$n, nn_rhd$missmech, nn_rhd$distnum, nn_rhd$varmethod, nn_rhd$k, nn_rhd$rho)
# N x n x missmech x distnum x varmethod x k x rho
table(nn_rhd$sampling)
(2*3*6*2*2*3*10)*13 #PPS
(2*3*6*3*1*3*10)*13 #SRS
dim(nn_rhd)
(2*3*6*3*1*3*10)*13 + (2*3*6*2*2*3*10)*13
ftable(wrhd$sampling, wrhd$N, wrhd$n, wrhd$missmech, wrhd$distnum, wrhd$varmethod, wrhd$k, wrhd$rho)
dim(wrhd)
(2*3*6*2*2*1*10)*13
ftable(bmi$sampling, bmi$N, bmi$n, bmi$missmech, bmi$distnum, bmi$varmethod, bmi$k, bmi$rho)
table(bmi$sampling)
(2*3*6*2*1*2*10)*13 #PPS
(2*3*6*3*1*2*10)*13 #SRS
dim(bmi)
(2*3*6*2*1*2*10)*13 + (2*3*6*3*1*2*10)*13

# combine
all <- as_tibble(bind_rows(nn_rhd,wrhd,bmi))
# all <- as_tibble(bind_rows(nn_rhd,wrhd))
# all <- as_tibble(nn_rhd)

# ##### REMOVE BMI-FPC (k=88) from PPS
# table(all$k, all$sampling)
# all <- all %>% subset(!(k==88 & sampling=="pps"))
# table(all$k, all$sampling)

##### REMOVE WR (varmethod=1) from PPS (use WOR estimator)
table(all$varmethod, all$sampling)
all <- all %>% subset(varmethod!=1)
table(all$varmethod, all$sampling)

# get summary stats for populations (to get true mean)
truth <- read_csv("./data_pop/pop_sumstats.csv", col_types=cols())
long <- truth %>% select(c("N","rho")|starts_with("meany")) %>% pivot_longer(starts_with("meany"),names_to="dist",values_to="truemean") %>% mutate(stat="bias")
long$distnum <- ifelse(long$dist=="meany_norm", 1, ifelse(long$dist=="meany_unif", 2, 3))
all <- left_join(all, long) %>% select(-"dist")

# Relative bias
all$relbias <- 100*all$est/all$truemean
all$relbias_lower <- 100*all$lower/all$truemean
all$relbias_upper <- 100*all$upper/all$truemean

# Separate variables for MAR vs MCAR and percent missing
table(all$missmech)
all$missrate <- ifelse(all$missmech %in% c("mcar25","mar25pos","mar25neg"), 25, 50)
all$mar <- ifelse(all$missmech %in% c("mcar25","mcar50"), 0, 1)

#Factors
all <- all %>% mutate(N_f=factor(N, levels=c(1000,100000), labels=c("N=1,000", "N=100,000")),
                      n_f=factor(n, levels=c(100,200,400), labels=c("n=100","n=200","n=400")), 
                      missrate_f=factor(missrate, levels=c(25,50), labels=c("25% missing","50% missing")),
                      mar_f=factor(mar, levels=c(0,1), labels=c("MCAR","MAR")),
                      missmech_f=factor(missmech, 
                                        levels=c("mcar25","mcar50","mar25pos","mar25neg","mar50pos","mar50neg"), 
                                        labels=c("MCAR/25", "MCAR/50", "MARpos/25", "MARneg/25", "MARpos/50", "MARneg/50")),
                      dist=factor(distnum, levels=1:3, labels=c("Normal","Uniform","Lognormal")),
                      k_f=factor(k, levels=c(0,10,1,5,99,88), labels=c("RHD","wRHD","NN1","NN5","BMI-NN1","BMI-NN1/FPC")),
                      samp_f=factor(sampling, levels=c("srs","pps"), labels=c("SRS","PPS")),
                      varmethod_f=factor(varmethod, levels=c(0,2,9), labels=c("WOR","WOR(Brewer)","n/a")))

save(all, file="./results/all.RData")
