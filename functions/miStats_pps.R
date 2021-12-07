####################
# NN MI Simulation
# Post-MI estimation of the mean for PPS
# Author: Rebecca Andridge
# Last Modified: 5/1/2021
####################

# Function to calculate estimated mean and variance with Rubin's rules

# MIDATA = matrix of imputed data, with columns = imputations
# PIK = vector of inclusion probs from PPS

miStats_pps <- function(MIDATA, PIK)
{
  D <- ncol(MIDATA)
  DAT <- as.data.frame(cbind(MIDATA, PIK))
  # with-replacement approximation
  des <- svydesign(id=~1, probs=~PIK, data=DAT)
  svyest <- svymean(~MIDATA, design=des)
  thetaHats <- coef(svyest)
  Ws <- diag(vcov(svyest))
  thetaHat <- mean(thetaHats)
  W <- mean(Ws)
  B <- (1/(D-1))*sum((thetaHats-thetaHat)^2)
  T <- W + ((D+1)/D)*B
  df <- (D-1)*(1 + (D/(D+1))*(W/B))^2
  FMI <- ((D+1)/D)*B/T
  lb <- thetaHat - qt(0.025, df, lower.tail=F)*sqrt(T)
  ub <- thetaHat + qt(0.025, df, lower.tail=F)*sqrt(T)
  stats_wr <- as.data.frame(cbind(mean=thetaHat, var=T, lb, ub, df, FMI, B=B, W=W))
  # brewer WOR approximation
  des <- svydesign(id=~1, fpc=~PIK, data=DAT, pps="brewer")
  svyest <- svymean(~MIDATA, design=des)
  thetaHats <- coef(svyest)
  Ws <- diag(vcov(svyest))
  thetaHat <- mean(thetaHats)
  W <- mean(Ws)
  B <- (1/(D-1))*sum((thetaHats-thetaHat)^2)
  T <- W + ((D+1)/D)*B
  df <- (D-1)*(1 + (D/(D+1))*(W/B))^2
  FMI <- ((D+1)/D)*B/T
  lb <- thetaHat - qt(0.025, df, lower.tail=F)*sqrt(T)
  ub <- thetaHat + qt(0.025, df, lower.tail=F)*sqrt(T)
  stats_wor <- as.data.frame(cbind(mean=thetaHat, var=T, lb, ub, df, FMI, B=B, W=W))
  # combina
  stats <- cbind(varmethod=1:2, rbind(stats_wr,stats_wor))
  return(stats)
}
