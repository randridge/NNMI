####################
# NN MI Simulation
# Post-MI estimation of the mean for SRS (includes FPC)
# Author: Rebecca Andridge
# Last Modified: 6/11/2020
####################

# Function to calculate estimated mean and variance with Rubin's rules

# MIDATA = matrix of imputed data, with columns = imputations

miStats_srs <- function(MIDATA, FPC)
{
  D <- ncol(MIDATA)
  thetaHats <- apply(MIDATA, 2, mean)
  Ws <- apply(MIDATA, 2, function(x) var(x)/nrow(MIDATA)) * FPC
  thetaHat <- mean(thetaHats)
  W <- mean(Ws)
  B <- (1/(D-1))*sum((thetaHats-thetaHat)^2)
  T <- W + ((D+1)/D)*B
  df <- (D-1)*(1 + (D/(D+1))*(W/B))^2
  FMI <- ((D+1)/D)*B/T
  lb <- thetaHat - qt(0.025, df, lower.tail=F)*sqrt(T)
  ub <- thetaHat + qt(0.025, df, lower.tail=F)*sqrt(T)
  stats <- as.data.frame(cbind(mean=thetaHat, var=T, lb, ub, df, FMI, B=B, W=W))
  return(stats)
}
