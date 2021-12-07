####################
# NN MI Simulation
# Post-FPBB estimation of the mean (without MI)
# Author: Rebecca Andridge
# Last Modified: 4/20/2021
####################

# Function to calculate estimated mean and variance with Rubin's rules

# ESTS = matrix (LxS) of estimates - L=# BB samples, S=# synth pops per BB sample

fpbbStats_noMI <- function(ESTS)
{
  L <- dim(ESTS)[1]
  Qbar_L <- mean(ESTS)       # overall pt est
  Q_l <- apply(ESTS,1,mean)  # means for each of the L BB samples (averaging over synth pops)
  V_L <- var(Q_l)            # variance of the means of the L synthetic pops
  V_TOTAL <- (1+1/L)*V_L     # overall variance
  DF <- L-1
  LB <- Qbar_L - qt(0.975, DF)*sqrt(V_TOTAL)
  UB <- Qbar_L + qt(0.975, DF)*sqrt(V_TOTAL)
  stats <- as.data.frame(cbind(mean=Qbar_L, var=V_TOTAL, se=sqrt(V_TOTAL), lb=LB, ub=UB, df=DF))
  return(stats)
}
