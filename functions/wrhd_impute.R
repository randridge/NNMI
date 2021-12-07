####################
# NN MI Simulation
# Function to impute using weighted RHD with ABB
# Author: Rebecca Andridge
# Last Modified: 4/23/2021
####################

# Y = vector with missing values
# W = sample weights
# D = number of imputations

wrhd_impute <- function(Y,W,D)
{
  m <- is.na(Y)             # missing indicator
  y_donors <- Y[m==0]       # y values for complete observations
  w_donors <- W[m==0]       # weights for complete observations
  id_donors <- 1:sum(1-m)   # sequential ID values for complete observations
  yimp <- matrix(Y, nrow=length(Y), ncol=D)  # matrix to hold imputed values
  for (d in 1:D)
  {
    id_donors_abb <- sample(id_donors, replace=T)      # ABB: bootstrap donors
    y_donors_abb <- y_donors[id_donors_abb]            # y values for donors
    w_donors_abb <- w_donors[id_donors_abb]            # weight values for donors
    prob_donors_abb <- w_donors_abb/sum(w_donors_abb)  # probabilities proportional to W
    # Fill in matrix of imputed values
    yimp[m==1,d] <- sample(y_donors_abb, prob=prob_donors_abb, size=sum(m), replace=T)
  }
  return(yimp)
}