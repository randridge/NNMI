####################
# NN MI Simulation
# Function to impute using NN or RHD with ABB
# Author: Rebecca Andridge
# Last Modified: 4/22/2021
####################

# Y = vector with missing values
# X = fully observed covariate
# K = size of donor pool (0=whole cell, i.e., RHD)
# D = number of imputations

nn_impute <- function(Y,X,K,D)
{
  m <- is.na(Y)             # missing indicator
  x_recipients <- X[m==1]   # x values for incomplete observations
  x_donors <- X[m==0]       # x values for complete observations
  y_donors <- Y[m==0]       # y values for complete observations
  id_donors <- 1:sum(1-m)   # sequential ID values for complete observations
  yimp <- matrix(Y, nrow=length(Y), ncol=D)  # matrix to hold imputed values
  for (d in 1:D)
  {
    id_donors_abb <- sample(id_donors, replace=T)  # ABB: bootstrap donors
    x_donors_abb <- x_donors[id_donors_abb]        # x values for donors
    y_donors_abb <- y_donors[id_donors_abb]        # y values for donors
    # Impute based on value of K
    if (K==0) # use all respondents as donor pool (RHD) - respondent Y values are bootstrapped
    {
      idx <- sample(id_donors, size=sum(m), replace=T) 
    } else {  # impute either K=1 closest or random draw from K closest
      idx <- mice:::matchindex(x_donors_abb, x_recipients, k=K) # randomly select donors from pool of K closest donors
    }
    # Fill in matrix of imputed values
    yimp[m==1,d] <- y_donors_abb[idx]
  }
  return(yimp)
}