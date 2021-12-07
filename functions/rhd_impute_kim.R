####################
# NN MI Simulation
# Function to impute using RHD modifications with ABB
# Author: Rebecca Andridge
# Last Modified: 6/11/2020
####################

# Y = vector with missing values
# X = fully observed covariate
# D = number of imputations

rhd_impute_kim <- function(Y,D)
{
  # KIM 2002 adjustment - change the size of the resampled donor pool
  # Determine size to resample
  m <- is.na(Y)
  n <- length(Y)
  r <- sum(1-m)
  r_kim <- round((r-1)*(n-r-1)*(n-2) / ((n-1)*(n-r-1)+n+r-1))  # adjusted size
  # Impute
  yimp <- matrix(Y, nrow=length(Y), ncol=D)
  for (d in 1:D)
  {
    # Bootstrap observed values, taking only r_kim of them
    y_donors_abb <- sample(Y[!m], size=r_kim, replace=TRUE)
    # Randomly select donors from entire resampled pool
    yimp[m==1,d] <- sample(y_donors_abb, size=sum(m), replace=TRUE)
  }
  # Return matrix of imputations (cols=MI data, rows=missing observations)
  return(yimp)
}