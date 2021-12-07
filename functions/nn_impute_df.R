####################
# NN MI Simulation
# Function to impute using NN (or variant) with ABB
# Author: Rebecca Andridge
# Last Modified: 4/28/2021
####################

# DATA = dataframe with Y,X where Y has missing values and X does not
# K = size of donor pool

nn_impute_df <- function(DATA,K)
{
  m <- is.na(DATA$y)             # missing indicator
  x_recipients <- DATA$x[m==1]   # x values for incomplete observations
  x_donors <- DATA$x[m==0]       # x values for complete observations
  y_donors <- DATA$y[m==0]       # y values for complete observations
  id_donors <- 1:sum(1-m)        # sequential ID values for complete observations
  id_donors_abb <- sample(id_donors, replace=T)  # ABB: bootstrap donors
  x_donors_abb <- x_donors[id_donors_abb]        # x values for donors
  y_donors_abb <- y_donors[id_donors_abb]        # y values for donors
  # Impute based on value of K
  if (K==0) # use all respondents as donor pool
  {
    idx <- sample(id_donors, size=sum(m), replace=T)
  } else {
    idx <- mice:::matchindex(x_donors_abb, x_recipients, k=K) # randomly select donors from pool of K closest donors
  }
  # Fill in matrix of imputed values
  impdat <- DATA
  impdat$y[m==1] <- y_donors_abb[idx]
  return(impdat)
}