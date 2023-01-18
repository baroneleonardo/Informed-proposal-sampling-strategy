
random_partition <- function(n){
  
  ### Set random initial partition
  time_indexes = 1:n
  k = sample(1:round(n/10, digits = 1), size = 1) # number of change-points
  changepoint_positions = sort( sample(1:(n-1), size = k, replace = F) ) #sample change-points locations
  
  # convert the previous object in rho_n_0, the vector of cardinalities of the groups
  rho_n_0 = c(changepoint_positions, n)
  rho_n_0[2:(k+1)] = rho_n_0[2:(k+1)] - rho_n_0[1:k]  
  return(rho_n_0)
  
}