# ------------------------------------------------

posteriors_vector <- function(y,rho_n){              # INPUT: DATA y AND PARTITION rho_n (THE ACTUAL ONE OR THE PROPOSAL...)
  
  rho_temp <- rho_n                         # TEMPORAL VARIABLE FOR NEIGHBOURHOODS
  k = length(rho_n)                         # NUMBER OF GROUPS
  # max_temp = numeric(0)
  post_vector <- as.numeric()
  count = 0
  
  # 1ST CASE: ONLY THE FIRST GROUP
  
  if( k==1 & rho_n[1] !=1 ){                # IF IT HAS ONLY ONE ELEMENT, DO NOTHING
    for (l in 1:(rho_n[1]-1)){                # FOR EACH ELEMENT WE DO A SPLIT
      
      rho_temp <- c(l,rho_n[1]-l)           # SPLIT IN TEMPORAL VARIABLE
      
      # COMPUTATION OF POSTERIOR
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      # max_temp = max(max_temp, post_rho_temp)
      count = count + 1
      post_vector[count] = post_rho_temp
      
    }
  }
  
  # 2ND CASE: ONLY THE FIRST AND LAST GROUPS
  
  if ( k==2 ){
    
    for (l in 1:rho_n[1]){                 # FOR EACH ELEMENT EXCEPT THE LAST ONE IN FIRST GROUP DO A SPLIT, 
      
      ifelse(l != rho_n[1], {rho_temp <- c(l,rho_n[1]-l,rho_n[2])},
             {rho_temp <- c(rho_n[1]+rho_n[2])})
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      # max_temp = max(max_temp, post_rho_temp)
      count = count + 1
      post_vector[count] = post_rho_temp
    }
    
    if(rho_n[2]!=1){
      for (l in 1:(rho_n[2]-1)){             # FOR EVERY ELEMENT EXCEPT LAST ONE IN LAST GROUP DOO A SPLIT
        
        rho_temp <- c(rho_n[1],l,rho_n[2]-l)
        
        post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
        # max_temp = max(max_temp, post_rho_temp)
        count = count + 1
        post_vector[count] = post_rho_temp
        
      }
    }
  }
  
  # 3RD CASE: FIRST GROUP, LAST GROUP AND INTERMEDIARY GROUPS
  
  if( k>2 ){  
    for (j in 1 : length(rho_n)){         
      for (l in 1 : rho_n[j]){            
        
        # DIFFERENTIATION FOR THE FIRST GROUP
        if( j==1 ){
          
          ifelse(l != rho_n[j], {rho_temp <- c(l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])}, # Split
                 {rho_temp <- c(rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})            # Merge for the last 
          
          post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
          # max_temp = max(max_temp, post_rho_temp)
          count = count + 1
          post_vector[count] = post_rho_temp
          
        }
        
        else{
          
          # DIFFERENTIATION FOR THE LAST GROUP
          if (j == k){
            
            if( l!=rho_n[j] ){
              
              rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l)
              
              post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
              # max_temp = max(max_temp, post_rho_temp)
              count = count + 1
              post_vector[count] = post_rho_temp
            }
          }
          
          # INTERMEDIARY CASES
          else{                      # SPLIT PART
            
            ifelse(l != rho_n[j], {rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
                   
                   # MERGE PART WITH DIFFERENTIATION FOR THE SECOND-LAST GROUP
                   {
                     ifelse(j == (length(rho_n) - 1),{rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)])},
                            {rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})
                   })
            
            
            post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
            # max_temp = max(max_temp, post_rho_temp)
            count = count + 1
            post_vector[count] = post_rho_temp
          }
        }
      }
    }
  }
  
  return (post_vector)                       
  
}

# Z_g FUNCTION ----------------------------------------------------------------

Z_sqrt_x <- function (y, rho_n, post_vector){
  
  k = length(rho_n)
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n) # MATRIX OF DATA SPLITTED FOR POSTERIOR COMPUTATION
  post_rho = posterior(k,gamma_k,rho_n)            # POSTERIOR OF GIVEN rho_n
  
  
  sum = 0
  max_temp = max(post_vector)
  
  for (count2 in 1:length(post_vector)){
    sum = sum + exp( 0.5*( post_vector[count2] - max_temp))
  }
  
  out = 0.5*(max_temp - post_rho) + log(sum)              # FINAL VALUE OF Z_sqrt_x
  return (out) 
  
}

# Q DISTRIBUTION --------------------------------------------------------------

Q_distribution <- function (y, rho_n, post_vector){      # Return a sample (1 elem)  from the distribution Q(n,n^i)
  
  k = length(rho_n)
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n) 
  post_rho = posterior(k,gamma_k,rho_n)
  
  q_dist <- as.numeric()
  
  for (i in 1:length(post_vector)){
    q_dist[i] = 0.5*(post_vector[i] - post_rho)   # Computing the q_probability vector (still in log)
  }
  
  q_dist = exp(q_dist)                            # Returning in non-log probability 
  s = sum(q_dist)                          
  q_dist = q_dist/s                               # Normalization of the vector
  
  elem = sample(1:(n-1), 1, prob = q_dist)        # Returning one sample whit prob. the q_distribution
  
  return(elem)
}


# Check if merge or split ----------------------------------------------------

merge_or_split <- function(rho_n, elem){ 
  
  n_elem = 0
  
  for (j in 1:(length(rho_n) - 1)){
    
    n_elem = n_elem + rho_n[j]
    
    if (elem == n_elem){return(our_merge(j,rho_n))}
    
    if (elem < n_elem){return(our_split(j, elem, rho_n))}
  }
  
  return(our_split(length(rho_n), elem, rho_n))
}


# Q_fraction ------------------------------------------------------------------


Q_fraction <- function(y, rho_n, rho_n_proposal, post_vector_1, post_vector_2, post_rho,  post_rho_proposal){
  
  out = Z_sqrt_x(y, rho_n, post_vector_1) - Z_sqrt_x(y, rho_n_proposal, post_vector_2) + post_rho - post_rho_proposal
  
  return (out)
}



# ALPHA FUNCTION ----------------------------------------------------------------------

our_alpha <- function(y, rho_n_proposal, rho_n, post_vector_1, post_vector_2, m_0){ # INPUT: DATA y,NEW POSSIBLE PARTITION rho_n_proposal
  
  gamma_k_proposal <- gamma_splitting_MULTIVARIATE(y,rho_n_proposal)
  k_proposal <- length(rho_n_proposal)
  
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  k <- length(rho_n)
  
  post_rho = posterior(k, gamma_k, rho_n)
  post_rho_proposal = posterior(k_proposal, gamma_k_proposal, rho_n_proposal)
  
  a1 = post_rho_proposal - post_rho + Q_fraction(y, rho_n, rho_n_proposal, post_vector_1, post_vector_2, post_rho, post_rho_proposal)
  a2 = log(1)
  
  alpha = min(a1, a2)
  
  return(alpha)
}

# SPLIT function --------------------------------------------------------------

our_split <- function( j,elem,rho_n ){
  
  output <- list()
  l = elem
  
  if (j!=1){
    
    for ( i in 1:(j-1) ){ l = l - rho_n[i] }
  }
  
  if(j != 1){
    
    # Split the j-group in [l; length(j-group)-l]
    ifelse(j != length(rho_n), {rho_n <- c(rho_n[1:(j-1)],l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
           {rho_n <- c(rho_n[1:(j-1)],l,rho_n[j]-l)}) # If I'm in the last group
  }
  
  # If I'm in the first group
  if(j == 1){
    
    ifelse(j != length(rho_n), {rho_n <- c(l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
           {rho_n <- c(l,rho_n[j]-l)})
  }
  
  output = rho_n      # First group of the split
  
  return(output)
  
}

# MERGE function --------------------------------------------------------------

our_merge <- function(j,rho_n){
  
  output <- list()
  c <- length(rho_n)
  
  if (c != 1){
    
    if (c != 2){
      
      ifelse(j != 1, {
        ifelse(j == (length(rho_n) - 1),{rho_n <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)])},
               {rho_n <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})
      },
      {rho_n <- c(rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})
    }
    
    if (c == 2){
      j <- 1
      rho_n <- c(rho_n[j] + rho_n[(j+1)])
    }
  } 
  
  output = rho_n
  
  return(output)
  
}

# Create a RANDOM INITIAL PARTITION from the dataset --------------------------

random_partition <- function(n){
  
  # Set random initial partition
  time_indexes = 1:n
  k = sample(1:round(n/10, digits = 1), size = 1) # Number of change-points
  changepoint_positions = sort( sample(1:(n-1), size = k, replace = F) ) # Sample change-points locations
  
  # Convert the previous object in rho_n_0, the vector of cardinalities of the groups
  rho_n_0 = c(changepoint_positions, n)
  rho_n_0[2:(k+1)] = rho_n_0[2:(k+1)] - rho_n_0[1:k]  
  return(rho_n_0)
  
}




