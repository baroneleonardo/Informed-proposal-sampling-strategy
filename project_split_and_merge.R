
our_split <- function( j,elem,rho_n ){
  
  output <- list()
    
    l = elem
    
    if (j!=1){
      
      for ( i in 1:(j-1) ){
        l = l - rho_n[i] 
      }
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

