# Abbiamo un po' iniziato a scrivere il codice

# Z_g(n) Normalization Costant -----------------------------------------------

Z_g <- function(y,rho_n, g){   # Da mettere il log
  
  k = length(rho_n)
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  post_rho = posterior(k,gamma_k,rho_n)
  rho_temp <- rho_n
  out = 0
  
  
  
  
  # CASO UN SOLO GRUPPO
  if( k==1 ){
    for (l in 1:rho_n[1]-1){
      
      rho_temp <- c(l,rho_n[1]-l)
      
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      out = out + g(post_rho_temp/post_rho)
    }
  }
  
  if ( k==2 ){
    for (l in 1:rho_n[1]){
      ifelse(l != length(rho_n[1]), {rho_temp <- c(l,rho_n[1]-l,rho_n[2])},
             {rho_temp <- c(rho_n[1]+rho_n[2])})
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      out = out + g(post_rho_temp/post_rho)
    }
    for (l in 1:rho_n[2]-1){
      rho_temp <- c(rho_n[1],l,rho_n[2]-l)
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      out = out + g(post_rho_temp/post_rho)
    }
  }
  
  #  CASO CON GRUPPI INTERMEDI
  if( k>2 ){  
    for (j in 1 : length(rho_n)){
      for (l in 1 : rho_n[j]){
        
        #  SPLIT E MERGE DIFFERENTE PER IL PRIMO GRUPPO
          if( j==1 ){
              
              ifelse(l != length(rho_n[j]), {rho_temp <- c(l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])}, # Split
                     {rho_temp <- c(rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])}) # Merge for the last 
              
              post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
              out = out + g(post_rho_temp/post_rho)
          }
          
          else{
            if (j == k){
              
              if( l!=length(rho_n(j)) ){
                
                rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l)
                post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
                out = out + g(post_rho_temp/post_rho)
              }
            }
            else{
              ifelse(l != length(rho_n[j]), {rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
                     
                     {ifelse(j == (length(rho_n) - 1),{rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)])},
                             {rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})})
              
              post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
              out = out + g(post_rho_temp/post_rho)
              
            }
          
            
          }
      }
    }
  }
  out = 1/(n-1) * out
  return (out)
  
  }
