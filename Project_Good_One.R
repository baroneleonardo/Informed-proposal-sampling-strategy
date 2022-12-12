# FIRST SKETCH


# FIRST SKETCH Z_g(n) FUNCTION USING g = sqrt(x) ------------------------------------------------

Z_sqrt_x <- function(y,rho_n){           # INPUT: DATA y AND PARTITION rho_n (THE ACTUAL ONE OR THE PROPOSAL...)
  
  rho_temp <- rho_n                         # TEMPORAL VARIABLE FOR NEIGHBOURHOODS
  k = length(rho_n)                         # NUMBER OF GROUPS
  out = 0                                   # OUTPUT INIZIALIZATION (IT'S A SUMMATORY)
  
  # 1ST CASE: ONLY THE FIRST GROUP
  
  if( k==1 & rho_n[1] !=1 ){                # IF IT HAS ONLY ONE ELEMENT, DO NOTHING
    for (l in 1:rho_n[1]-1){                # FOR EACH ELEMENT WE DIVIDE pho_n IN 2 GROUPS
      
      rho_temp <- c(l,rho_n[1]-l)           # SPLIT IN TEMPORAL VARIABLE
      
      # COMPUTATION OF POSTERIOR
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      out = out + exp(post_rho_temp)        # UPDATE OF THE SUMMATORY
    }
  }
  
  # 2ND CASE: ONLY THE FIRST AND LAST GROUPS
  
  if ( k==2 ){
    for (l in 1:rho_n[1]){                 # FOR EACH ELEMENT EXCEPT THE LLAST ONE IN FIRST GROUP DO A SPLIT, 
                                           # ON THE LAST ELEMENT DO A MERGE
      ifelse(l != length(rho_n[1]), {rho_temp <- c(l,rho_n[1]-l,rho_n[2])},
             {rho_temp <- c(rho_n[1]+rho_n[2])})
      post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
      out = out + exp(post_rho_temp)       # UPDATE OF THE SUMMATORY
    }
    if(rho_n[2]!=1){
      for (l in 1:rho_n[2]-1){             # FOR EVERY ELEMENT EXCEPT LAST ONE IN LAST GROUP DOO A SPLIT
        rho_temp <- c(rho_n[1],l,rho_n[2]-l)
        post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
        out = out + exp(post_rho_temp)     # UPDATE OF THE SUMMATORY
      }
    }
  }
  
  # 3RD CASE: FIRST GROUP, LAST GROUP AND INTERMEDIARY GROUPS
  
  if( k>2 ){  
    for (j in 1 : length(rho_n)){         # FOR EACH GROUP...
      for (l in 1 : rho_n[j]){            # FOR EACH ELEMENT...
        
        # DIFFERENTIATION FOR THE FIRST GROUP
        
        if( j==1 ){
          
          ifelse(l != length(rho_n[j]), {rho_temp <- c(l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])}, # Split
                 {rho_temp <- c(rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})            # Merge for the last 
          
          post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
          out = out + exp(post_rho_temp)  # UPDATE OF THE SUMMATORY
        }
        
        else{
          
          # DIFFERENTIATION FOR THE LAST GROUP
          
          if (j == k){
            
            if( l!=length(rho_n[j]) ){
              
              rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l)
              post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
              out = out + exp(post_rho_temp) # UPDATE OF THE SUMMATORY
            }
          }
          
          # INTERMEDIARY CASES
          
          else{                      # SPLIT PART
            ifelse(l != length(rho_n[j]), {rho_temp <- c(rho_n[1:(j-1)],l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
                   
                   # MERGE PART WITH DIFFERENTIATION FOR THE SECOND-LAST GROUP
                   {ifelse(j == (length(rho_n) - 1),{rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)])},
                           {rho_temp <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})})
            
            post_rho_temp = posterior(length(rho_temp), gamma_splitting_MULTIVARIATE(y,rho_temp), rho_temp)
            out = out + exp(post_rho_temp)  # UPDATE OF THE SUMMATORY
            
          }
          
          
        }
      }
    }
  }
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n) # MATRIX OF DATA SPLITTED FOR POSTERIOR COMPUTATION
  post_rho = posterior(k,gamma_k,rho_n)            # POSTERIOR OF GIVEN rho_n (IT'S ALREADY IN LOG!!!)
  out = 0.5*(1 - post_rho) + log(out)              # FINAL VALUE OF Z_sqrt_x
  return (out)                       
  
}

# FIRST SKETCH alpha FUNCTION ----------------------------------------------------------------------

alpha <- function(y,rho_n_proposal,rho_n,m_0){ # INPUT: DATA y,NEW POSSIBLE PARTITION rho_n_proposal
                                               # AND THE OLD ONE rho_n, MEAN m_0
  # COMPUTATION OF POSTERIORS
  gamma_k_proposal <- gamma_splitting_MULTIVARIATE(y,rho_n_proposal)
  k_proposal <- length(rho_n_proposal)
  
  gamma_k <- gamma_splitting_MULTIVARIATE(y,rho_n)
  k <- length(rho_n)
  
  post_rho = posterior(k, gamma_k, rho_n)
  post_rho_proposal = posterior(k_proposal, gamma_k_proposal, rho_n_proposal)
  
  # TWO POSSIBILITY IN LOG(SEE FUNCTION Q_f BELOW)
  a1 = post_rho_proposal - post_rho + Q_f(y, rho_n, rho_n_proposal, g, post_rho, post_rho_proposal)
  a2 = log(1)
  
  # DECISION
  alpha = min(a1, a2)
  
  return(alpha)
}

# FIRST SKETCH Q_f FUNCTION ------------------------------------------------------------------
#WRONGGGGGGG
Q_f <- function(y, rho_n, rho_n_proposal, post_rho, post_rho_proposal){
  
  out = Z_sqrt_x(y, rho_n) - Z_sqrt_x(y, rho_n_proposal) + post_rho - post_rho_proposal
  
  return (out)
}