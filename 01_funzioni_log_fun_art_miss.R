
# SPLIT FUNCTION ----------------------------------------------------------

split <- function(k,rho_n){
  
  output <- list()
  
  z <- which(rho_n > 1) # Return a vector of the position of the group with element > 1
  
  if (length(z) != 0 ){
    
   # If I have more than one group I choose with sample where to split 
    
  ifelse(length(z) != 1, {j <- sample(z,size = 1)}, {j <- z}) # we need to consider the case of a single partition sample
  
    l <- sample(1:(rho_n[j]-1),1) # Choose an element in the j-group (exepet the last)
    
    if(j != 1){
      
      # Split the j-group in [l; length(j-group)-l]
      
      ifelse(j != length(rho_n), {rho_n <- c(rho_n[1:(j-1)],l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
           {rho_n <- c(rho_n[1:(j-1)],l,rho_n[j]-l)}) # If I'm in the last group
    
    k <- k + 1
      
    }
    # If I'm in the first group
    if(j == 1){
      
      ifelse(j != length(rho_n), {rho_n <- c(l,rho_n[j]-l,rho_n[(j+1):length(rho_n)])},
             {rho_n <- c(l,rho_n[j]-l)})
      
    }
  }
  
  output[[1]] = rho_n
  output[[2]] = j      # First group of the split
  
  return(output)
  
}

# MERGE FUNCTION ----------------------------------------------------------

merge <- function(k,rho_n){
  
  output <- list()
  
  c <- length(rho_n)
  
  if (c != 1){
    
    if (c != 2){
    
   j <- sample(1:(k-1),1)
  
   u <- runif(1) # ???
  
      ifelse(j != 1, {
        ifelse(j == (length(rho_n) - 1),{rho_n <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)])},
                            {rho_n <- c(rho_n[1:(j-1)],rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})
                      },
        
               {rho_n <- c(rho_n[j] + rho_n[(j+1)], rho_n[(j+2):length(rho_n)])})
  
      k <- k - 1 
  }
  
     if (c == 2){
      
       j <- 1
      
       rho_n <- c(rho_n[j] + rho_n[(j+1)])
        
       k <- length(rho_n) 
        
     }
    
  } 
  
  output[[1]] = rho_n
  output[[2]] = j
  
  return(output)
  
}

# SHUFFLE FUNCTION --------------------------------------------------------

shuffle <- function(k,rho_n){
  
  i <- sample(1:(k-1),1)
  
  j <- sample(1:(rho_n[i] + rho_n[i + 1] - 1),1)
    
  rho_n[i+1] <- (rho_n[i] + rho_n[i+1] - j)
    
  rho_n[i] <- j 
  
  return(rho_n)
}

# VARIOUS FUNCTIONS ----------------------------------------------------------

#Funzione Indicatrice

indicator_1 <- function(k,n) ifelse(0 < k & k < n,1,0)

indicator_2 <- function(k) ifelse(k == 1,1,0)

#Shifted Gamma

r_shifted_gamma <- function(u,o,sigma_shift){
  rgamma(1,u,o) + sigma_shift
}

#Absolute Stirling Number 1st kind

abs_stirling_number_1st <- function(r,k1){
  
  ifelse(k1 > 0 & k1 <= r, {abs_str_num <- ((-1)^(r-k1))*as.numeric(Stirling1(r,k1))}, {abs_str_num <- 0})
  
  if(r == 0 & k1 == 0){abs_str_num <- 1}
  
  return(abs_str_num)
  
}

#Variation of Information for final partition

VI <- function(c1, c2){
  n <- length(c1)
  
  T1 <- table(c1) / n
  T2 <- table(c2) / n
  
  H1 <- - sum(T1 * log(T1))
  H2 <- - sum(T2 * log(T2))
  
  T12 <- table(c1, c2) / n
  I12 <- sum(T12 * log(T12 / T1 %*% t(T2)), na.rm = T)
  
  return((H1 + H2 - 2 * I12) / log(n))
}



# PRIOR FUNCTIONS ON THE PARAMETERS --------------------------------------------



# GAMMA 

full_conditional_gamma <- function(n, k, rho_n, sigma, theta, gamma, gamma_k){
  
  
  ##--##
  
  d     =  nbasis = 3
  phi_0 = phi_0
  nu_0  =  nu_0 
  m_0   = m_0
  k_0   =  k_0

  ##--##  
  
  vec_num <- as.numeric()
  
  for(j in 1:(k-1)){
    
    vec_num[j] = theta + j*sigma
    
  }
  

  
  num <- lfactorial(n) + (k*d/2)*log(k_0) + (k*nu_0/2) * log(det(phi_0)) + sum(log(vec_num))
  
  den <- lfactorial(k) + log_pochhammer((theta+1),(n-1)) + (n*d/2) * log(pi) + ((n-k)/2)*log(1-gamma^2) + k * lgamma_d_a(d,(nu_0/2))
  
  vec_dx <- as.numeric()
  
  for(j in 1:k){
    
    n_j = rho_n[j]
    
    k_nj = compute_k_n(k_0, gamma, n_j)
      
    gamma_k_j = as.matrix(as.data.frame(gamma_k[j]))

    S_nj = compute_s_n(d, m_0, k_0, nu_0, phi_0, gamma, n_j, gamma_k_j, nbasis)
    
    ###
    
    num_vec_dx = lgamma_d_a(d,((nu_0 + n_j)/2)) + log_pochhammer((1-sigma),(n_j-1))  
    
    den_vec_dx = lfactorial(n_j) + (d/2)*(k_nj) + ((nu_0 + n_j)/2) * log(det(S_nj))
    
    vec_dx[j] = num_vec_dx - den_vec_dx
    
  }
  
  res <- num - den + sum(vec_dx) 
  
  return(res)
  
  
}





# SIGMA

full_conditional_sigma <- function(sigma,theta,k,rho_n,a,b,c,d){
  
  vec_prod_1 <- as.numeric(0)
  
  for(i in 1:(k-1)){
    vec_prod_1[i] = log(theta + i*sigma)
  }
  
  produttoria_1 <- sum(vec_prod_1)
  
  vec_prod_2 <- as.numeric(0)
  
  for(i in 1:k){
    vec_prod_2[i] = log_pochhammer((1-sigma),(rho_n[i]-1))
  }
  
  produttoria_2 <- sum(vec_prod_2)
  
  res <- (a-1) * log(sigma) + (b-1) * log(1-sigma) + (c-1)*log(theta + sigma) + log(exp(-d*sigma)) + produttoria_1 + produttoria_2
  
  return(res)
}




# THETA 

full_conditional_theta <- function(prior_c,prior_d,candidato,k){
  
  vecc <- as.numeric()
  
  z <- rbeta(1,candidato + 2, n)
  
  f <- rexp(1,candidato + 1)
  
  for(x in 0:(k+1)){
    
    omega_j_NUM <-((n-sigma)*(n+1-sigma)*as.numeric(abs_stirling_number_1st(k-1,x)) + (2*n + 1 - 2*sigma)*sigma*as.numeric(abs_stirling_number_1st(k-1,x-1)) + sigma^2 * abs_stirling_number_1st(k-1,x-2)) * gamma(prior_c+x)
    
    omega_j_DEN <- (sigma*(prior_d + f - log(z)))^x
    
    omega_j <- omega_j_NUM/omega_j_DEN
    
    vecc[x+1] <- omega_j
    
  }
  
  vecc = vecc/sum(vecc)
  
  u = runif(1)
  
  component <- min(which(cumsum(vecc) > u))
  
  theta_j <- r_shifted_gamma(prior_c+ (component-1), prior_d + f -log(z), sigma)
  
  return(theta_j)
  
}
