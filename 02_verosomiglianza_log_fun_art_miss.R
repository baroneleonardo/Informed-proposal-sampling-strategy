# POCHAMMER

log_pochhammer <- function(x,n){   
  
  num_vec <- as.numeric()
  num_vec[1] = x
  
  if (n!=0){
    for(i in 1:(n-1)){num_vec[i+1] <- (x+i)} 
  }
  
  if (n==0){num_vec[1] <- 1}
  
  num <- sum(log(num_vec))

  return(num)
}


pochammer <- function(x,n){
  
  num_vec <- as.numeric()
  
  num_vec[1] = x
  
  if (n!=0){
    for(i in 1:(n-1)){
      num_vec[i+1] <- (x+i)
    } 
  }
  
  if (n==0){num_vec[1] <- 1}
  
  num <- prod(num_vec)
  
  return(num)
  
} 

# 

prod_1_prior <- function(k,theta,sigma){   
  
  num <- log(theta + 1*sigma)
  
  for(i in 2:(k-1)){
    
    num <- num + log(theta + i*sigma)
    
  }
  return(num)
}

prod_2_prior <- function(k,sigma,rho_n){
  vec <- as.numeric()
  for (j in 1:k){
    vec[j] <- log_pochhammer((1-sigma),(n_s(j,rho_n)-1)) - lfactorial(n_s(j,rho_n))
  }
  out <- sum(vec)
  return(out)
}


# MULTIVARIATE GAMMA FUNCTION

gamma_d_a <- function(d,a){
  
  vettore_gamma <- as.numeric()
  
  for (j in 1:d){vettore_gamma[j] <- gamma(a + ((1-j)/2))}
  
  ret <- pi^((d*(d-1))/4) * prod(vettore_gamma)
  
  return(ret)
}

# LOG MULTIVARIATE GAMMA FUNCTION

lgamma_d_a <- function(d,a){
  
  vettore_gamma <- as.numeric()
  
  for (j in 1:d){vettore_gamma[j] <- lgamma(a + ((1-j)/2))}
  
  ret <- log(pi^((d*(d-1))/4)) + sum(vettore_gamma)
  
  return(ret)
}

# S MATRIX 

theta_matrix <- function(j,theta,rho_n){
  
  dim <- n_s(j,rho_n)
  
  if (dim != 1){
    
  matrix <- matrix(data = 0, nrow = dim, ncol = dim)
  
  diag(matrix) <- 1 + theta^2
  
  matrix[1,1] <- 1
  matrix[dim,dim] <- 1
  
    for (i in 1:(dim-1)){
      matrix[i+1,i] <- -theta
      matrix[i,i+1] <- -theta
    }    
  }
  
  if (dim == 1){matrix <- 1}    
  
  return(matrix)
}

# PRIOR ON THE PARTITIONS --------------------------------------------------

prior <- function(n,k,theta,sigma,rho_n){
  
  a1 <- lfactorial(n) - lfactorial(k)
  
  a2 <- prod_1_prior(k,theta,sigma) - log_pochhammer((theta + 1), (n - 1))
  
  a3 <- prod_2_prior(k,sigma,rho_n)
  
  return(a1 + a2 + a3)
}


# LIKELIHOOD --------------------------------------------------------------

likelihood <- function(gamma_k,k,rho_n){
  
  vettore_veros <- as.numeric(0)
  
  nbasis = n_col # Insert the number of columns of Y 
  
  k_0 = k_0 
  
  V_0 = nu_0 
  
  phi_0 = phi_0
  
  m_0 = m_0
  
  for(j in 1:k){
    
    d <- ncol(as.data.frame(gamma_k[j]))
    
    gamma_k_j = as.matrix(as.data.frame(gamma_k[j]))
    
    n_gamma <- nrow(as.data.frame(gamma_k[j]))
    
    vettore_veros[j] = inner_likelihood(d,m_0,k_0,V_0,phi_0,gamma, n_gamma, gamma_k_j, nbasis)
    
  }
  
  likelihood = sum(vettore_veros)
  
  return(likelihood)
  
}

# POSTERIOR ---------------------------------------------------------------


posterior <- function(k,gamma_k,rho_n,m_0){
  
  post = prior(n,k,theta,sigma,rho_n) + likelihood(gamma_k,k,rho_n)
  
  return(post)
  
}
  
  