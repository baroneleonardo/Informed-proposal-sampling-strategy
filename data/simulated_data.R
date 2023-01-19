# # DATA SIMULATION --------

 gamma_sim_1 = 0.5
 gamma_sim_2 = 0.5
 gamma_sim_3 = 0.5

 mu_1 <- c(-2,-2,-2)
 mu_2 <- c(3,3,3)
 mu_3 <- c(-1,-1,-1)


 sigma_1 = sigma_2 = sigma_3 = matrix(0, nrow = 3, ncol = 3)


 diag(sigma_1) = c(0.5,0.5,0.5)

 #sigma_1[1,2]  = 0.2
 #sigma_1[1,3]  = 0.3
 #sigma_1[2,1]  = 0.2
 #sigma_1[3,1]  = 0.3


 diag(sigma_2) <- c(0.5,0.5,0.5)

 diag(sigma_3) <- c(0.5,0.5,0.5)
 #sigma_3[2,3]  = 0.1
 #sigma_3[3,2]  = 0.1

 data_scenario_1 <- as.data.frame(matrix(nrow = 300, ncol = 3))

 data_scenario_1[1,] = mu_1

 for(i in 2:100){

   data_scenario_1[i,] = gamma_sim_1*data_scenario_1[i-1,] + (1-gamma_sim_1)*mu_1 + mvrnorm(n = 1,mu = rep(0,3), Sigma =  sigma_1)

 }


 data_scenario_1[101,] = mu_2

 for(i in 102:200){

   data_scenario_1[i,] = gamma_sim_2*data_scenario_1[i-1,] + (1-gamma_sim_2)*mu_2 + mvrnorm(n = 1,  mu = rep(0,3), Sigma =  sigma_2)

 }

 data_scenario_1[201,] = mu_3

 for(i in 202:300){

   data_scenario_1[i,] = gamma_sim_3*data_scenario_1[i-1,] + (1-gamma_sim_3)*mu_3 + mvrnorm(n = 1, mu = rep(0,3), Sigma =  sigma_3)

}