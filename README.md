# Project of Bayesian Statistics Course 

In 2017, Zanella (*"Informed Proposals for Local MCMC in Discrete Spaces"*)  proposed an original approach naturally applicable to discrete spaces without destroying their nature using a particular class of locally balanced pointwise informed proposal functions. This approach showed improvements in efficiency compared to uninformed schemes.
In our work, we combine Zanella  informed model with the multivariate change point detection model of Corradin, Danese and Ongaro (*"Bayesian nonparametric change point detection for multivariate time series with missing observations"*) and we test this new strategy on both simulated and real data.


The process is divided in files as following: 

  - files *"01_funzioni_log_fun_art_miss"*, *02_verosomiglianza_log_fun_art_miss"* and *"03_alpha_log_fun_art_miss"* contain the R functions 
  - *"wade.cpp"* and *"cpp_funz.cpp"* contain the C++ functions 
  - *"043_algoritmo_log_fun_art_miss_50rip"* is the main code 
  - *"project_functions.R"* contains the code of our project
  
The dataset are available in *"data"* folder:
  - *"dpc-covid19-ita-regioni"*" contains data from the COVID19 application
  - *"read_data"* contains the code to read and transform COVID19 data
  - *"exchange_rate"* contains data from the EU/CHF, EU/GBP, EU/USD exchange rate
  - *"exchange_data"* contains the code to read and transform EXCHANGE RATE
  -*"simulated_data"* for synthetic data 
