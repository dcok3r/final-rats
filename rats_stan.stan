
// data block: these variables must be fed to stan as a list from R
data {
  int<lower=0> N;    // sample size
  int<lower=0> K;   // number of parameters
  vector[N] y;       // response vector
  matrix [N,K] X;  // model matrix of predictors. nrow=N ncol=K
  
  int<lower=0> n_sp;    // number of sp
  array [N] int<lower=1,upper=n_sp> species;   // index of sp for individual

  // predictors
  //array [N] int <lower=0,upper=2> distance;    
  //array [N] int <lower=0,upper=1> macrotis; 
  //array [N] int <lower=0,upper=1> male;
  //interaction terms
  //array [N] int <lower=0,upper=2> macrotis_distance;    
  //array [N] int <lower=0,upper=1> macrotis_male; 
  //array [N] int <lower=0,upper=2> male_distance;
}

transformed data {

}

parameters {
  vector[n_sp] alpha; // intercepts for each sp
  vector[K] beta;
  real<lower=0.0> sigma;
}

transformed parameters {
  //vector[N] mu; // predicted mean values
  //for (x in 1:N) {
    //mu[x] = ;
                
      // X*beta          
    //}
}

model {
  alpha~normal(0,1);
  beta~normal(0,1);
 
  sigma ~ exponential(1);
 
  y ~ normal(alpha[species] + X*beta, sigma); 
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | X[n] * beta, sigma);
  }
}
