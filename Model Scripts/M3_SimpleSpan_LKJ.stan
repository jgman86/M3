data {
  int <lower=0> N;  // number of subjects
  int <lower=0> K;  // categories 
  int <lower=0> J;  // Dims of Cov Matrix
  int R[K];         // number of responses per category
  int count[N,K];   // observed data
  real scale_b;     // set scaling for background noise
}

parameters {
  // subject parameters
 
  corr_matrix[J] Omega;
  vector<lower=0>[J] sigma;
  vector<lower=0>[J] subj_pars[N];
  vector<lower=0>[J] hyper_pars;
 
}

transformed parameters{

  // activations
  
  vector[K] acts[N];
  real SummedActs[N];
  
  // probabilities
  vector[K] probs[N];
  
   


// loop over subjects to compute activations and probabilites

  for (i in 1:N){ 
    

    acts[i,1] = scale_b + subj_pars[i,1] + subj_pars[i,2];
    acts[i,2] = scale_b + subj_pars[i,2];
    acts[i,3] = scale_b;
    
    SummedActs[i] = R[1] * acts[i,1] + R[2] * acts[i,2] + R[3] * acts[i,3];
    
    probs[i,1] = (R[1] * acts[i,1]) ./ (SummedActs[i]);  
    probs[i,2] = (R[2] * acts[i,2]) ./ (SummedActs[i]);
    probs[i,3] = (R[3] * acts[i,3]) ./ (SummedActs[i]);
    
    
  }
}

model{
 
  
  Omega ~ lkj_corr(2);
  sigma ~ gamma(1,0.01);
    
  hyper_pars[1] ~ normal(20,10);
  hyper_pars[2] ~ normal(2,10);
  
  // Loop over subjects
  
  // Draw subject parameters from multivariate normal with estimated covariance matrix omega
  
  subj_pars[,] ~ multi_normal(hyper_pars,quad_form_diag(Omega, sigma));

  
  for(i in 1:N){
    
    // draw data from probabilities determined by MMM parms
    count[i,]  ~ multinomial(probs[i,]);
  }
}
