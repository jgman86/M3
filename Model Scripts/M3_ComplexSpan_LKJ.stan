// Implemented in Stan by Jan Goetttmann 
// M3 - Model after Oberauer & Lewandowsky, 2019
// Questions regarding the mode code to Jan.Goettmann@psychologie.uni-heidelberg.de



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
  vector [J] hyper_pars;
  vector [J] subj_pars[N];
  
  
}

transformed parameters{
  
  
  // Transform f Parameter
  
  real f[N] = inv_logit(subj_pars[,3]);
  real mu_f = inv_logit(hyper_pars[3]);

  
  // activations
  vector[K] acts[N];
  real SummedActs[N];
  
  // probabilities
  
  vector[K] probs[N];
  
  
  // Transformations
  
  
  // loop over subjects to compute activations and probabilites
  
  for (i in 1:N){
    acts[i,1] = scale_b + subj_pars[i,1] + subj_pars[i,2]; // Item in Position
    acts[i,2] = scale_b + subj_pars[i,2];        // Item in Other Position
    acts[i,3] = scale_b + f[i]* (subj_pars[i,1]+subj_pars[i,2]);// Distractor in Position
    acts[i,4] = scale_b + f[i]*subj_pars[i,2]; // Distractor in other Position
    acts[i,5] = scale_b; // non presented Lure
    
    SummedActs[i] = R[1] * acts[i,1] + R[2] * acts[i,2] + R[3] * acts[i,3]+ R[4] * acts[i,4]+ R[5] * acts[i,5];
    
    probs[i,1] = (R[1] * acts[i,1]) ./ (SummedActs[i]);  
    probs[i,2] = (R[2] * acts[i,2]) ./ (SummedActs[i]);
    probs[i,3] = (R[3] * acts[i,3]) ./ (SummedActs[i]);
    probs[i,4] = (R[4] * acts[i,4]) ./ (SummedActs[i]);
    probs[i,5] = (R[5] * acts[i,5]) ./ (SummedActs[i]);
    
    
  }
  
}


model{
  
  
  // priors for hyper parameters
  
  Omega ~ lkj_corr(2);
  sigma ~ gamma(1,0.01);
  
  hyper_pars[1] ~ normal(20,10);
  hyper_pars[2] ~ normal(2,10);
  hyper_pars[3] ~ normal(0,10);

  

  
  // Sample Parameters
  

  subj_pars[,] ~ multi_normal(hyper_pars,quad_form_diag(Omega, sigma));
  
  
  // Loop over subjects
  for(i in 1:N){
    
    // draw data from probabilities determined by MMM parms
    count[i,]  ~ multinomial(probs[i,]);
  }
}

