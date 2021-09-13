// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de

data {
  int <lower=0> N;               // number rownumber
  int <lower=0> K;              // Number of Response categories 
  int <lower=0> J;             // Dims of Cov Matrix
  int <lower=0> Con1;          // Levels of Condition 1 
  int <lower=0> Con2;         // Levels of Condition 2
  int R[K];                  // Number of responses per category
  int count[N*Con1*Con2,K]; // Observed Data
  real scale_b;            // Set scaling for background noise
  
  vector[N*Con1*Con2] t_eU;     // Cue-Word-Intervall - extended Updating benefit
  vector[N*Con2*Con2] t_rm;    // New Word-Cue Intervall - removal benefit after after new presented Item for old item
}

parameters {
  // Defining vector for hyper and subject parameters 
  
  cholesky_factor_corr[J] Lcorr;
  vector<lower=0>[J] sigma;
  vector [J] hyper_pars;
  vector [J] subj_pars[N];
  
}

transformed parameters{
  // Transform d Parameter
  
  real d[N] = inv_logit(subj_pars[,5]);
  real mu_d = inv_logit(hyper_pars[5]);
  real sig_d = sd(d);
  
  
  // Activations
  real acts_IIP[N*Con1*Con2];
  real acts_IOP[N*Con1*Con2];
  real acts_DIP[N*Con1*Con2];
  real acts_DIOP[N*Con1*Con2];
  real acts_NPL[N*Con1*Con2];


  // probabilities
  vector[K] probs[N*Con1*Con2];
  real SummedActs[N*Con1*Con2];

  
  // loop over subjects and conditions to compute activations and probabilites
  
  for (i in 1:N){ // for each subject
    for(j in 1:Con1*Con2) {
     
      
    acts_IIP[j + (i-1)*Con1*Con2] = scale_b + (1+subj_pars[i,4]*t_eU[j])*(subj_pars[i,2] + subj_pars[i,1]); // Item in Position                      
    acts_IOP[j + (i-1)*Con1*Con2] = scale_b + (1+subj_pars[i,4]*t_eU[j])*subj_pars[i,2];        // Item in Other Position
    acts_DIP[j + (i-1)*Con1*Con2] = scale_b + exp(-subj_pars[i,3]*t_rm[j])*d[i]*(1+subj_pars[i,4]*t_eU[j])*(subj_pars[i,2]+subj_pars[i,1]);// Old Item in Position
    acts_DIOP[j + (i-1)*Con1*Con2] = scale_b + exp(-subj_pars[i,3]*t_rm[j])*d[i]*(1+subj_pars[i,4]*t_eU[j])*subj_pars[i,2]; // Old item in other Position
    
    
    acts_NPL[j + (i-1)*Con1*Con2] = scale_b; // non presented Lure
    
    SummedActs[j + (i-1)*Con1*Con2] = R[1] * acts_IIP[j + (i-1)*Con1*Con2] + R[2] * acts_IOP[j + (i-1)*Con1*Con2] + R[3] * acts_DIP[j + (i-1)*Con1*Con2]+
                      R[4] * acts_DIOP[j + (i-1)*Con1*Con2]+ R[5] * acts_NPL[j + (i-1)*Con1*Con2];
    
    probs[j + (i-1)*Con1*Con2,1] = (R[1] * acts_IIP[j + (i-1)*Con1*Con2]) ./ (SummedActs[j + (i-1)*Con1*Con2]);  
    probs[j + (i-1)*Con1*Con2,2] = (R[2] * acts_IOP[j + (i-1)*Con1*Con2]) ./ (SummedActs[j + (i-1)*Con1*Con2]);
    probs[j + (i-1)*Con1*Con2,3] = (R[3] * acts_DIP[j + (i-1)*Con1*Con2]) ./ (SummedActs[j + (i-1)*Con1*Con2]);
    probs[j + (i-1)*Con1*Con2,4] = (R[4] * acts_DIOP[j + (i-1)*Con1*Con2]) ./ (SummedActs[j + (i-1)*Con1*Con2]);
    probs[j + (i-1)*Con1*Con2,5] = (R[5] * acts_NPL[j + (i-1)*Con1*Con2]) ./ (SummedActs[j + (i-1)*Con1*Con2]);
      }
  }

}

model {
  
  // priors for hyper parameters
  hyper_pars[1] ~ normal(20,10); // c
  hyper_pars[2] ~ normal(2,10); // a
  hyper_pars[3] ~ normal(1,10); // r
  hyper_pars[4] ~ normal(1,10); // EU
  hyper_pars[5] ~ normal(0,10); // d

  
  Lcorr ~ lkj_corr_cholesky(2);
  sigma ~ gamma(1,0.01);
  
  
  // Loop over subjects
  subj_pars[,] ~ multi_normal_cholesky(hyper_pars,Lcorr);
  
  
  // Loop over subjects
  
  for(i in 1:N){
    // Draw subject parameters from truncated normal
  
    
    for (j in 1:Con1*Con2) {
    
      // draw data from probabilities determined by MMM parms
     
    
      count[j + (i-1)*Con1*Con2,]  ~ multinomial(probs[j + (i-1)*Con1*Con2,]);  
      
    }
  }
}

generated quantities 
{
matrix[J,J] Omega;
matrix[J,J] Sigma;

Omega = multiply_lower_tri_self_transpose(Lcorr);
Sigma = quad_form_diag(Omega,sigma);
  
}

