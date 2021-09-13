library(brms)
library(tidyverse)
library(bayesplot)
library(progress)
library(tictoc)
library(HDInterval)
library(MCMCglmm)
library(rstan)
library(psych)

# start fresh
rm(list=ls())   # clean up workspace

# set working directory to source file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# set stan options
options(mc.cores = parallel::detectCores())
options(buildtools.check = function(action) TRUE)

# Load Functions

source("M3_functions.R")

# Specify Simulation Factors----

N <- c(5)
K <- c(8,16,32,64)
nRetrievals <- 100
minFT <- 0.5
maxFT <- 2
nFT <- c(2,4) # 2,4,10 Conditions between 0.2 and 2
SampleSize <- 100
reps2con <- 100
n_sim <- length(N)*length(K)*length(nRetrievals)*length(nFT)*reps2con
n_obs <- 200*200 + 200*400



# Set Range for Parameter Means
range_muC <- c(1,100)
range_muA <- c(0,0.5)
#range_muF <- c(0,1) # fix to 0.5
range_muE <-c(0,0.5)
range_muR <- c(0,25) # range 0 - 25 empirical derived
eta <- 5 # Simulated N = 10000 with eta = 5, 95 % of all values lie within 0 -+ 0.56


sigC <- c(0.125,0.5)
sigA <- c(0.125,0.5)
sigF <- c(0.0001,0.1)
sigE <- c(1,2)
sigR <- c(0.125,0.5) # abhänbgig von removal parameter -> analog zu c und a
sigB <- c(0.0001, 0.1)

# Test Cons (500 Rets, 16 NPLs, 3 IOPs, 4 FT Cons)
i <- 3
j <- 5
k <- 4
l <- 1
n <- 2


# Prepare Simulation 
pb <- progress_bar$new(total = n_sim,
                       format = "(:spin) [:bar] :percent :eta")

m <- 1

# Preallocate Data Objects for iterative data ----
benchmarks_subject_pars<- matrix(NaN,nrow=n_obs, ncol = 75)
colnames(benchmarks_subject_pars) <- c("N","K","Retrievals","n_conFT","FT","meanPC","minPC","maxPC","hyper_mu_c","hyper_sig_c",
                                       "hyper_mu_a","hyper_sig_a","relCA","hyper_mu_f","sig_f","hyper_mu_e","hyper_sig_e","hyper_mu_r","hyper_sig_r",
                                       
                                       "emp_cor_ca","est_cor_ca","lower_95_HDI_cor_ca","upper_95_HDI_cor_ca",
                                       "emp_cor_cf","est_cor_cf","lower_95_HDI_cor_cf","upper_95_HDI_cor_cf",
                                       "emp_cor_af", "est_cor_af","lower_95_HDI_cor_af","upper_95_HDI_cor_af",
                                       "emp_cor_ce", "est_cor_ce","lower_95_HDI_cor_ce","upper_95_HDI_cor_ce",
                                       "emp_cor_cr", "est_cor_cr","lower_95_HDI_cor_cr","upper_95_HDI_cor_cr",
                                       
                                       "emp_cor_ae", "est_cor_ae","lower_95_HDI_cor_ae","upper_95_HDI_cor_ae",
                                       "emp_cor_ar", "est_cor_ar","lower_95_HDI_cor_ar","upper_95_HDI_cor_ar",
                                       
                                       "emp_cor_fe", "est_cor_fe","lower_95_HDI_cor_fe","upper_95_HDI_cor_fe",
                                       "emp_cor_fr", "est_cor_fr","lower_95_HDI_cor_fr","upper_95_HDI_cor_fr",
                                       
                                       
                                       "cor_subj_muA","cor_subj_modeA","cor_subj_muC","cor_subj_modeC",
                                       "cor_subj_muF","cor_subj_modeF","cor_subj_muE","cor_subj_modeE","cor_subj_muR","cor_subj_modeR",
                                       
                                       
                                       "mu_subj_resA","sig_subj_resA","mu_subj_resC","sig_subj_resC", 
                                       "mu_subj_resF", "sig_subj_resF", "mu_subj_resE","sig_subj_resE","mu_subj_resR","sig_subj_resR")

benchmarks_hyper_pars <- matrix(NaN,nrow= n_obs/reps2con, ncol=23)
colnames(benchmarks_hyper_pars) <- c("N","K","Retrievals","n_conFT","FT","meanPC","minPC","maxPC","hyper_cor_a", 
                                     "hyper_cor_c","hyper_cor_f","hyper_cor_e","hyper_cor_r",
                                     "mu_res_hyper_a","sigma_res_hyper_a","mu_res_hyper_c",
                                     "sigma_res_hyper_c","mu_res_hyper_f","sigma_res_hyper_f",
                                     "mu_res_hyper_e","sigma_res_hyper_e","mu_res_hyper_r","sigma_res_hyper_r")


hyper_pars <- matrix(NaN,  nrow = n_obs, ncol=26)
colnames(hyper_pars) <-c("Repetition","N","K","Retrievals","n_conFT","FT","meanPC","minPC","maxPC",
                         "real_hyper_mu_a","est_hyper_mu_a","rhat_hyper_mu_a","real_hyper_mu_c",
                         "est_hyper_mu_c","rhat_hyper_mu_c", "real_hyper_mu_f_normal", "est_hyper_mu_f_normal", 
                         "real_hyper_log_mu_f","est_hyper_log_mu_f","rhat_hyper_mu_f",
                         "real_hyper_mu_e","est_hyper_mu_e","rhat_hyper_mu_e",
                         "real_hyper_mu_r","est_hyper_mu_r","rhat_hyper_mu_r")

tic()
# Simulation ----

for(i in 1:length(N)) 
{
  for (j in 1:length(K)) 
  {
      for(l in 1:length(nFT))
      {
        for (n in 1:reps2con)  
        {
          
          # Set Seed for pseudo randomized processes
          time <- Sys.time()
          seed <- as.integer(Sys.time()) %% 100000
          
          # set seed for random number generator
          set.seed(seed)
          
          # Set Up Data Object
          
          simulation.ss <- list()
          
          
          # Save Simulation Conditions----
          conN <- N[i]
          conK <- K[j]
          conRet <- 100
          con_nFT <- nFT[l]
          conFT <- seq(from = minFT, to = maxFT, length.out = con_nFT) # eventuell log scale 0.2 0.8 2.4 oder so?
          
          
          # Sample Hyper-Parameter Means with C as fixpoint ----
          relCA <- runif(1, min = range_muA[1],max = range_muA[2])
          Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
          Mean_Apar <- Mean_Cpar*(relCA)
          Mean_Epar <- runif(1, min =range_muE[1], max = range_muE[2])
          Mean_Rpar <- runif(1, min= range_muR[1], max = range_muR[2])
          Mean_bpar <- 0.1
          Mean_Fpar <- 0.50
          log_mu_f <- log(Mean_Fpar/(1-Mean_Fpar))
          
          hyper_mus <- c(Mean_Cpar,Mean_Apar,log_mu_f,Mean_Epar,Mean_Rpar, Mean_bpar)
          
          
          # Sample Variances and Set Covariances----
          
          sig_c <- runif(1, min = sigC[1], max = sigC[2])*Mean_Cpar
          sig_a <- runif(1, min = sigA[1], max = sigA[2])*Mean_Apar 
          sig_e <- runif(1, min = sigE[1], max= sigE[2])
          sig_r <- runif(1, min = sigR[1], max= sigR[2])*Mean_Rpar
          sig_b <- 0.001
          sig_f <- 0.001
          
          sigs <-c(sig_c,sig_a,sig_f,sig_e,sig_r,sig_b)
          Sig <- diag(length(hyper_mus))
          
          Sig[1,1] <- (sig_c)^2
          Sig[2,2] <- (sig_a)^2
          Sig[3,3] <- (sig_f)^2
          Sig[4,4] <- (sig_e)^2
          Sig[5,5] <- (sig_r)^2
          Sig[6,6] <- (sig_b)^2
          
          
          # Set Correlations for Parameters ----
          
          # Sample Covariance Matrix Sigma
          
          omega <- rlkjcorr(1,length(hyper_mus),eta)
          
          # Little Hack for fixing coveriance of b to zer0
          
          #omega[6,1:5] = omega[1:5,6] = 0 #fix cov of b to 0
          #omega[3,1:5] = omega[1:5,3] = 0 #fix cov of f to 0
          
          
          Sigma <- cor2cov(omega,sigs)
          
          
          # Sample Parameters from MVN ----
          
          parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=Sigma,
                                      lower=c(0,0,-Inf,0,0,0),upper = c(Inf,Inf,Inf,Inf,Inf,Inf))
          # Merge Parameters to one Matrix
          colnames(parms) <- c("conA","genA","f","e","r","baseA")
          parms[,6] <- 0.1
          parms[,3] <- 0

          emp_cor_ca <-round(cor(parms[,1],parms[,2]),3)
          #emp_cor_cf <-round(cor(parms[,1],parms[,3]),3)
          #emp_cor_af <-round(cor(parms[,2],parms[,3]),3)  
          emp_cor_ce <-round(cor(parms[,1],parms[,4]),3)
          emp_cor_cr <-round(cor(parms[,1],parms[,5]),3)
          emp_cor_ae <-round(cor(parms[,2],parms[,4]),3)  
          emp_cor_ar <-round(cor(parms[,2],parms[,5]),3)  
          #emp_cor_fe <-round(cor(parms[,3],parms[,4]),3)  
          #emp_cor_fr <-round(cor(parms[,3],parms[,5]),3)  
          
          
          parms[,3] <- 1 / (1+exp(-parms[,3]))
          
          
          
          # Simulate Data for Estimation ----
          
          ParmsFT <- matrix(rep(parms,each =con_nFT), nrow = length(parms[,1])*con_nFT, ncol = ncol(parms), byrow = F)
          colnames(ParmsFT) <- c("conA","genA","f","e","r","baseA")
          FT <- rep(conFT,length.out = nrow(ParmsFT))
          
          data <- simData_CSpanEE(ParmsFT,as.vector(respOpt_Cspan(conN,conK)), conRet,FT)
          
          # Calculate Percentage Correct per Condition
          
          
          # Generate Stan Data ----
          
          
          stan.dat <- list(count = data[,4:8], 
                           K = 5,
                           R = as.vector(respOpt_Cspan(conN,conK)),
                           J = length(sigs)-2,
                           N = length(unique(data[,"ID"])),
                           Con = length(unique(data[,"Freetime"])),
                           Freetime = unique(data[,"Freetime"]),
                           scale_b = 0.1,
                           fixed_f = .5)
          
          # Fit-The-Shit----
tic()
                    fit_M3 <- stan(file = "M3_ComplexSpan_EE_LKJ_Cholesky_BF.stan",data=stan.dat,
                         warmup = 2000, iter = 4000,
                         chains = 4,refresh = 100, init=init_fun, control=list(metric="unit_e"))
toc()
          # Extract Posterior Parameters from fit object---- 
          
          post_samples <- rstan::extract(fit_M3, pars=c("subj_pars","Omega"), inc_warmup = F)
          
          # Hyper Parameter
          hyper_means <- as.vector(get_posterior_mean(fit_M3, par="hyper_pars")[,5])
          
          mean_est_hyper_c <- hyper_means[1]
          mean_est_hyper_a <- hyper_means[2]
          # mean_est_hyper_f_log <- hyper_means[3]
          # mean_est_hyper_f_normal <- (1 / (1+exp(-hyper_means[3])))
          mean_est_hyper_e <- hyper_means[3]
          mean_est_hyper_r <- hyper_means[4]
          
          # Subject Parameters
          
          means_c <- colMeans(post_samples$subj_pars[,,1])
          means_a <- colMeans(post_samples$subj_pars[,,2])
          #means_f <- colMeans(post_samples$f)
          means_e <- colMeans(post_samples$subj_pars[,,3])
          means_r <- colMeans(post_samples$subj_pars[,,4])
          
          # cor(means_c, parms[,1])
          # cor(means_a, parms[,2])      
          # cor(means_e, parms[,4])      
          # cor(means_r, parms[,5])      
         
          
          
          # Correlations
          
          est_cor_ca <- rstan::get_posterior_mean(fit_M3, par=c("Omega[1,2]"))[5]
          #est_cor_cf <- rstan::get_posterior_mean(fit_M3, par=c("Omega[1,3]"))[5]
         # est_cor_af  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[2,3]"))[5]
          est_cor_ce  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[1,4]"))[5]
          est_cor_cr  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[1,5]"))[5]
          est_cor_ae  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[2,4]"))[5]
          est_cor_ar  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[2,5]"))[5]
         # est_cor_fe  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[3,4]"))[5]
          #est_cor_fr  <- rstan::get_posterior_mean(fit_M3, par=c("Omega[3,5]"))[5]
          
          
          
          # Calculate Highest Density Interval----
          
          HDI_c <- t(round(hdi(post_samples$subj_pars[,,1], credMass = .95),3))
          HDI_a <- t(round(hdi(post_samples$subj_pars[,,2], credMass = .95),3))
         # HDI_f <- t(round(hdi(post_samples$f, credMass = .95),3))
          HDI_e <- t(round(hdi(post_samples$subj_pars[,,4], credMass = .95),3))
          HDI_r <- t(round(hdi(post_samples$subj_pars[,,5], credMass = .95),3))
          
          HDI_cor_ca <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[1,2]"),inc_warmup=F),credMass= 0.95))),2)
         # HDI_cor_cf <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[1,3]"),inc_warmup=F),credMass= 0.95))),2)
         # HDI_cor_af <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[2,3]"),inc_warmup=F),credMass= 0.95))),2)
          HDI_cor_ce <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[1,4]"),inc_warmup=F),credMass= 0.95))),2)
          HDI_cor_cr <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[1,5]"),inc_warmup=F),credMass= 0.95))),2)
          HDI_cor_ae <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[2,4]"),inc_warmup=F),credMass= 0.95))),2)
          HDI_cor_ar <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[2,5]"),inc_warmup=F),credMass= 0.95))),2)
          #HDI_cor_fe <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[3,4]"),inc_warmup=F),credMass= 0.95))),2)
          #HDI_cor_fr <- round(as.vector(unlist(hdi(extract(fit_M3, pars=c("Omega[3,5]"),inc_warmup=F),credMass= 0.95))),2)
          
          
          HDI_include_c <- c()
          HDI_include_a <- c()
          #HDI_include_f <- c()
          HDI_include_e <- c()
          HDI_include_r <- c()
          
          HDI_include_cor_ca <- isTRUE(between(est_cor_ca, HDI_cor_ca[1], HDI_cor_ca[2]))
         # HDI_include_cor_cf <- isTRUE(between(est_cor_cf, HDI_cor_ca[1], HDI_cor_ca[2]))
         # HDI_include_cor_af <- isTRUE(between(est_cor_af, HDI_cor_ca[1], HDI_cor_ca[2]))
          HDI_include_cor_ce <- isTRUE(between(est_cor_ce, HDI_cor_ce[1], HDI_cor_ce[2]))
          HDI_include_cor_cr <- isTRUE(between(est_cor_cr, HDI_cor_cr[1], HDI_cor_cr[2]))
          HDI_include_cor_ae <- isTRUE(between(est_cor_ae, HDI_cor_ae[1], HDI_cor_ae[2]))
          HDI_include_cor_ar <- isTRUE(between(est_cor_ar, HDI_cor_ar[1], HDI_cor_ar[2]))
         # HDI_include_cor_fe <- isTRUE(between(est_cor_fe, HDI_cor_fe[1], HDI_cor_fe[2]))
          #HDI_include_cor_fr <- isTRUE(between(est_cor_fr, HDI_cor_fr[1], HDI_cor_fr[2]))
          
          
          for (f in 1:SampleSize)
          {
            HDI_include_c[f]<-isTRUE(between(means_c[f], HDI_c[f,1], HDI_c[f,2]))
          }
          
          
          for (f in 1:SampleSize)
          {
            HDI_include_a[f]<-isTRUE(between(means_a[f], HDI_a[f,1], HDI_a[f,2]))
            
          }
          
          for (f in 1:SampleSize)
          {
            HDI_include_f[f]<-isTRUE(between(means_f[f], HDI_f[f,1], HDI_f[f,2]))
            
          }
          
          for (f in 1:SampleSize)
          {
            HDI_include_e[f]<-isTRUE(between(means_e[f], HDI_e[f,1], HDI_e[f,2]))
            
          }
          
          for (f in 1:SampleSize)
          {
            HDI_include_r[f]<-isTRUE(between(means_r[f], HDI_r[f,1], HDI_r[f,2]))
            
          }
          
          # Calculate Posterior Mode ----
          
          mode_c <- posterior.mode(post_samples$subj_pars[,,1])
          mode_a <- posterior.mode(post_samples$subj_pars[,,2])
         # mode_f <- posterior.mode(post_samples$f)
          mode_e <- posterior.mode(post_samples$subj_pars[,,4])
          mode_r <- posterior.mode(post_samples$subj_pars[,,5])
          
          # Extract rhat ----
          
          rhats_subj_pars <-as.data.frame(t(rhat(fit_M3, pars=c("subj_pars"))))
          rhats_f <-as.data.frame(t(rhat(fit_M3, pars=c("f"))))
          
          rhats_subj_pars <- rhats_subj_pars %>% gather(value = rhat) %>%
            separate(key, sep=",", into = c("Subj","par")) %>% 
            mutate(par = as.factor(par), par = recode(par,"1]" = "c","2]"="a", "4]" = "e", "5]" = "r" )) %>%
            pivot_wider(names_from = par, values_from=rhat)
          
          rhats_f<- rhats_f %>% gather(value = rhat)
          
          rhats_a<- rhats_subj_pars$a
          rhat_max_a <- max(rhats_a)
          
          rhats_c<- rhats_subj_pars$c
          rhat_max_c <- max(rhats_c)
          
          rhats_f<- rhats_f$rhat
          rhat_max_f <- max(rhats_f)
          
          rhats_e <- rhats_subj_pars$e
          rhat_max_e <- max(rhats_e)
          
          rhats_r <- rhats_subj_pars$r
          rhat_max_r <- max(rhats_r)
          
          
          
          rhat_cor_ca <- as.vector(rhat(fit_M3, pars=c("Omega[1,2]")))
          rhat_cor_cf <- as.vector(rhat(fit_M3, pars=c("Omega[1,3]")))
          rhat_cor_af <- as.vector(rhat(fit_M3, pars=c("Omega[2,3]")))
          rhat_cor_ce <- as.vector(rhat(fit_M3, pars=c("Omega[1,4]")))
          rhat_cor_cr <- as.vector(rhat(fit_M3, pars=c("Omega[1,5]")))
          rhat_cor_ae <- as.vector(rhat(fit_M3, pars=c("Omega[2,4]")))
          rhat_cor_ar <- as.vector(rhat(fit_M3, pars=c("Omega[2,5]")))
          rhat_cor_fe <- as.vector(rhat(fit_M3, pars=c("Omega[3,4]")))
          rhat_cor_fr <- as.vector(rhat(fit_M3, pars=c("Omega[3,5]")))
          
          rhat_hyper_c <- as.vector(rhat(fit_M3, pars=c("hyper_pars")))[1]
          rhat_hyper_a <- as.vector(rhat(fit_M3, pars=c("hyper_pars")))[2]
          rhat_hyper_f <- as.vector(rhat(fit_M3, pars=c("hyper_pars")))[3]
          rhat_hyper_e <- as.vector(rhat(fit_M3, pars=c("hyper_pars")))[4]
          rhat_hyper_r <- as.vector(rhat(fit_M3, pars=c("hyper_pars")))[5]
          
          
          
          # Sampling Correlation between Chains
          
          mcmc_cor <- round(fisherz2r(mean(fisherz(apply(post_samples$subj_pars,MARGIN=2, FUN = "cor")[2,]))),2)
          
          
          # Calculate Residuals ----
          
          resC <- means_c - parms[,1]
          resA <- means_a - parms[,2]
          resF <- means_f - parms[,3]
          resE <- means_e - parms[,4]
          resR <- means_r - parms[,5]
          
     
          
          # Save Data in Object ----
          
          simulation <- data.frame("ID" = data[,1],
                                   "N"= conN,
                                   "K"= conK,
                                   "nFT" = con_nFT,
                                   "FT" = FT,
                                   "Retrievals" = conRet,
                                   "Repetition" = n,
                                   "PC" = (data[,3]),
                                   "RelCa" = relCA,
                                   "mcmc_cor" = mcmc_cor,
                                   
                                   "real_cor_ca"=emp_cor_ca,
                                   "est_cor_ca" = est_cor_ca,
                                   "lower_95_HDI_cor_ca" = HDI_cor_ca[1], 
                                   "upper_95_HDI_cor_ca" = HDI_cor_ca[2], 
                                   "rhat_cor_ca" = rhat_cor_ca,
                                   
                                   "real_cor_cf"=emp_cor_cf,
                                   "est_cor_cf" = est_cor_cf,
                                   "lower_95_HDI_cor_cf" = HDI_cor_cf[1], 
                                   "upper_95_HDI_cor_cf" = HDI_cor_cf[2], 
                                   "rhat_cor_cf" = rhat_cor_cf,
                                   
                                   "real_cor_af"=emp_cor_af,
                                   "est_cor_af" = est_cor_af,
                                   "lower_95_HDI_cor_af" = HDI_cor_af[1], 
                                   "upper_95_HDI_cor_af" = HDI_cor_af[2], 
                                   "rhat_cor_af" = rhat_cor_af,
                                   
                                   "real_cor_ce"=emp_cor_ce,
                                   "est_cor_ce" = est_cor_ce,
                                   "lower_95_HDI_cor_ce" = HDI_cor_ce[1], 
                                   "upper_95_HDI_cor_ce" = HDI_cor_ce[2], 
                                   "rhat_cor_ce" = rhat_cor_ce,
                                   
                                   "real_cor_cr"=emp_cor_cr,
                                   "est_cor_cr" = est_cor_cr,
                                   "lower_95_HDI_cor_cr" = HDI_cor_cr[1], 
                                   "upper_95_HDI_cor_cr" = HDI_cor_cr[2], 
                                   "rhat_cor_cr" = rhat_cor_cr,
                                   
                                   "real_cor_ae"=emp_cor_ae,
                                   "est_cor_ae" = est_cor_ae,
                                   "lower_95_HDI_cor_ae" = HDI_cor_ae[1], 
                                   "upper_95_HDI_cor_ae" = HDI_cor_ae[2], 
                                   "rhat_cor_ae" = rhat_cor_ae,
                                   
                                   "real_cor_ar"=emp_cor_ar,
                                   "est_cor_ar" = est_cor_ar,
                                   "lower_95_HDI_cor_ar" = HDI_cor_ar[1], 
                                   "upper_95_HDI_cor_ar" = HDI_cor_ar[2], 
                                   "rhat_cor_ar" = rhat_cor_ar,
                                   
                                   "real_cor_fe"=emp_cor_fe,
                                   "est_cor_fe" = est_cor_fe,
                                   "lower_95_HDI_cor_fe" = HDI_cor_fe[1], 
                                   "upper_95_HDI_cor_fe" = HDI_cor_fe[2], 
                                   "rhat_cor_fe" = rhat_cor_fe,
                                   
                                   "real_cor_fr"=emp_cor_fr,
                                   "est_cor_fr" = est_cor_fr,
                                   "lower_95_HDI_cor_fr" = HDI_cor_fr[1], 
                                   "upper_95_HDI_cor_fr" = HDI_cor_fr[2], 
                                   "rhat_cor_fr" = rhat_cor_fr,
                                   
                                   "hyper_mu_c"= Mean_Cpar,
                                   "hyper_sig_c" = sig_c,
                                   "hyper_est_c" = mean_est_hyper_c,
                                   "rhat_hyper_mu_c" = rhat_hyper_c,
                                   
                                   "hyper_mu_a"= Mean_Apar,
                                   "hyper_sig_a" = sig_a,
                                   "hyper_est_a" = mean_est_hyper_a,
                                   "rhat_hyper_mu_a" = rhat_hyper_a,
                                   
                                   "hyper_mu_f" = Mean_Fpar,
                                   "hyper_sig_f" = sig_f,
                                   "hyper_est_f_normal" = mean_est_hyper_f_normal,
                                   "rhat_hyper_mu_f" = rhat_hyper_f,
                                   
                                   "hyper_mu_e"= Mean_Epar,
                                   "hyper_sig_e" = sig_e,
                                   "hyper_est_e" = mean_est_hyper_e,
                                   "rhat_hyper_mu_e" = rhat_hyper_e,
                                   
                                   "hyper_mu_r"= Mean_Rpar,
                                   "hyper_sig_r" = sig_r,
                                   "hyper_est_r" = mean_est_hyper_r,
                                   "rhat_hyper_mu_r" = rhat_hyper_r,
                                   
                                   "mu_est_c" = means_c,
                                   "mode_est_c" = mode_c,
                                   "lower_95_HDI_c" = HDI_c[,1],
                                   "upper_95_HDI_c" = HDI_c[,2],
                                   "HDI_include_est_c" = HDI_include_c, 
                                   "resC" = resC,
                                   "max_rhat_c" = rhat_max_c,
                                 
                                   "mu_est_a" = means_a,
                                   "mode_est_a" = mode_a,
                                   "lower_95_HDI_a" = HDI_a[,1],
                                   "upper_95_HDI_a" = HDI_a[,2],
                                   "HDI_include_est_a" = HDI_include_a, 
                                   "resA"  = resA,
                                   "max_rhat_a" = rhat_max_a,
                                   
                                   "mu_est_f" = means_f,
                                   "mode_est_f" = mode_f,
                                   "lower_95_HDI_f" = HDI_f[,1],
                                   "upper_95_HDI_f" = HDI_f[,2],
                                   "HDI_include_est_f" = HDI_include_f, 
                                   "resF" = resF,
                                   "max_rhat_f" = rhat_max_f,
                                   
                                   "mu_est_e" = means_e,
                                   "mode_est_e" = mode_e,
                                   "lower_95_HDI_e" = HDI_e[,1],
                                   "upper_95_HDI_e" = HDI_e[,2],
                                   "HDI_include_est_e" = HDI_include_e, 
                                   "resE" = resE,
                                   "max_rhat_e" = rhat_max_e,
                                   
                                   
                                   "mu_est_r" = means_r,
                                   "mode_est_r" = mode_r,
                                   "lower_95_HDI_r" = HDI_r[,1],
                                   "upper_95_HDI_r" = HDI_r[,2],
                                   "HDI_include_est_r" = HDI_include_r, 
                                   "resR" = resR,
                                   "max_rhat_r" = rhat_max_r,
                                   
                                   
                                   "mu_real_c" = parms[,1],
                                   "mu_real_a" = parms[,2],
                                   "mu_real_f" = parms[,3],
                                   "mu_real_e" = parms[,4],
                                   "mu_real_r" = parms[,5])
          
          
          
          if(m == 1){
            start_row <- 1
            end_row <- SampleSize*con_nFT 
          }else{
            start_row <- end_row + 1
            end_row <- (start_row - 1) + SampleSize*con_nFT
          }
          
          benchmarks_subject_pars[start_row:end_row,] <-cbind(conN, conK,conRet,con_nFT, conFT,
                                                              meanPC = tapply(data[,3],data[,2],  "mean"),
                                                              minPC = tapply(data[,3],data[,2],  "min"),
                                                              maxPC = tapply(data[,3],data[,2],  "max"),
                                                              Mean_Cpar,sig_c,Mean_Apar,sig_a,relCA,
                                                              Mean_Fpar,sig_f, Mean_Epar, sig_e, Mean_Rpar, sig_r,
                                                              emp_cor_ca,est_cor_ca,
                                                              HDI_cor_ca[1],HDI_cor_ca[2],
                                                              emp_cor_cf, est_cor_cf,
                                                              HDI_cor_cf[1],HDI_cor_cf[2],
                                                              emp_cor_af, est_cor_af, 
                                                              HDI_cor_af[1],HDI_cor_af[2],
                                                              emp_cor_ce, est_cor_ce,
                                                              HDI_cor_ce[1],HDI_cor_ce[2],
                                                              emp_cor_cr, est_cor_cr,
                                                              HDI_cor_cr[1],HDI_cor_cr[2],
                                                              emp_cor_ae, est_cor_ae,
                                                              HDI_cor_ae[1],HDI_cor_ae[2],
                                                              emp_cor_ar, est_cor_ar,
                                                              HDI_cor_ar[1],HDI_cor_ar[2],
                                                              emp_cor_fe, est_cor_fe,
                                                              HDI_cor_fe[1],HDI_cor_fe[2],
                                                              emp_cor_fr, est_cor_fr,
                                                              HDI_cor_fr[1],HDI_cor_fr[2],
                                                              cor(simulation$mu_est_a,simulation$mu_real_a),
                                                              cor(simulation$mode_est_a, simulation$mu_real_a),
                                                              cor(simulation$mu_est_c,simulation$mu_real_c),
                                                              cor(simulation$mode_est_c, simulation$mu_real_c),
                                                              cor(simulation$mu_est_f,simulation$mu_real_f),
                                                              cor(simulation$mode_est_f, simulation$mu_real_f),
                                                              
                                                              cor(simulation$mu_est_e,simulation$mu_real_e),
                                                              cor(simulation$mode_est_e, simulation$mu_real_e),
                                                              cor(simulation$mu_est_r,simulation$mu_real_r),
                                                              cor(simulation$mode_est_r, simulation$mu_real_r),
                                                              
                                                              mean(simulation$mu_est_a-simulation$mu_real_a), 
                                                              sd(simulation$mu_est_a-simulation$mu_real_a),  
                                                              mean(simulation$mu_est_c-simulation$mu_real_c),
                                                              sd(simulation$mu_est_c-simulation$mu_real_c),
                                                              mean(simulation$mu_est_f-simulation$mu_real_f),
                                                              sd(simulation$mu_est_f-simulation$mu_real_f),
                                                              
                                                              mean(simulation$mu_est_e-simulation$mu_real_e),
                                                              sd(simulation$mu_est_e-simulation$mu_real_e),
                                                              mean(simulation$mu_est_r-simulation$mu_real_r),
                                                              sd(simulation$mu_est_r-simulation$mu_real_r)) 

          hyper_pars[start_row:end_row,] <- cbind(n,conN,conK,conRet,con_nFT,conFT,
                                                  meanPC = tapply(data[,3],data[,2],  "mean"),
                                                  minPC = tapply(data[,3],data[,2],  "min"),
                                                  maxPC = tapply(data[,3],data[,2],  "max"),
                                                  Mean_Apar, 
                                                  mean_est_hyper_a, 
                                                  rhat_hyper_a,
                                                  Mean_Cpar,
                                                  mean_est_hyper_c,
                                                  rhat_hyper_c,
                                                  Mean_Fpar,
                                                  mean_est_hyper_f_normal,
                                                  log_mu_f,
                                                  mean_est_hyper_f_log,
                                                  rhat_hyper_f,
                                                  Mean_Epar, 
                                                  mean_est_hyper_e, 
                                                  rhat_hyper_e,
                                                  Mean_Rpar, 
                                                  mean_est_hyper_r, 
                                                  rhat_hyper_r)
          
          
          
          simulation.ss[[paste0("N_",conN,"_K_",conK,"_Retrievals_",conRet,"n_conFT_",con_nFT,"_Set_", n)]] <- simulation
          simulation.ss[[paste0("Benchmark Subject Pars")]] <- benchmarks_subject_pars[start_row:end_row,]
          simulation.ss[[paste0("Data")]] <- data
          saveRDS(simulation.ss, file=paste0("DataTest/ComplexSpan_LKJ_EE_","N_",conN,"_K_",conK,"_Retrievals_",conRet, "n_conFT_",con_nFT,"_Set_", n,".RDS"))
          
          
          
          if(m == n_sim)
            
          {
            hyper_pars <- as.data.frame(hyper_pars)
            
            benchmarks_hyper_pars <- hyper_pars %>% 
              group_by(N,K,Retrievals,n_conFT,FT) %>% 
              summarise(meanPC = mean(meanPC),minPC = mean(minPC), maxPC = mean(maxPC),
                        cor_hyper_muC = cor(real_hyper_mu_c,est_hyper_mu_c),
                        cor_hyper_muA = cor(real_hyper_mu_a,est_hyper_mu_a),
                        cor_hyper_muF_normal = cor(real_hyper_mu_f_normal,est_hyper_mu_f_normal),
                        cor_hyper_muE = cor(real_hyper_mu_e,est_hyper_mu_e),
                        cor_hyper_muR = cor(real_hyper_mu_r,est_hyper_mu_r),
                        mean_hyper_resC = mean(est_hyper_mu_c - real_hyper_mu_c),
                        sig_hyper_resC = sd(est_hyper_mu_c - real_hyper_mu_c),
                        mean_hyper_resA = mean(est_hyper_mu_a - real_hyper_mu_a),
                        sig_hyper_resA = sd(est_hyper_mu_a - real_hyper_mu_a),
                        mean_hyper_resF = mean(est_hyper_mu_f_normal - real_hyper_mu_f_normal),
                        sig_hyper_resF = sd(est_hyper_mu_f_normal - real_hyper_mu_f_normal),
                        mean_hyper_resE = mean(est_hyper_mu_e - real_hyper_mu_e),
                        sig_hyper_resE = sd(est_hyper_mu_e - real_hyper_mu_e),
                        mean_hyper_resR = mean(est_hyper_mu_r - real_hyper_mu_r),
                        sig_hyper_resR = sd(est_hyper_mu_r - real_hyper_mu_r))


            
            saveRDS(benchmarks_subject_pars, file="DataTest/ComplexSpan_LKJ_EE_Benchmarks_subject_pars.RDS")
            saveRDS(benchmarks_hyper_pars, file="DataTest/ComplexSpan_LKJ_EE_Benchmarks_hyper_pars.RDS")
            saveRDS(hyper_pars, file ="DataTest/ComplexSpan_LKJ_EE_hyper_pars.RDS")
          }
          
          
          m <- m + 1
          
          gc(full = T)
          
          pb$tick()
          
        }
        
      }
      
    }
    
  }

