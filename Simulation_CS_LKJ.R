#!/usr/bin/env Rscript

library(brms)
library(rstan)
library(tidyverse)
library(bayesplot)
library(progress)
library(tictoc)
library(HDInterval)
library(MCMCglmm)
library(optparse) # to create a command line interface

library(sys) # to infer system environment variables
library(compiler) # to accelerate using a just in time compiler

# start fresh
rm(list=ls())   # clean up workspace

#CM: not needed here, the working directory needs to be
#    defined on the command line or will be the submit dir.
# set working directory to source file location
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set stan options
options(mc.cores = parallel::detectCores())

# we can create a default output directory path, here:
default.out = paste("/home/", Sys.getenv("USER"), sep="", collapse=NULL)

option_list <- list(
    make_option(c("-o", "--outdir"), action="store",
                  default=default.out), 
    make_option(c("-N"), type="integer",
                  help="helpmessage goes here",
                  metavar="number"),
    make_option(c("-K"), type="integer",
                  help="helpmessage goes here",
                  metavar="number"),
    make_option(c("-n", "--nRetrievals"), type="integer",
                  help="helpmessage goes here",
                  metavar="number"),

                    )


# Load functions

source("M3_functions.R")

# Specify Simulation Factors----

N <- seq(1:5) 
K <- c(1,2,4,8,20)
nRetrievals <- c(100,250,500,1000)
SampleSize <- 100
reps2con <- 1 #100 # times to fit the conditions
n_sim <- length(N)*length(K)*length(nRetrievals)*reps2con

# Test Cons
i <- 4
j <- 5
k <- 4
n <- 2

# draw Hyper-Parameter Means and variances from uniform distributions ----

range_muC <- c(1,100)
range_muA <- c(0,0.5)
range_muF <- c(0,1)

sigC <- c(0.5,5)
sigA <- c(0.2,5)
sigF <- c(0.1,0.3)
sigB <- c(0.0001, 0.1)

# Prepare Simulation 
pb <- progress_bar$new(total = length(N)*length(K)*length(nRetrievals)*reps2con,
                       format = "(:spin) [:bar] :percent :eta")

m <- 1

# Preallocate Data Objects for iterative data ----
benchmarks_subject_pars<- matrix(NaN,nrow=n_sim, ncol = 23)
colnames(benchmarks_subject_pars) <- c("N","K","Retrievals","Correlation","meanPC","minPC","maxPC","hyper_mu_c","hyper_sig_c",
                                       "hyper_mu_a","hyper_sig_a","relCA","hyper_mu_f","sig_f", "cor_subj_muA","cor_subj_muC","cor_subj_muF",
                                       "mu_subj_resA","sig_subj_resA","mu_subj_resC","sig_subj_resC", "mu_subj_resF", "sig_subj_resF")

benchmarks_hyper_pars <- matrix(NaN,nrow= n_sim/reps2con, ncol=10)
colnames(benchmarks_hyper_pars) <- c("Repetition","hyper_cor_a", "hyper_cor_c","hyper_cor_f","mu_res_hyper_a",
                                     "sigma_res_hyper_a","mu_res_hyper_c","sigma_res_hyper_c","mu_res_hyper_f",
                                     "sigma_res_hyper_f")


hyper_pars <- matrix(NaN,  nrow = n_sim, ncol=13)
colnames(hyper_pars) <-c("Repetition","N","K","Retrievals","meanPC","minPc","maxPC",
                         "real_hyper_mu_a","est_hyper_mu_a","real_hyper_mu_c",
                         "est_hyper_mu_c", "real_hyper_mu_f", "est_hyper_mu_f")



tic()
# Simulation ----

for(i in 1:length(N)) 
{
  for (j in 1:length(K)) 
  {
    for(k in 1:length(nRetrievals)) 
    {
      for (n in 1:reps2con)  
      {
        
        # Set Up Data Object
        
        simulation.ss <- list()
        
        
        # Save Simulation Conditions----
        conN <- N[i]
        conK <- K[j]
        conRet <- nRetrievals[k]
        par_cor <- runif(1, min = 0, max =.5)
        
        
        # Sample Hyper-Parameter Means with C as fixpoint ----
        relCA <- runif(1, min = range_muA[1],max = range_muA[2])
        Mean_Cpar <- 8.5 #Mean_Cpar <- runif(1, min =range_muC[1], max= range_muC[2])
        Mean_Apar <- 4.2 #Mean_Apar <- Mean_Cpar*(relCA)
        Mean_Fpar <- runif(1, min =range_muF[1], max= range_muF[2])
        Mean_Fpar<-Mean_Fpar/(1-Mean_Fpar)
        Mean_bpar <- 0.1
        
        hyper_mus <- c(Mean_Cpar,Mean_Apar,Mean_Fpar,Mean_bpar)
        
        
        # Sample Variances and Set Covariances----
        
        sig_c <- runif(1, min = sigC[1],max = sigC[2])
        sig_a <- runif(1, min = sigA[1],max = sigA[2]) 
        sig_f <- 0.3 #sig_f <- runif(1, min = sigF[1],max= sigF[2])  # Wie genau ziehen, da sig von logit ?
        
        sig_b <- 0.001
        
        HyperPar_Cov <- diag(length(hyper_mus))
        
        HyperPar_Cov[1,1] <- (sig_c)^2
        HyperPar_Cov[2,2] <- (sig_a)^2
        HyperPar_Cov[3,3] <- (sig_f)^2
        HyperPar_Cov[4,4] <- (sig_b)^2
        
        
        # Set Correlations for Parameters ----
        
        HyperPar_Cov[1,2] = HyperPar_Cov[2,1] <- par_cor*(sqrt(HyperPar_Cov[1,1])*sqrt(HyperPar_Cov[2,2]))
        HyperPar_Cov[1,3] = HyperPar_Cov[3,1] <- par_cor*(sqrt(HyperPar_Cov[1,1])*sqrt(HyperPar_Cov[3,3]))
        HyperPar_Cov[2,3] = HyperPar_Cov[3,2] <- par_cor*(sqrt(HyperPar_Cov[2,2])*sqrt(HyperPar_Cov[3,3]))
        
        
        
        # Sample Parameters from MVN ----
        
        parms <- tmvtnorm::rtmvnorm(n=SampleSize, mean= hyper_mus, sigma=HyperPar_Cov, 
                                    lower=c(0,0,-Inf,0),upper = c(Inf,Inf,Inf,Inf))
        # Merge Parameters to one Matrix
        colnames(parms) <- c("conA","genA","f","baseA")
        parms[,4] <- 0.1
        parms[,3] <-dlogis(parms[,3])
        
        cor_emp_CA <-round(cor(parms[,1],parms[,2]),3)
        cor_emp_CF <-round(cor(parms[,1],parms[,3]),3)
        cor_emp_AF <-round(cor(parms[,2],parms[,3]),3)
        
        # Simulate Data for Estimation ----
        
        data <- simData_CSpan(parms,as.vector(respOpt_Cspan(conN,conK)),conRet)
        
        # Generate Stan Data ----
        a<-respOpt_Cspan(conN,conK)
        stan.dat <- list(count = data, 
                         K = ncol(data), N = SampleSize,
                         R = a[1,],
                         scale_b = 0.1)
        
        # Fit-The-Shit----
        
        fit_M3 <- stan(file = "M3_ComplexSpan.stan",data=stan.dat,
                       warmup = 1000, iter = 5000,
                       chains = 4,refresh = 1, cores=4)
        
        # Extract Posterior Parameters from fit object----
        
        post_samples <- rstan::extract(fit_M3, pars=c("c","a","f"), inc_warmup = F)
        
        means_a <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("a"))[,5])
        means_c <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("c"))[,5])
        means_f <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("f"))[,5])
        
        mean_est_hyper_a <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("mu_a"))[,5])
        mean_est_hyper_c <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("mu_c"))[,5])
        mean_est_hyper_f <- as.vector(rstan::get_posterior_mean(fit_M3, pars=c("mu_f"))[,5])
        
        
        
        # Calculate Highest Density Interval----
        
        HDI_c <- t(round(hdi(post_samples$c,credMass = .95),3))
        HDI_a <- t(round(hdi(post_samples$a, credMass = .95),3))
        HDI_f <-t(round(hdi(post_samples$f, credMass = .95),3))
        
        HDI_include_c <- c()
        HDI_include_a <- c()
        HDI_include_f <- c()
        
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
        
        # Calculate Posterior Mode ----
        
        mode_c <- posterior.mode(post_samples$c)
        mode_a <- posterior.mode(post_samples$a)
        mode_f <- posterior.mode(post_samples$f)
        
        
        # Extract rhat ----
        rhats_a<- as.data.frame(rhat(fit_M3, pars=c("a")))
        rhat_max_a <- max(rhats_a)
        
        rhats_c<- as.data.frame(rhat(fit_M3, pars=c("c")))
        rhat_max_c <- max(rhats_c)
        
        rhats_f<- as.data.frame(rhat(fit_M3, pars=c("f")))
        rhat_max_f <- max(rhats_f)
        
        # Calculate Residuals ----
        
        resC <- means_c - parms[,1]
        resA <- means_a - parms[,2]
        resF <- means_f - parms[,3]
        
        
        # Save Data in Object ----
        n <-1
        
        simulation <- data.frame("ID" = seq(1:SampleSize),
                                 "N"= conN,
                                 "K"= conK,
                                 "Retrievals" = conRet,
                                 "Cor_CA"=cor_emp_CA,
                                 "Cor_CF"= cor_emp_CF,
                                 "Cor_AF" = cor_emp_AF,
                                 "Repetition" = n,
                                 "RelCa" = relCA,
                                 "PC" = (data[,1] / nRetrievals),
                                 "hyper_mu_c"= Mean_Cpar,
                                 "hyper_sig_c" = sig_c,
                                 "hyper_mu_a"= Mean_Apar,
                                 "hyper_sig_a" = sig_a,
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
                                 "resF"  = resF,
                                 "max_rhat_f" = rhat_max_f,
                                 "mu_real_c" = parms[,1],
                                 "mu_real_a" = parms[,2],
                                 "mu_real_f" = parms[,3])
        
        
        benchmarks_subject_pars[m,] <-      cbind(conN, conK,conRet,par_cor,
                                                  meanPC = round(mean(data[,1] / conRet),2),
                                                  minPC = round(min(data[,1] / conRet),2),
                                                  maxPC = round(max(data[,1] / conRet),2),
                                                  Mean_Cpar,sig_c,Mean_Apar,sig_a,relCA,
                                                  Mean_Fpar, sig_f,
                                                  cor(simulation$mu_est_a,simulation$mu_real_a),
                                                  cor(simulation$mu_est_c,simulation$mu_real_c),
                                                  cor(simulation$mu_est_f,simulation$mu_real_f),
                                                  mean(simulation$mu_real_a-simulation$mu_est_a), 
                                                  sd(simulation$mu_real_a-simulation$mu_est_a),  
                                                  mean(simulation$mu_real_c-simulation$mu_est_c),
                                                  sd(simulation$mu_real_c-simulation$mu_est_c),
                                                  mean(simulation$mu_real_f-simulation$mu_est_f),
                                                  sd(simulation$mu_real_f-simulation$mu_est_f)) 
        
        
        hyper_pars[m,] <- cbind(n,conN,conK,conRet,
                                meanPC = round(mean(data[,1] / conRet),2),
                                minPC = round(min(data[,1] / conRet),2),
                                maxPC = round(max(data[,1] / conRet),2),
                                Mean_Apar, 
                                mean_est_hyper_a, 
                                Mean_Cpar,
                                mean_est_hyper_c,
                                Mean_Fpar,
                                mean_est_hyper_f)
        
        
        
        simulation.ss[[paste0("N_",conN,"_K_",conK,"_Retrievals_",conRet,"_Cor_",round(par_cor,2),"_Set_", n)]] <- simulation
        simulation.ss[[paste0("Benchmark Subject Pars")]] <- benchmarks_subject_pars[m,]
        saveRDS(simulation.ss, file=paste0("SimpleSpan_","N_",conN,"_K_",conK,"_Retrievals_",conRet,"_Cor_",round(cor_emp,2),"_Set_", n,".RDS"))
        
        
        
        if(m == n_sim)
          
        {
          hyper_pars <- as.data.frame(hyper_pars)
          
          benchmarks_hyper_pars <- hyper_pars %>% 
            group_by(N,K,Retrievals) %>% 
            summarise(cor_hyper_muC = cor(real_hyper_mu_c,est_hyper_mu_c),
                      cor_hyper_muA = cor(real_hyper_mu_a,est_hyper_mu_a),
                      cor_hyper_muF = cor(real_hyper_mu_f, est_hyper_mu_f),
                      mean_hyper_resC = mean(real_hyper_mu_c - est_hyper_mu_c),
                      sig_hyper_resC = sd(real_hyper_mu_c - est_hyper_mu_c),
                      mean_hyper_resA = mean(real_hyper_mu_a - est_hyper_mu_a),
                      sig_hyper_resA = sd(real_hyper_mu_a - est_hyper_mu_a),
                      mean_hyper_resF = mean(real_hyper_mu_f - est_hyper_mu_f),
                      sig_hyper_resF = sd(real_hyper_mu_f - est_hyper_mu_f)) 
          
          
          saveRDS(benchmarks_subject_pars, file=paste("SimpleSpan_Benchmarks_subject_pars",".RDS"))
          saveRDS(benchmarks_hyper_pars, file=paste("SimpleSpan_Benchmarks_hyper_pars",".RDS"))
          
        }
        
        
        m <- m + 1
        
        gc(full = T)
        
        pb$tick()
        
      }
      
    }
    
  }
  
}



toc()

print(fit_M3, pars=c("mu_a","mu_c","mu_f"))

post_samples <- rstan::extract(fit_M3, pars=c("a","c","f"), inc_warmup = F)

k <- 1
mat <- matrix(0,ncol=3,nrow = SampleSize)
for (i in 1:SampleSize)
{
  colnames(mat) <- c("meanC","meanA","meanF")                               
  
  meanC = mean(post_samples$c[,i])
  meanA = mean(post_samples$a[,i])
  meanF = mean(post_samples$f[,i])
  
  mat[k,] <- cbind(meanC,meanA,meanF)
  k <- 1 + k
  
  
}

real_pars<-parms                                                                                        

# Correlate Parameters----                                                                                                                  


paste("Cor C = ", round(cor(mat[,1],real_pars[,1]),2))
paste("Cor A = ", round(cor(mat[,2],real_pars[,2]),2))
paste("Cor F = ", round(cor(mat[,3],real_pars[,3]),2))
