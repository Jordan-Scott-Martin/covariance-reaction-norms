#code for running analysis in tutorial markdown file

setwd("...")

source("CRN functions.R")

#simulate data
stan.dl = sim_CRN_QG(N = 600, Nc = 30, npred = 3, ntrait = 3, l_es = 0.5, u_es = 0.9)

library(shinystan)
library(cmdstanr)
library(rstan)

#directory for cmdstan installation
set_cmdstan_path("C:/_Install/stan/Library/bin/cmdstan")

#compile model
CRN_mod = cmdstan_model(stan_file = "CRN_tutorial_mod.stan", 
                        stanc_options = list("O1"))

#estimate model
est = CRN_mod$sample(
  data = stan.dl, 
  iter_sampling = 1000, 
  iter_warmup = 1000, 
  init = 0.01, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.8, 
  max_treedepth = 10,
  refresh = 10)

saveRDS(est, "CRN_fit.RDS")
launch_shinystan(est)
