#############################################################
#prep workspace

#set directory
setwd("...")

#load packages
library(SBC)
library(rstan)
library(rethinking)
library(plyr)
library(mvtnorm)
library(cmdstanr)
library(ggplot2)

set_cmdstan_path("...")

#SBC cmdstan #################################################

#cache for SBC in RStan
cache_dir = "./_basic_usage_cmdstan_SBC_cache"
if(!dir.exists(cache_dir)) {dir.create(cache_dir) }

#utility function for canonical partial correlations
lkj_to_chol_corr <- function(constrained_reals, ntrait) {
  x = matrix(0, nrow=ntrait, ncol=ntrait)
  x[1,1] = 1

  z_counter = 1
  for(i in 2:ntrait) {
    for(j in 1:ntrait) {
      if(i==j) {
        x[i,j] = sqrt(1 - (sum(x[i,1:j]^2)))
      } else if(i > j) {
        x[i,j] = constrained_reals[z_counter]*sqrt(1 - (sum(x[i,1:j]^2)))
        z_counter = z_counter + 1
      }
    }
  }

  return(x)
}

#load Gaussian NLS selection model
SBC_mod1 = cmdstan_model(stan_file = "SBC_CRN_mod.stan", stanc_options = list("O1"))

#function for generating datasets and parameter values
generator_function1 = function(N, Nx1, ntrait){
  
  
  #generate values for the data and model parameters defined in the model
  #see Stan file for details
  
  #generate covariates
  x1 = rnorm(Nx1, 0, 1) #standardized continuous predictor
  x1 = rep(x1, each = N/Nx1) #10 individuals / continuous context
  x2 = c(0,1) #binary predictor
  x2 = rep(rep(x2, each = N/Nx1/2), Nx1) #5 each discrete context / cont. context
  df = data.frame(cbind(x1,x2))
  
  #make context-level design matrix
  contexts = unique(df) #unique environmental contexts
  X = data.frame(model.matrix(1:nrow(contexts) ~ x1*x2, data = contexts))
  rownames(X) = 1:nrow(X)
  X
  
  #index contexts per measurement
  c_id = list()
  for(c in 1:nrow(contexts)){
    c_id[[c]] = as.integer(rownames(match_df(df, contexts[c,]))) 
  }
  c_id  = data.frame(context=rep(seq_along(c_id),lengths(c_id)), row=unlist(c_id))
  df$c_id = c_id[order(c_id$row),"context"]
  
  #CRN parameters
  npred = ncol(X) #number of mean and variance params / response
  ncor = (ntrait*(ntrait-1))/2; #number of unique correlations
  
  #effect sizes
  #B_m = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait)
  #B_v = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait) #variance RN
  #B_cpc = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ncor) #correlation RN
  
  B_m = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait) #mean RN
  B_v = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait) #variance RN
  B_cpc = matrix(rnorm(npred*ncor, 0,1), nrow = npred, ncol = ncor) #correlation RN
  
  #generate means, standard deviations, and correlations
  mu = as.matrix(X) %*% (B_m)
  G_sd = sqrt(exp(as.matrix(X) %*% B_v))
  G_cpc = tanh(as.matrix(X) %*% B_cpc)
  
  G_cor = list()
  for(r in 1:nrow(G_cpc)){
    LR = lkj_to_chol_corr(G_cpc[r,], ntrait)
    R = LR %*% t(LR)
    G_cor[[r]] = R[lower.tri(R)]
  }
  
  #construct G correlation and covariance matrices
  Gcor = list()
  Gcov = list()
  for(c in 1:nrow(X)){
    mat = diag(1, ntrait)
    mat[lower.tri(mat)] = G_cor[[c]]
    mat[upper.tri(mat)] = mat[lower.tri(mat)]
    Gcor[[c]] = mat
    Gcov[[c]] = diag(G_sd[c,]) %*% mat %*% diag(G_sd[c,])
  }
  
  #relatedness matrix
  A = rlkjcorr(1, N, 1)
  LA = t(chol(A))
  seq = N/length(Gcov)
  G_z = rmvnorm(N, mean=rep(0,ntrait), sigma=diag(rep(1,ntrait)))
  G_z = (G_z - apply(G_z, 2, mean))/apply(G_z, 2, sd)
  G_s = list()
  p = 1
  for(i in 1:length(Gcor)){
    G_s[[i]] = LA[p:(p+seq-1),p:(p+seq-1)] %*% G_z[p:(p+seq-1),] %*% chol(Gcov[[i]])
    p = p + seq
  }
  G_vals = Reduce(rbind.data.frame, G_s)
  colnames(G_vals) = paste0("a",1:ntrait)
  
  #add context-specific values to dataframe
  df$id = 1:N
  df = cbind(df, G_vals)
  
  #generate residuals
  E = rlkjcorr(1, ntrait, 10)
  E_sd = rexp(ntrait, 2)
  E_cv = diag(E_sd) %*% E %*% diag(E_sd)
  
  res = rmvnorm(N, mean=rep(0,ntrait), sigma=E_cv)
  
  #generate responses
  z = mu[df$c_id,] + df[,c(paste0("a",1:ntrait))] + res
  colnames(z) = paste0("z",1:ntrait)
  df = cbind(df, z)
  
  #matrix indexing subjects in contexts
  c_list = list()
  for(c in 1:nrow(X)){
    c_list[[c]] = unique(df[df$c_id==c,"id"])
  }
  cmat = matrix(t(sapply(c_list, function(x) x[1:max(lengths(c_list))])),
                ncol = max(lengths(c_list)))
  cmat = t(apply(cmat, 1, function(x) `length<-`(na.omit(x), length(x))))
  cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) ) #N/context
  cmat[is.na(cmat)] = 0 #remove NAs
  
  #new index corresponding to subject position in cmat
  temp = t(cmat)
  corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))
  df$idc = match(paste0(df$id, df$c_id,sep="."), paste0(corder$id,corder$c,sep="."))
  
  #stan data list
  list(
    variables=list(B_m = B_m, B_v = B_v, B_cpc = B_cpc, 
                   sd_E = E_sd, L_E = t(chol(E)), Z_G = G_z),
    
    generated=list(N = nrow(df), C = nrow(X), I = max(df$id),
                   D = ntrait, P = ncol(X), 
                   id = df$id, c_id = df$c_id, idc = df$idc,
                   X = X, A = A, cmat = cmat, cm = ncol(cmat),cn = cmat_n,
                   cnt = sum(cmat_n), z = df[,c(paste0("z",1:ntrait))]))

}

#generate objects for SBC
set.seed(9)
n_sims = 100  # Number of sim per run

SBCf1 = SBC_generator_function(generator_function1, N = 100, Nx1 = 10, ntrait = 3)
df1 = generate_datasets(SBCf1, n_sims)
backend1 = SBC_backend_cmdstan_sample(SBC_mod1, iter_sampling = 500, iter_warmup = 1000, 
                                      init = 0.01, chains = 5,  parallel_chains = 5, 
                                      adapt_delta = 0.80, max_treedepth = 10) 
#conduct SBC
results1 =  compute_SBC(datasets = df1, backend = backend1, cache_mode = "results", 
                        cache_location = file.path(cache_dir, "results"),
                        keep_fits = FALSE)
saveRDS(results1, "results_SBC_CRN.RDS")

#load results
results1 = readRDS("results_SBC_CRN.RDS")

#plot results
bdiff = plot_ecdf_diff(results1, variables = paste0("B_v","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)
cdiff = plot_ecdf_diff(results1, variables = paste0("B_cpc","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)

bdiff$facet$params$nrow=4
bdiff$facet$params$ncol=3
bdiff = bdiff + theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
              axis.ticks.y=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_text(size=12,face="bold"),
              axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.line = element_line(linewidth = 1),
              
              panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                        linetype="solid"),
              strip.text = element_blank(),
              strip.background = element_blank(),
              panel.background= element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
              guides(fill="none", color ="none")

cdiff$facet$params$nrow=4
cdiff$facet$params$ncol=3
cdiff = cdiff + theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
                      axis.ticks.y=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.x=element_text(size=12,face="bold"),
                      axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.line = element_line(linewidth = 1),
                      
                      panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                                linetype="solid"),
                      strip.text = element_blank(),
                      strip.background = element_blank(),
                      panel.background= element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides(fill="none", color ="none")

library(cowplot)
cbp = plot_grid(bdiff, cdiff, ncol = 2, rel_widths = c(1,1))
save_plot("fig 2 crn.png", cbp, base_height = 3.5, base_width = 6)

