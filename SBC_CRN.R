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
library(Matrix)

#directory for cmdstan installation
set_cmdstan_path("...")


#Test model #################################################

#load Gaussian NLS selection model
SBC_mod1 = cmdstan_model(stan_file = "SBC_CRN_mod.stan", 
                         stanc_options = list("O1"))

res = SBC_mod1$sample(data = df1$generated[[1]], iter_sampling = 2000, iter_warmup = 1000,
                      parallel_chains = 5, init = 0.001, adapt_delta = 0.99)

shinystan::launch_shinystan(res)

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
generator_function1 = function(N, Nc, npred, ntrait){
  
  
  #generate values for the data and model parameters defined in the model
  #see Stan file for details
  
  #generate covariate matrix
  df = data.frame(rmvnorm(Nc, mean = rep(0,npred)))
  df = df[rep(seq_len(nrow(df)), each = N/Nc), ]
  
  #make context-level design matrix
  contexts = unique(df) #unique environmental contexts
  suppressWarnings({
    X = data.frame(model.matrix(formula(c("~",paste0(colnames(df), collapse = " * "))), 
                                data = contexts))
  })
  rownames(X) = 1:nrow(X)
  
  #index contexts per measurement
  suppressMessages({
    c_id = list()
    for(c in 1:nrow(contexts)){
      c_id[[c]] = as.integer(rownames(match_df(df, contexts[c,]))) 
    }
    c_id  = data.frame(context=rep(seq_along(c_id),lengths(c_id)), row=unlist(c_id))
    df$c_id = c_id[order(c_id$row),"context"]
  })
  
  #CRN parameters
  npred = npred + (npred*(npred-1))/2 + 1
  ncor = (ntrait*(ntrait-1))/2; #number of unique correlations
  
  #effect sizes
  B_m = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait) #mean RN
  B_v = matrix(rnorm(npred*ntrait, 0,1), nrow = npred, ncol = ntrait) #variance RN
  B_cpc = matrix(rnorm(npred*ncor, 0,1), nrow = npred, ncol = ncor) #correlation RN
  
  #generate means, standard deviations, and correlations
  mu = as.matrix(X) %*% B_m
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
  sd_E = 1
  E_cv = diag(rep(sd_E,ntrait)) %*% E %*% diag(rep(sd_E,ntrait))
  
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
                   sd_E = rep(sd_E,ntrait), L_E = t(chol(E)), Z_G = G_z),
    
    generated=list(N = nrow(df), C = nrow(X), I = max(df$id),
                   D = ntrait, P = ncol(X), sd_E = sd_E, 
                   id = df$id, c_id = df$c_id, idc = df$idc,
                   X = X, A = A, cmat = cmat, cm = ncol(cmat),cn = cmat_n,
                   cnt = sum(cmat_n), z = df[,c(paste0("z",1:ntrait))]))

}

#generate objects for SBC
set.seed(9)
n_sims = 100  # Number of sim per run

start = Sys.time()
SBCf1 = SBC_generator_function(generator_function1, N = 200, Nc = 10, npred = 2, ntrait = 3)
df1 = generate_datasets(SBCf1, n_sims)
backend1 = SBC_backend_cmdstan_sample(SBC_mod1, iter_sampling = 1000, iter_warmup = 1000, 
                                      init = 0.001, chains = 5,  parallel_chains = 5, 
                                      adapt_delta = 0.99, max_treedepth = 10)
#conduct SBC
results1 =  compute_SBC(datasets = df1, backend = backend1, cache_mode = "results", 
                        cache_location = file.path(cache_dir, "results"),
                        keep_fits = FALSE)
saveRDS(results1, "results_SBC_CRN.RDS")
end = Sys.time()
end - start

#load results
results1 = readRDS("results_SBC_CRN.RDS")

#subset to remove problematic simulations (<1% divergence)
keep = results1$backend_diagnostics$sim_id[
  results1$backend_diagnostics$n_divergent < .001*1000*5]
results1nd = results1[keep]

#plot results
plot_rank_hist(results1, variables = paste0("B_v","[",rep(1:4,3),",",rep(1:3,each=4),"]"), 
               prob = 0.95)
plot_ecdf_diff(results1, variables = paste0("B_v","[",rep(1:4,3),",",rep(1:3,each=4),"]"), 
               prob = 0.95, K = "max")
plot_ecdf_diff(results1, variables = paste0("B_cpc","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)


bdiff = plot_ecdf_diff(results1, variables = paste0("B_v","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)
cdiff = plot_ecdf_diff(results1, variables = paste0("B_cpc","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)

bdiff$facet$params$nrow=4
bdiff$facet$params$ncol=3
bdiff = 
  bdiff + coord_cartesian(ylim=c(-0.25,0.25)) +
              scale_x_continuous(breaks = c(0, 0.5, 1))+
              scale_y_continuous(breaks = c(-0.2, 0, 0.2))+
              theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
              axis.text.x=element_text(size=7),
              axis.text.y=element_text(size=7),
              axis.line = element_line(linewidth = 1),
              panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                        linetype="solid"),
              strip.text = element_blank(),
              strip.background = element_blank(),
              panel.background= element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
              guides(fill="none", color ="none")

cdiff$facet$params$nrow=4
cdiff$facet$params$ncol=3
cdiff = 
  cdiff + coord_cartesian(ylim=c(-0.25,0.25))+
            scale_x_continuous(breaks = c(0, 0.5, 1))+
            scale_y_continuous(breaks = c(-0.2, 0, 0.2))+
            theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
                      axis.text.x=element_text(size=7),
                      axis.text.y=element_text(size=7),
                      axis.line = element_line(linewidth = 1),
                      panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                                linetype="solid"),
                      strip.text = element_blank(),
                      strip.background = element_blank(),
                      panel.background= element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
  guides(fill="none", color ="none")

library(cowplot)
cbp = plot_grid(bdiff, cdiff, ncol = 2, rel_widths = c(1,1))
save_plot("fig 2 crn_v2.png", cbp, base_height = 2.5, base_width = 5.5)

