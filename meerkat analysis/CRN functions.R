list.of.packages = c("reshape2", "ggplot2","tidybayes","plyr","mvtnorm","rethinking", "devtools", "posterior")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if("rethinking" %in% new.packages) devtools::install_github("rmcelreath/rethinking")

library(reshape2)
library(ggplot2)
library(tidybayes)
library(plyr)
library(Matrix)
library(mvtnorm)
library(rethinking)
library(posterior)

#this function generates single measure/subject data for a
#quantitative genetic CRN analysis using a supplied matrix 
#of environmental predictors (continuous and/or discrete)
#The function assumes balanced sampling and so requires 
#input values that divide evenly across conditions
#N = total number of individuals (even number)
#X = design matrix of covariates
#l_es = lower effect size (absolute value)
#u_es = upper effect size (absolute value)
sim_CRN_QG = function(N, Nc, npred, ntrait, l_es, u_es, standl = TRUE){
  try(if(N/Nc %% 2==1) stop("N/Nc must be even for balanced sampling."))
  try(if(l_es < 0) stop("l_es is an absolute value. Please input a positive number."))
  try(if(u_es < 0) stop("u_es is an absolute value. Please input a positive number."))
  #generate values for the data and model parameters defined in the model
  #see Stan file for details
  
  #generate covariate matrix
  df = data.frame(rmvnorm(Nc, mean = rep(0,npred)))
  df = df[rep(seq_len(nrow(df)), each = N/Nc), ]
  
  #make context-level design matrix
  contexts = unique(df) #unique environmental contexts
  suppressWarnings({
  X = data.frame(model.matrix(formula(c("~",paste0(colnames(df), collapse = " + "))), 
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
  npred = npred + 1
  ncor = (ntrait*(ntrait-1))/2; #number of unique correlations

  B_m = matrix(runif(npred*ntrait, min = l_es, max = u_es)
          * sample(c(-1,1), npred*ntrait, replace = TRUE, prob = c(0.5,0.5)),
               nrow = npred, ncol = ntrait) #mean RN
  B_v = matrix(runif(npred*ntrait, min = l_es, max = u_es) 
               * sample(c(-1,1), npred*ntrait, replace = TRUE, prob = c(0.5,0.5)),
               nrow = npred, ncol = ntrait) #variance RN
  B_cpc = matrix(runif(npred*ncor, min = l_es, max = u_es)
                 * sample(c(-1,1), npred*ncor, replace = TRUE, prob = c(0.5,0.5)),
                 nrow = npred, ncol = ncor) #correlation RN
  
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
  E_cv = diag(rep(1,ntrait)) %*% E %*% diag(rep(1,ntrait))
  
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
  df$idc = match(paste(df$id, df$c_id,sep="."), paste(corder$id,corder$c,sep="."))
  
  if(standl == TRUE){
  
  #stan data list
  return(list(N = nrow(df), C = nrow(X), I = max(df$id),
       D = ntrait, P = ncol(X), 
       id = df$id, c_id = df$c_id, idc = df$idc,
       X = X, A = A, cmat = cmat, cm = ncol(cmat),cn = cmat_n,
       cnt = sum(cmat_n), z = df[,c(paste0("z",1:ntrait))],
       B_m = B_m, B_v = B_v, B_cpc = B_cpc, Gcov = Gcov,
       sd_E = rep(1,ntrait), E = t(chol(E)) %*% chol(E)))
  }
  if(standl == FALSE){
    standf = data.frame(id = df$id, df[,c(paste0("z",1:ntrait))], X[,-1])
    rownames(standf) = 1:nrow(standf)
    return(standf)
  }
}

#this function simplifies the extraction of posterior samples
#from cmdstan objects to match traditional rstan functions
extract = function(fit_obj) {
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

#functions below are necessary for predicting environmental effects
#on trait variances, correlations, and covariances using the 
#posterior MCMC samples from a CRN Stan model

#this function generates Cholesky factorized correlation matrices
#from posterior samples of canonical partial correlations
#the function is used within other functions below
lkj_to_chol_corr = function(constrained_reals, ntrait) {
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

#this function generates posterior samples
#of trait variances predicted by environmental covariates
#x = input prediction matrix
#traits = vector of names (must follow order of model)
#v_b = matrix of beta coefficients for variances
v_f = function(x, traits, v_b){
  v = list()
  for(r in 1:nrow(v_b)){
    v[[r]] = exp(as.matrix(x) %*% v_b[r,,] )
  }
  v = melt(v)
  v = v[, -which(colnames(v) %in% c("L1"))]
  colnames(v) = c("npred","element", "var")
  ntrait = length(traits)
  key = data.frame(old_name = seq(1:ntrait), new_name = traits)
  v$element = setNames(key$new_name, key$old_name)[v$element]
  v = v[order(v$npred),]
  v = cbind(v, x[rep(seq_len(nrow(x)), each = length(traits) * nrow(v_b)), ])
  return(v)
}

#this function generates posterior samples
#of canonical partial correlations 
#predicted by environmental covariates
#x = input prediction matrix
#post_cpc_b = matrix of beta coefficients for partial correlations
cpc_f = function(x, cpc_b){
  cpc = list()
  for(r in 1:nrow(cpc_b)){
    cpc[[r]] = tanh( as.matrix(x) %*% cpc_b[r,,])}
  return(cpc)
}

#this function generates posterior samples
#of correlations predicted by environmental covariates
#x = input prediction matrix
#traits = vector of names (must follow order of model)
#cpc = output of cpc_f (predicted partial correlations)
cor_f = function(x, traits, cpc){
  ntrait = length(traits)
  ncor = (ntrait*(ntrait-1))/2
  gcor_cpc = list()
  for(r in 1:length(cpc)){
    cor_p = list()
    cpc_p = cpc[[r]]
    for(p in 1:nrow(cpc_p)){
      lcor =  lkj_to_chol_corr(cpc_p[p,], ntrait)
      cor_p[[p]] = lcor %*% t(lcor)
    }
    gcor_cpc[[r]] = cor_p
  }
  cor_p = list()
  for(i in 1:length(gcor_cpc)){
    temp = list()
    for(p in 1:nrow(x)){
    m = gcor_cpc[[i]][[p]]
    temp[[p]] = data.frame(element1 = which(upper.tri(m),arr.ind=TRUE)[,1],
                           element2 = which(upper.tri(m),arr.ind=TRUE)[,2],
                           correlation = m[upper.tri(m)], npred = p) 
    }
    df = do.call(rbind.data.frame, c(temp, row.names = NULL))
    cor_p[[i]] = df
    }
  cor_p = do.call(rbind.data.frame, cor_p)
  cor_pl = melt(cor_p, id.vars = c("npred","element1","element2"))
  cor_pl = cor_pl[order(cor_pl$npred),]
  
  key = data.frame(old_name = seq(1:ntrait), new_name = traits)
  cor_pl$element1 = setNames(key$new_name, key$old_name)[cor_pl$element1]
  cor_pl$element2 = setNames(key$new_name, key$old_name)[cor_pl$element2]
  cor_pl$element = apply(cor_pl[,c("element1","element2")], 1, paste0, collapse = "_")
  cor_pl = cbind(cor_pl[order(cor_pl$npred),], x[rep(seq_len(nrow(x)), each = ncor * length(cpc)), ], row.names = NULL)
  return(cor_pl)
}

#this function generates posterior samples
#of covariances predicted by environmental covariates
#x = input prediction matrix
#traits = vector of names (must follow order of model)
#v_b = matrix of beta coefficients for variances
#cpc = output of cpc_f (predicted partial correlations)
cov_f = function(x, traits, v_b, cpc){
  ntrait = length(traits)
  ncor = (ntrait*(ntrait-1))/2
  v = list()
  for(r in 1:nrow(v_b)){
    v[[r]] = exp(as.matrix(x) %*% v_b[r,,] )
    }
  gcovl = list()
  for(r in 1:length(cpc)){
    cov_p = list()
    cpc_p = cpc[[r]]
    v_p = v[[r]]
    sds = sqrt(v_p)
    for(p in 1:nrow(cpc_p)){
      lcor =  lkj_to_chol_corr(cpc_p[p,], ntrait)
      cor_p = lcor %*% t(lcor)
      cov_p[[p]] = diag(sds[p,]) %*% cor_p %*% diag(sds[p,])
    }
    gcovl[[r]] = cov_p
  }
  cov_p = list()
  for(i in 1:length(gcovl)){
    temp = list()
    for(p in 1:nrow(x)){
    m = gcovl[[i]][[p]]
    temp[[p]] = data.frame(element1 = which(upper.tri(m),arr.ind=TRUE)[,1],
                           element2 = which(upper.tri(m),arr.ind=TRUE)[,2],
                           covariance = m[upper.tri(m)], npred = p) 
    }
    df = do.call(rbind.data.frame, c(temp, make.row.names = F))
    cov_p[[i]] = df
    }
  cov_p = do.call(rbind.data.frame, c(cov_p, row.names = NULL))
  cov_pl = melt(cov_p, id.vars = c("npred","element1","element2"))
  cov_pl = cov_pl[order(cov_pl$npred),]
  
  key = data.frame(old_name = seq(1:ntrait), new_name = traits)
  cov_pl$element1 = setNames(key$new_name, key$old_name)[cov_pl$element1]
  cov_pl$element2 = setNames(key$new_name, key$old_name)[cov_pl$element2]
  cov_pl$element = apply(cov_pl[,c("element1","element2")], 1, paste0, collapse = "_")
  
  cov_pl = cbind(cov_pl[order(cov_pl$npred),], x[rep(seq_len(nrow(x)), each = ncor * length(cpc)), ], row.names = NULL)
  return(cov_pl)
}

#this function generates beta coefficients
#for environmental effects on variances in long format
#x = environmental covariate matrix (incl. intercept)
#traits = vector of names (must follow order of model)
#post_v_b = matrix of beta coefficients for variances
v_beta_f = function(x, traits, v_b){
  v_list = alply(v_b,1)
  v_list = lapply(v_list, function(s) {
                  colnames(s) = traits
                  rownames(s) = colnames(x)
                  return(s)
                  })
  pv_b = melt(v_list)
  colnames(pv_b) = c("v_predictor", "trait","value", "L1")
  pv_b = reshape(pv_b,
           idvar = c("trait","L1"),
           timevar = "v_predictor",
           direction = "wide")
  pv_b = pv_b[, -which(colnames(pv_b) %in% c("L1"))]
  colnames(pv_b) = c("element", colnames(x))
  return(pv_b)
  }
  
#this function generates beta coefficients
#for environmental effects on correlations
#from a CRN model specified with partial correlations
#x = environmental covariate matrix
#traits = vector of names (must follow order of model)
#post_cpc_b = matrix of beta coefficients for partial correlations
cor_beta_f = function(x, traits, cpc_b){
  ntrait = ncol(cpc_b[1,,])
  pred.n = colnames(x)[-1]
  x.p = matrix(0, nrow = (ncol(x)-1)*2, ncol = ncol(x))
  x.p[,1] = 1
  for(c in 2:ncol(x.p)){x.p[2*(c-1),c] = 1}
  
  mv_cpc = cpc_f(x.p, cpc_b)
  mv_cor = cor_f(x.p, traits, mv_cpc)
  key = data.frame(old_name = seq(1:nrow(x.p)), new_name = paste0(rep(pred.n, each=2), c(0,1)))
  mv_cor$npred = setNames(key$new_name, key$old_name)[mv_cor$npred]
  cors = unique(mv_cor$element)
  
  predfl = list()
  for(j in 1:length(cors)){
  predl = list()
  temp = 1
  for(p in 1:length(pred.n)){
  predl[[p]] = atanh(mv_cor[mv_cor$npred==key$new_name[2*p] & 
                            mv_cor$element==cors[[j]] ,"value"]) -
                atanh(mv_cor[mv_cor$npred==key$new_name[temp] &
                             mv_cor$element==cors[[j]],"value"])
  temp = temp + 2
  }
  predf = data.frame(predl)
  colnames(predf) = pred.n
  predf$element = cors[[j]]
  predf$X.Intercept. = atanh(mv_cor[mv_cor$npred==key$new_name[[1]] & 
                                      mv_cor$element==cors[[j]] ,"value"])
  predfl[[j]] = predf
  }
  predcor = do.call(rbind.data.frame, c(predfl, make.row.names = F))
  predcor = predcor[,c("element","X.Intercept.",pred.n)]
  return(predcor)
}


