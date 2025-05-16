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

#function for simulating relatedness (A) matrix
pedfun = function(popmin, popmax, ngenerations,
                  epm, nonb, nids, I, missing=FALSE){
  
  # get list of individuals and their generations
  gener=1:ngenerations
  
  genern = rep(1:ngenerations, times = nids)
  ID = 1:sum(nids)
  
  # runs on generation-by-generation basis
  for(i in 1:ngenerations){
    
    id=ID[which(genern==i)]
    dam=rep(NA, nids[i])
    sire=rep(NA, nids[i])
    
    # randomly allocates sex (0 = male, 1 = female)
    sex=sample(c(0,1), length(id), replace=TRUE)
    
    # for first generation, no dams or sires are known 
    # so remain NA
    
    if(i==1){
      
      # combine into single data frame
      pedigree=data.frame(id=id, dam=dam, sire=sire, 
                          generation=i, sex=sex)
      
    }
    
    else if(i>1){
      
      # for all generations after first
      # list of all possible dams and sires
      # from previous generation
      pdams=pedigree$id[which(pedigree$generation==(i-1) &
                                pedigree$sex==1)]
      psires=pedigree$id[which(pedigree$generation==(i-1) &
                                 pedigree$sex==0)]
      
      # determine number of pairs
      # depending on how many males and females
      # and the proportion of the population that is non-breeding
      npairs=min(length(pdams), length(psires)) - 
        round(min(length(pdams), length(psires))*nonb)
      
      # selects breeding males and females
      pdams=pedigree$id[which(pedigree$generation==(i-1) & 
                                pedigree$sex==1)]
      psires=pedigree$id[which(pedigree$generation==(i-1) & 
                                 pedigree$sex==0)]
      
      if(length(pdams)<npairs | length(psires)<npairs){
        npairs=min(length(pdams), length(psires))
      }
      
      # selects pairs from possible dams and sires
      pairs=data.frame(dam=sample(pdams, npairs, replace=FALSE),
                       sire=sample(psires, npairs, replace=FALSE))
      # gives each offspring their parental pair
      pairid=as.numeric(sample(rownames(pairs), 
                               length(id), replace=TRUE))
      
      # gives each offspring their sex
      sex=sample(c(0,1), length(id), replace=TRUE)
      
      # put into dataframe format
      addped=data.frame(id=id, 
                        dam=pairs$dam[pairid], 
                        sire=pairs$sire[pairid],
                        generation=i, 
                        sex=sex)
      
      
      # deals with extra-pair mating (if included)
      if(!is.null(epm)){
        
        # for each individual, sample if they are extra pair
        # if 0 not extra pair
        # if 1 sire resampled from breeding population
        # if 2 dam resampled
        ext=sample(c(0,1,2), nrow(addped), 
                   replace=TRUE, 
                   prob = c(1-epm, epm/2, epm/2))
        for(j in 1:nrow(addped)){
          if(ext[j]>0){
            if(ext[j]==1){
              addped$sire[j]=sample(psires,1,replace=TRUE)
            }else if (ext[j]==2){
              addped$dam[j]=sample(pdams,1,replace=TRUE)
            }
            
          }
        }
      }
      
      
      # add new generation to the whole pedigree
      pedigree=rbind(pedigree, addped)
    }
    
  }
  
  ped = pedigree
  
  # make id's non-numeric
  ped$id=paste("ID",ped$id, sep="")
  ped$dam[which(!is.na(ped$dam))]=paste("ID",ped$dam[which(!is.na(ped$dam))], sep="")
  ped$sire[which(!is.na(ped$sire))]=paste("ID",ped$sire[which(!is.na(ped$sire))], sep="")
  ped$id=as.character(ped$id)
  ped$dam=as.character(ped$dam)
  ped$sire=as.character(ped$sire)
  
  IDs = sample(ped[ped$generation==ngenerations, "id"], I, replace=FALSE)
  ped = MCMCglmm::prunePed(ped, keep = IDs, make.base=TRUE)
  inv.phylo = MCMCglmm::inverseA(ped[,c("id","dam","sire")])
  A = solve(inv.phylo$Ainv)
  A = cov2cor(A)
  A = (A + t(A))/2 # Not always symmetric after inversion
  A = as.matrix(A)
  rownames(A) = rownames(inv.phylo$Ainv)
  colnames(A) = rownames(inv.phylo$Ainv)
  
  #subset to final generation
  A_sub=A[IDs,IDs]
  A_mat = as.matrix(nearPD(A_sub)$mat)
  A_mat = cov2cor(A_mat)
  return(A_mat)
}

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
  Gcor2 = list()
  Gcov = list()
  for(c in 1:nrow(X)){
    mat = diag(1, ntrait)
    mat[lower.tri(mat)] = G_cor[[c]]
    mat[upper.tri(mat)] = mat[lower.tri(mat)]
    Gcor2[[c]] = mat
    Gcov[[c]] = diag(G_sd[c,]) %*% Gcor2[[c]] %*% diag(G_sd[c,])
  }
  
  #relatedness matrix
  #population properties
  popmin = N
  popmax = N*1.25
  ngenerations = 10
  nids = sample(popmin:popmax, ngenerations, replace=TRUE) #N / generation
  epm = sample(seq(0.15, 0.25,by=0.05),1) #extra-pair mating
  nonb = sample(seq(0.3,0.6,by=0.05),1) #proportion of non-breeding / generation
  A = pedfun(popmin=popmin, popmax=popmax, ngenerations=ngenerations,
             epm=epm, nonb=nonb, nids=nids, I=N, missing=FALSE)
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
  sd_E = rep(1,ntrait)
  E_cv = diag(sd_E) %*% E %*% diag(sd_E)
  
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
  
  #stan data list
  list(
    variables=list(B_m = B_m, B_v = B_v, B_cpc = B_cpc, 
                   sd_E = sd_E, L_E = t(chol(E)), Z_G = G_z),
    
    generated=list(N = nrow(df), C = nrow(X), I = max(df$id),
                   D = ntrait, P = ncol(X), sd_E = sd_E, 
                   id = df$id, c_id = df$c_id, idc = df$idc,
                   X = X, A = A, cmat = cmat, cm = ncol(cmat),cn = cmat_n,
                   cnt = sum(cmat_n), z = df[,c(paste0("z",1:ntrait))]))

}

#generate objects for SBC
set.seed(9)
n_sims = 200  # Number of sim per run

start = Sys.time()
SBCf1 = SBC_generator_function(generator_function1, N = 300, Nc = 10, npred = 2, ntrait = 3)
df1 = generate_datasets(SBCf1, n_sims)
saveRDS(df1, "df1.RDS")
backend1 = SBC_backend_cmdstan_sample(SBC_mod1, iter_sampling = 500, iter_warmup = 1000, 
                                      init = 0.001, chains = 4,  parallel_chains = 4, 
                                      adapt_delta = 0.90, max_treedepth = 10)
#conduct SBC
results1 =  compute_SBC(datasets = df1, backend = backend1, cache_mode = "results", 
                        cache_location = file.path(cache_dir, "results"),
                        keep_fits = TRUE)
saveRDS(results1, "results_SBC_CRN.RDS")
end = Sys.time()
end - start

#load results
results1 = readRDS("results_SBC_CRN.RDS")

#plot results
bdiff = 
  plot_ecdf_diff(results1, variables = paste0("B_v","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)
cdiff = 
  plot_ecdf_diff(results1, variables = paste0("B_cpc","[",rep(1:4,3),",",rep(1:3,each=4),"]"), prob = 0.95)

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
save_plot("fig s1 results.png", cbp, base_height = 2.5, base_width = 5.5)

