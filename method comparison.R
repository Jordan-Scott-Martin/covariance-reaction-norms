#set working directory
setwd(...)

#load functions and packages
source("CRN functions.R")
library(cmdstanr)

#directory for cmdstan installation
set_cmdstan_path("...")

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
generator_function = function(N, years, es){
  
  #temporal variation
  hx = 1 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = years) + 0.1*1:years
  lx = 0 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = years) + 0.1*1:years
  x = c(hx,lx)
  x = scale(x)
  plot(x[1:length(hx)], col = "red", ylim = c(-2,2), type = "l")
  lines(x[(length(hx)+1):length(x)], col = "blue")
  
  s.temp = data.frame(hx, lx, year = 1:years )
    y = sample(1:years, N, replace =T)
  df = data.frame(x1 = 1, x2 = scale(as.vector(t(s.temp[y, 1:2]))),
                  year = rep(y,each=2), season = rep(c("summer","winter"),N))
  
  #make context-level design matrix
  contexts = unique(df[,1:2]) 
  X = contexts
  rownames(X) = 1:nrow(X)
  
  #index contexts per measurement
  suppressMessages({
    c_id = list()
    for(c in 1:nrow(X)){
      c_id[[c]] = as.integer(rownames(match_df(df, X[c,]))) 
    }
    c_id  = data.frame(context=rep(seq_along(c_id),lengths(c_id)), row=unlist(c_id))
    df$c_id = c_id[order(c_id$row),"context"]
  })
  
  #CRN parameters
  np = 2
  ncor = 1
  
  #effect sizes
  ves = es[1]
  ces = es[2]
  B_m = matrix(c(0, 0, 0, 0), nrow = np, ncol = 2)
  B_v = matrix(c(-1.2, ves, -1.2, ves), nrow = np, ncol = 2) 
  B_cpc = matrix(c(0, ces), nrow = np, ncol = 1) 
  
  #generate means, standard deviations, and correlations
  mu = as.matrix(X) %*% B_m
  G_sd = sqrt(exp(as.matrix(X) %*% B_v))
  G_sdtv = sqrt(exp(matrix(c(1,1,0,1), nrow = 2,ncol = 2) %*% B_v))
  G_vdt = c(G_sdtv[2,1]^2-G_sdtv[1,1]^2, G_sdtv[2,2]^2-G_sdtv[1,2]^2)
  G_cpc = tanh(as.matrix(X) %*% B_cpc)
  G_ctv = tanh(matrix(c(1,1,0,1), nrow = 2,ncol = 2) %*% B_cpc)
  G_cdt = G_ctv[2,1]-G_ctv[1,1]
  G_cvdt = G_ctv[2,1] * G_sdtv[2,1] * G_sdtv[2,2] -
    G_ctv[1,1] * G_sdtv[1,1] * G_sdtv[1,2] 
  
  G_cor = list()
  for(r in 1:nrow(G_cpc)){
    LR = lkj_to_chol_corr(G_cpc[r,], 2)
    R = LR %*% t(LR)
    G_cor[[r]] = R[lower.tri(R)]
  }
  
  #construct G correlation and covariance matrices
  Gcor2 = list()
  Gcov = list()
  for(c in 1:nrow(X)){
    mat = diag(1, 2)
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
  
  #cross-environment correlations
  prho = 0.8
  cmat = matrix(0, nrow = 20, ncol = 20)
  for (i in 1:20) {
    for (j in i:20) {
      if (i %% 2 == 1 && j %% 2 == 1) {  # both odd
        cmat[i, j] = 0.8
        cmat[j, i] = 0.8
      } else if (i %% 2 == 0 && j %% 2 == 0) {  # both even
        cmat[i, j] = 0.8
        cmat[j, i] = 0.8
      }
    }
  }
  diag(cmat) = 1
  G_z = rmvnorm(N, mean = rep(0, 20), sigma = cmat)
  G_zs = lapply(1:10, function(i) {
    G_z[, ((2*i - 1):(2*i))]})
  
  G_s = list()
  for(i in 1:length(Gcor2)){
    G_s[[i]] = LA %*% G_zs[[i]] %*% chol(Gcov[[i]])
  }
  
  df$id = rep(1:N, each = 2)
  df$a1 = 0
  df$a2 = 0
   for(i in 1:nrow(df)){
    df[i,c("a1","a2")] = G_s[[df$c_id[i]]][df$id[i],]
  }
  
  #generate residuals
  E = rlkjcorr(1, 2, 10)
  sd_E = rep(sqrt(0.7),2)
  E_cv = diag(sd_E) %*% E %*% diag(sd_E)
  res = rmvnorm(N*2, mean=rep(0,2), sigma=E_cv)
  
  #generate responses
  z = mu[df$c_id,] + df[,c(paste0("a",1:2))] + res
  colnames(z) = paste0("z",1:2)
  df = cbind(df, z)
  
  #matrix indexing subjects in contexts
  c_list = list()
  for(c in 1:nrow(X)){
    c_list[[c]] = unique(df[df$c_id==c,"id"])
  }
  cmat = matrix(t(sapply(c_list, function(x) x[1:max(lengths(c_list))])),
                ncol = max(lengths(c_list)))
  cmat = t(apply(cmat, 1, function(x) `length<-`(na.omit(x), length(x))))
  cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) ) 
  cmat[is.na(cmat)] = 0 
  
  #new index corresponding to subject position in cmat
  temp = t(cmat)
  corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))
  df$idc = match(paste(df$id, df$c_id,sep="."), paste(corder$id,corder$c,sep="."))
  
  #character state values
  df$ebin = ifelse(df$season=="winter",0,1)
  dflow = df[df$ebin==0,]
  dfhigh = df[df$ebin==1,]
  
  #stan data list
  list(
    variables=list(B_m = B_m[2,], B_v = B_v[2,], B_cpc = B_cpc[2], 
                   sd_E = sd_E, L_E = t(chol(E)), Z_G = G_z,
                   G_vdt = G_vdt, G_cdt = G_cdt, G_cvdt = G_cvdt, 
                   es = c(es[1], es[2], es[2]*es[1]) ),
    
    generated=list(N = nrow(df), C = nrow(X), I = max(df$id),
                   D = 2, P = ncol(X), sd_E = sd_E, 
                   id = df$id, c_id = df$c_id, idc = df$idc,
                   X = X, x = df$x2, Xpred = matrix(c(1,0,1,1), nrow = 2,ncol = 2),
                   A = A, cmat = cmat, cm = ncol(cmat),cn = cmat_n,
                   cnt = sum(cmat_n), z = df[,c(paste0("z",1:2))],
                   zcs = data.frame(dflow[,c(paste0("z",1:2))],
                                    dfhigh[,c(paste0("z",1:2))]),
                   acs = data.frame(dflow[,c(paste0("a",1:2))],
                                    dfhigh[,c(paste0("a",1:2))]))
    
  )
  
}

#function for predicting values from random regression
covar_func <- function(G){
  X<-matrix(c(
    1,0,0,0,
    0,0,1,0,
    1,1,0,0,
    0,0,1,1
  ),4,4,byrow=TRUE)
  return( X %*% G %*% t(X) )
}

#load models
CRN_mod = cmdstan_model(stan_file = "CRNsim_mod.stan", 
                        stanc_options = list("O1"))
cs_mod = cmdstan_model(stan_file = "cs_mod.stan", 
                       stanc_options = list("O1"))
rr_mod = cmdstan_model(stan_file = "rr_mod.stan", 
                       stanc_options = list("O1"))

#multiple iterations#########################################

#general settings
nsim = 50
N = c(300, 600)
es = list(c(0.3, 0.1), c(0.6, 0.2))
y = 5

#specific settings
#to replicate all results, rerun with changed settings
Np = N[[1]]
esp = es[[1]]
yp = y[[1]]

dls = list()
for(d in 1:nsim){ 
  dls[[d]] =  generator_function(N = Np, years = yp, es = esp)
}
predls = list()
ndiv = list()
for(d in 1:nsim){
  
#crn
est = CRN_mod$sample(
  data = dls[[d]]$generated, 
  iter_sampling = 500, 
  iter_warmup = 1000, 
  init = 0.01, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.9, 
  max_treedepth = 10,
  refresh = 10)

divg = data.frame(crn = est$diagnostic_summary()$num_divergent, cs = 0, rr = 0)

post = extract(est)
  v10 = exp(post$B_v[,1,1])
  v11 = exp(post$B_v[,1,1] + post$B_v[,2,1])
  v20 = exp(post$B_v[,1,2])
  v21 = exp(post$B_v[,1,2] + post$B_v[,2,2])
  c0 = tanh(post$B_cpc[,1,1])
  c1 = tanh(post$B_cpc[,1,1] + post$B_cpc[,2,1])

  crn.df = data.frame(
  v10 = log(v10),
  v11 = log(v11),
  v20 = log(v20),
  v21 = log(v21),
  c00 = atanh(c0),
  c11 = atanh(c1), 
  v1 = v11 - v10,
  v2 = v21 - v20,
  c1 = c1 - c0,
  cv = c1 * sqrt(v11) * sqrt(v21) - 
       c0 * sqrt(v10) * sqrt(v20),
  pv1 = sum(v11 - v10>0)/length(v11),
  pv2 = sum(v21 - v20>0)/length(v11 - v10),
  pc1 = sum(c1 - c0>0)/length(v11 - v10),
  pcv = sum(c1 * sqrt(v11) * sqrt(v21) - 
              c0 * sqrt(v10) * sqrt(v20) > 0)/length(v11 - v10),
  mod = "crn",
  itr = d)

crnl = reshape2::melt(crn.df, id.var = c("itr","mod"))

#cs
est2 = cs_mod$sample(
  data = dls[[d]]$generated, 
  iter_sampling = 500, 
  iter_warmup = 1000, 
  init = 0.01, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.9, 
  max_treedepth = 10,
  refresh = 10)

divg$cs = est2$diagnostic_summary()$num_divergent

post = extract(est2)
cs.df = data.frame(
  v10 = log(post$sd_G[,1]^2),
  v11 = log(post$sd_G[,3]^2),
  v20 = log(post$sd_G[,2]^2),
  v21 = log(post$sd_G[,4]^2),
  c00 = atanh(post$Gr[,1,2]),
  c11 = atanh(post$Gr[,3,4]),
  v1 = post$sd_G[,3]^2 - post$sd_G[,1]^2,
  v2 = post$sd_G[,4]^2 - post$sd_G[,2]^2,
  c1 = post$Gr[,3,4] - post$Gr[,1,2],
  cv = post$G[,3,4] - post$G[,1,2],
  pv1 = sum(post$sd_G[,3]^2 - post$sd_G[,1]^2>0)/length(v11 - v10),
  pv2 = sum(post$sd_G[,4]^2 - post$sd_G[,2]^2 > 0)/length(v11 - v10),
  pc1 = sum(post$Gr[,3,4] - post$Gr[,1,2]>0)/length(v11 - v10),
  pcv = sum(post$G[,3,4] - post$G[,1,2] > 0)/length(v11 - v10),
  mod = "cs",itr = d)
csl = reshape2::melt(cs.df, id.var = c("itr","mod"))

#rr
est3 = rr_mod$sample(
  data = dls[[d]]$generated, 
  iter_sampling = 500, 
  iter_warmup = 1000, 
  init = 0.01, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.8, 
  max_treedepth = 10,
  refresh = 10)

divg$rr = est3$diagnostic_summary()$num_divergent
ndiv[[d]] = divg

post = extract(est3)

v10 = list()
v11 = list()
v20 = list()
v21 = list()
c00 = list()
c11 = list()
v1 = list()
v2 = list()
c1 = list()
cv1 = list()
for(i in 1:nrow(post$u)){
  a = post$mat_G[i,,]
  xp = matrix(c(
    1,0,0,0,
    0,0,1,0,
    1,1,0,0,
    0,0,1,1
  ),4,4,byrow=TRUE)
  z = a %*% t(xp)
  v10[[i]]  = var(z[,1])
  v11[[i]]  = var(z[,3])
  v20[[i]]  = var(z[,2])
  v21[[i]]  = var(z[,4])
  c00[[i]] = cor(z[,1:2])[2,1]
  c11[[i]] = cor(z[,3:4])[2,1]
  v1[[i]]  = var(z[,3]) - var(z[,1])
  v2[[i]]  = var(z[,4]) - var(z[,2])
  c1[[i]] = cor(z[,3:4])[2,1] - cor(z[,1:2])[2,1]
  cv1[[i]] = cor(z[,3:4])[2,1] * sd(z[,3]) * sd(z[,4]) -
    cor(z[,1:2])[2,1] * sd(z[,1]) * sd(z[,2])
}
rr.df = data.frame(
  v10 = log(unlist(v10)),
  v11 = log(unlist(v11)),
  v20 = log(unlist(v20)),
  v21 = log(unlist(v21)),
  c00 = atanh(unlist(c00)),
  c11 = atanh(unlist(c11)),
  v1 = unlist(v1),
  v2 = unlist(v2),
  c1 = unlist(c1),
  cv = unlist(cv1),
  pv1 = sum(unlist(v1)>0)/nrow(post$u),
  pv2 = sum(unlist(v2)>0)/nrow(post$u),
  pc1 = sum(unlist(c1)>0)/nrow(post$u),
  pcv = sum(unlist(cv1)>0)/nrow(post$u),
  mod = "rr",itr = d)
rrl = reshape2::melt(rr.df, id.var = c("itr","mod"))

predl = rbind(crnl, csl, rrl)
predl$tv = ifelse(predl$variable=="v1", dls[[d]]$variables$G_vdt[[1]],
                  ifelse(predl$variable=="v2", dls[[d]]$variables$G_vdt[[2]],
                         ifelse(predl$variable=="c1", dls[[d]]$variables$G_cdt,
                                dls[[d]]$variables$G_cvdt)))

predl$es = ifelse(predl$variable=="v1", dls[[d]]$variables$es[[1]],
                  ifelse(predl$variable=="v2", dls[[d]]$variables$es[[1]],
                         ifelse(predl$variable=="c1", dls[[d]]$variables$es[[2]],
                                dls[[d]]$variables$es[[3]])))

predls[[d]] = predl
}
predlr = do.call(rbind, predls)
mytextp = paste0("predlr_N", Np, "_y", yp, "_es", paste(esp, collapse = "-"))
mytextd = paste0("ndiv_N", Np, "_y", yp, "_es", paste(esp, collapse = "-"))
saveRDS(predlr, paste0(mytextp, ".RDS"))
saveRDS(ndiv, paste0(mytextd, ".RDS"))

#results###########################################

library(ggplot2)
library(cowplot)

#load results
sN = readRDS("predlr_N300_y5_es0.3-0.1.RDS")
lN = readRDS("predlr_N600_y5_es0.3-0.1.RDS")
sN2 = readRDS("predlr_N300_y5_es0.6-0.2.RDS")
lN2 = readRDS("predlr_N600_y5_es0.6-0.2.RDS")

#combine
sN$N = "300"
sN$E = "e1"
lN$N = "600"
lN$E = "e1"

sN2$N = "300"
sN2$E = "e2"
lN2$N = "600"
lN2$E = "e2"

predlc = rbind(sN, lN, sN2, lN2)
predlc$bias = (predlc$tv - predlc$value) / (predlc$tv) 

#plot
bias =
ggplot(predlc[predlc$variable=="cv" & predlc$N %in% c("300","600") ,], 
       aes(x = N, y = bias, group = interaction(N,mod)))+
  scale_color_manual(values = c("#ff395d", "#aa2aae", "#541bff"))+
  geom_boxplot(aes(color = mod), position = position_dodge(0.85),
               outlier.shape = NA, outliers = F, fill = NA, size = 1)+
  geom_hline(aes(yintercept = tv), linetype = "dashed")+
  facet_wrap(.~ E, nrow = 1)+
  theme(    panel.background = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA)  )+
  guides(color = FALSE)

power =
ggplot(predlc[predlc$variable=="pcv" & predlc$N %in% c("300","600"),], 
       aes(x = N, y = value, group = interaction(N,mod)))+
  scale_color_manual(values = c("#ff395d", "#aa2aae", "#541bff"))+
  geom_boxplot(aes(color = mod), position = position_dodge(0.85),
               outlier.shape = NA, outliers = F, fill = NA, size = 1)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(.~ E, nrow = 1)+
  theme(    panel.background = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA)  )+
  guides(color = FALSE)

cbp = plot_grid(bias, power, align = "hv", nrow = 2)
save_plot("comp_mods.png", cbp, base_height = 7, base_width = 10, dpi = 600)

#temp sim
png("temp_sim.png", width = 10, height = 3, units = "in", res = 600)
#first two panels
{  
  par(mfrow = c(1,3), mar = c(2, 3, 2, 3), oma = c(1, 1, 1, 1))
  
  hx = 1 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  lx = 0 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  x1 = scale(c(hx,lx))
  p1 = acf(x1, plot = F)
  
  hx = 1 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  lx = 0 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  x2 = scale(c(hx,lx))
  p2 = acf(x2, plot = F)
  
  hx = 1 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  lx = 0 + arima.sim(model = list(ar = c(0.9, -0.5), ma = 0.9), mean = 0, sd = 0.1,
                     n = 5) + 0.1*1:5
  x3 = scale(c(hx,lx))
  p3 = acf(x3, plot = F)
  
  #temps
  plot(x1[1:5], col = "deeppink4", ylim = c(-2,2), type = "l", lwd = 2,
       xlab = "", ylab = "")
  lines(x1[6:10], col = "cornflowerblue", lwd = 2)
  lines(x2[1:5], col = "deeppink4", lty = 2, lwd = 2)
  lines(x2[6:10], col = "cornflowerblue", lty = 2, lwd = 2)
  lines(x3[1:5], col = "deeppink4", lty = 3, lwd = 2)
  lines(x3[6:10], col = "cornflowerblue", lty = 3, lwd = 2)
  
  #acf
  lag = 0:9
  lag2 = 0:9 - 0.2
  lag3 = 0:9 + 0.2
  cols = ifelse((lag %% 2) == 1, "cornflowerblue", "deeppink4")
  plot(p1$lag, p1$acf, type = "n", lty = 1, ylim = c(-0.5, 1),
       xaxt = "n", xlab = "", ylab = "ACF")
  for (i in 1:10) {
    lines(rep(lag[i], 2), c(0, p1$acf[i]), col = cols[i], lwd = 2)
    lines(rep(lag2[i], 2), c(0, p2$acf[i]), col = cols[i], lwd = 2, lty = 2)
    lines(rep(lag3[i], 2), c(0, p3$acf[i]), col = cols[i], lwd = 2, lty = 3)
    axis_lags = seq(0, 8, by = 2)
    axis_labels = 1:5
    axis(1, at = axis_lags, labels = axis_labels)
  }
 } 
#third panel
{
  set.seed(123)
  library(MASS)
  
  draw_ellipse_with_axes = function(center, covmat, radius = 1, col = "black", lty = 1, lwd = 2) {
    # Ellipse outline
    angles = seq(0, 2 * pi, length.out = 100)
    unit_circle = cbind(cos(angles), sin(angles))
    ell = radius * unit_circle %*% chol(covmat) + matrix(center, 100, 2, byrow = TRUE)
    lines(ell, col = col, lwd = lwd, lty = lty)
    
    # Principal axes
    eig = eigen(covmat)
    for (i in 1:2) {
      vec = eig$vectors[, i] * sqrt(eig$values[i]) * radius
      segments(center[1] - vec[1], center[2] - vec[2],
               center[1] + vec[1], center[2] + vec[2],
               col = adjustcolor(col, 0.8), lwd = 1.5, lty = lty)
    }
  }
  
  #Covariance matrices
  mu = c(0, 0)
  sigma1 = matrix(c(0.3, 0, 0, 0.3), 2)
  sigma2 = matrix(c(0.41, 0.04, 0.04, 0.41), 2)
  sigma3 = matrix(c(0.55, 0.11, 0.11, 0.55), 2)
  
  #par(mfrow = c(1, 2), mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1))
  plot(NA, xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Z1", ylab = "Z2")
  draw_ellipse_with_axes(mu, sigma1, col = "black", lty = 1)
#  plot(NA, xlim = c(-1, 1), ylim = c(-1, 1),
#       xlab = "Z1", ylab = "Z2")
  draw_ellipse_with_axes(mu, sigma2, col = "green3", lty = 2)
  draw_ellipse_with_axes(mu, sigma3, col = "seagreen", lty = 2)
}
dev.off()
