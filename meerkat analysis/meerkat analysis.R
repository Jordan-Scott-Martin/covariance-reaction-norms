#prep workspace#############################

#load data
setwd("...")

library(cmdstanr)
library(rstan)
library(shinystan)
library(reshape2)
library(ggplot2)
library(tidybayes)
library(plyr)
library(cowplot)

set_cmdstan_path("...")

#data prep
{
#relatedness matrices
invAfd = as.matrix(readRDS("invA_feed2.rds"))
invAbs = as.matrix(readRDS("invA_BS2.rds"))
invAgd = as.matrix(readRDS("invA_GD2.rds"))
colnames(invAfd) = rownames(invAfd)
colnames(invAbs) = rownames(invAbs)
colnames(invAgd) = rownames(invAgd)

#subset to common subject matrix
subj = intersect(rownames(invAgd),
        intersect(rownames(invAfd), rownames(invAbs)))
invA = invAfd[subj,subj]
A = solve(invA) #undo matrix inversion
A = cov2cor(A)

#empirical data
fd = read.csv("df_feed.csv")
bs = read.csv("df_babysit.csv")
gd = read.csv("df_guard.csv")

#subset to common subjects
fd = fd[fd$focal.id %in% subj,]
bs = bs[bs$focal.id %in% subj,]
gd = gd[gd$focal.id %in% subj,]
A = A[rownames(A) %in% subj, colnames(A) %in% subj]

#adjust rownames
rownames(fd) = seq(1:nrow(fd))
rownames(bs) = seq(1:nrow(bs))
rownames(gd) = seq(1:nrow(gd))

#age months to years
fd$age = round(fd$age/12)
bs$age = round(bs$age/12)
gd$age = round(gd$age/12)

#group size round to units of 5
mround = function(x,base) base*round(x/base)
fd$groupsize = mround(fd$groupsize, 5)
bs$groupsize = mround(bs$groupsize, 5)
gd$groupsize = mround(gd$groupsize, 5)

#find unique environmental contexts across repeated measures
match.vars = c("groupsize","age","sex.m","dominant")
fdu = data.frame(unique(fd[,match.vars]))
bsu = data.frame(unique(bs[,match.vars]))
gdu = data.frame(unique(gd[,match.vars]))

library(dplyr)
contexts = inner_join(fdu, inner_join(bsu, gdu, 
                  relationship = "many-to-many"), 
                  relationship = "many-to-many")

#numeric indices of environmental contexts used for
#aligning context-specific G matrices with repeated measures
library(plyr)
fd.cor.id = list() 
bs.cor.id = list()
gd.cor.id = list()

for(i in 1:nrow(contexts)){
  fd.cor.id[[i]] = as.integer(rownames(match_df(fd, contexts[i,])))
  bs.cor.id[[i]] = as.integer(rownames(match_df(bs, contexts[i,])))
  gd.cor.id[[i]] = as.integer(rownames(match_df(gd, contexts[i,])))
}
fd.cor.id  = data.frame(context=rep(seq_along(fd.cor.id),lengths(fd.cor.id)),
             row=unlist(fd.cor.id))
bs.cor.id  = data.frame(context=rep(seq_along(bs.cor.id),lengths(bs.cor.id)),
             row=unlist(bs.cor.id))
gd.cor.id  = data.frame(context=rep(seq_along(gd.cor.id),lengths(gd.cor.id)),
             row=unlist(gd.cor.id))

#drop rows lacking common environmental contexts
fd = fd[rownames(fd) %in% fd.cor.id$row,]
bs = bs[rownames(bs) %in% bs.cor.id$row,]
gd = gd[rownames(gd) %in% gd.cor.id$row,]

fd$context.id = fd.cor.id[order(fd.cor.id$row),"context"]
bs$context.id = bs.cor.id[order(bs.cor.id$row),"context"]
gd$context.id = gd.cor.id[order(gd.cor.id$row),"context"]

#create matrix of context-specific IDs
c_list = list() #which IDs have all variables in which contexts
for(c in 1:max(fd$context.id)){
  c_list[[c]] = unique(c(fd[fd$context.id==c, "focal.id"],
                         bs[bs$context.id==c, "focal.id"],
                         gd[gd$context.id==c, "focal.id"]))}
c_mat = matrix(t(sapply(c_list, function(x) x[1:max(lengths(c_list))])),
               ncol = max(lengths(c_list)))

#create numeric indices of random effects for Stan
key.focal = unique(subj)
new.fid = seq(1:length(key.focal))
focal_fd = new.fid[match(fd$focal.id, key.focal)]
focal_bs = new.fid[match(bs$focal.id, key.focal)]
focal_gd = new.fid[match(gd$focal.id, key.focal)]
id_fun = function(x){
  key.group = unique(x$group.id)
  key.season = unique(x$season.id)
  key.focalxseason = unique(x$focal.by.season)

  new.gid = seq(1:length(key.group))
  new.sid = seq(1:length(key.season))
  new.fxsid = seq(1:length(key.focalxseason))
  
  return(data.frame(
    group = new.gid[match(x$group.id, key.group)],
    season = new.sid[match(x$season.id, key.season)],
    fxs = new.fxsid[match(x$focal.by.season, key.focalxseason)])
    ) 
  }  

key.c = rownames(contexts)
new.cid = seq(1:nrow(contexts))
c_fd = new.cid[match(fd$context.id, key.c)]
c_bs = new.cid[match(bs$context.id, key.c)]
c_gd = new.cid[match(gd$context.id, key.c)]
rownames(contexts) = new.cid

fd = cbind(fd, focal_fd, id_fun(fd), c_fd)
bs = cbind(bs, focal_bs, id_fun(bs), c_bs)
gd = cbind(gd, focal_gd, id_fun(gd), c_gd)

#sort relatedness matrix to match ID order
Ab = A
dimnames(Ab)[[1]] = new.fid[match(dimnames(Ab)[[1]], key.focal)]
dimnames(Ab)[[2]] = new.fid[match(dimnames(Ab)[[2]], key.focal)]
Ab = as.matrix(Ab[order(as.numeric(row.names(Ab))), order(as.numeric(colnames(Ab)))])

#adjust c_matrix indices
cmat.id = matrix(new.fid[match(c_mat, key.focal)], nrow = nrow(c_mat), ncol = ncol(c_mat))
cmat.id = t(apply(cmat.id, 1, function(x) `length<-`(na.omit(x), length(x))))
cmat.n = apply(cmat.id, 1, FUN = function(x) sum(!is.na(x)) )
cmat.id[is.na(cmat.id)] = 0
cnt = sum(cmat.n)

#vectors indexing individual_context

#fd$cid_fid = paste(fd$focal_fd, fd$c_fd,sep="_")
#bs$cid_fid = paste(bs$focal_bs, bs$c_bs,sep="_")
#gd$cid_fid = paste(gd$focal_gd, gd$c_gd,sep="_")

temp = t(cmat.id)
corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat.id)), times = cmat.n))
fd$cid_fid = match(paste(fd$focal_fd, fd$c_fd,sep="_"), paste(corder$id,corder$c,sep="_"))
bs$cid_fid = match(paste(bs$focal_bs, bs$c_bs,sep="_"), paste(corder$id,corder$c,sep="_"))
gd$cid_fid = match(paste(gd$focal_gd, gd$c_gd,sep="_"), paste(corder$id,corder$c,sep="_"))

#standardize centered covariates for effect size comparison
contexts$age = (contexts$age - mean(contexts$age))/sd(contexts$age)
contexts$groupsize = (contexts$groupsize - mean(contexts$groupsize))/sd(contexts$groupsize)

fd$offset = fd$offset.min.cent/sd(fd$offset.min.cent)
gd$offset = gd$offset.min.cent/sd(gd$offset.min.cent)

#create matrix of environmental predictors (nrow = # unique contexts of measurement)
X = data.frame(model.matrix(rep(1,nrow(contexts)) ~ 
                              age*sex.m*dominant + + groupsize + I(groupsize^2), 
                              data = contexts))

#data list for Stan
stan.df = list(
  D = 3, #number of traits
  P = ncol(X),
  nref = 3, #number of non-focal ID random effects
  nfd = nrow(fd),
  nbs = nrow(bs),
  ngd = nrow(gd),
  C = nrow(contexts), #unique environmental contexts
  
  nfocal = nrow(Ab),
  focal_fd = fd$focal_fd,
  focal_bs = bs$focal_bs,
  focal_gd = gd$focal_gd,
  
  context_fd = fd$c_fd,
  context_bs = bs$c_bs,
  context_gd = gd$c_gd,
  
  cmat = cmat.id,
  cm = ncol(cmat.id),
  cnt = cnt,
  cn = cmat.n,
  cfid_fd = fd$cid_fid,
  cfid_bs = bs$cid_fid,
  cfid_gd = gd$cid_fid,
  
  ngroup_fd = length(unique(fd$group)),
  group_fd = fd$group,
  ngroup_bs = length(unique(bs$group)),
  group_bs = bs$group,
  ngroup_gd = length(unique(gd$group)),
  group_gd = gd$group,
  
  nseason_fd = length(unique(fd$season)),
  season_fd = fd$season,
  nseason_bs = length(unique(bs$season)),
  season_bs = bs$season,
  nseason_gd = length(unique(gd$season)),
  season_gd = gd$season,
  
  nfxs_fd = max(fd$fxs),
  fxs_fd = fd$fxs,
  nfxs_bs = max(bs$fxs),
  fxs_bs = bs$fxs,
  nfxs_gd = max(gd$fxs),
  fxs_gd = gd$fxs,
  
  A = Ab,
  X = X,
  
  offset_fd = fd$offset,
  offset_gd = gd$offset,
  
  feed = fd$feed.count,
  babysit = bs$babysit,
  bs_tot = bs$total.day,
  guard =  gd$guard.min 
  )
}
saveRDS(stan.df, "meerk_data.RDS")
stan.df = readRDS("meerk_data.RDS")

#estimate CRN model#############################

#compile model
mod = cmdstan_model(stan_file = "meerkat_CRN_mod.stan",
                    stanc_options = list("O1"))

est = mod$sample(
  data = stan.df,
  iter_sampling = 500,
  iter_warmup = 500,
  init = 0.01,
  chains = 6, 
  parallel_chains = 6,
  adapt_delta = 0.80,
  max_treedepth = 12,
  refresh = 10) 

saveRDS(est, "fit_meerkat_CRN.RDS")
fit = readRDS("fit_meerkat_CRN.RDS")
launch_shinystan(fit)

#predict environmental effects on (co)variances################

#custom functions for getting results from Stan model
source("CRN functions.R")

#name of traits in order of Stan matrices (1, 2, ... D)
traits = c("FD", "BS", "GD")

#load model and extract posterior samples

#load dataset used for analysis
stan.df = readRDS("meerk_data.RDS")
fit = readRDS("fit_meerkat_CRN.RDS")
post = extract(fit)

#summarize beta coefficients for variances###############

v_b = v_beta_f(x = stan.df$X, traits = traits, v_b = post$B_v)

#posterior median
aggregate(.~ element, data = v_b, FUN = function(x) 
  round(median(x),2))

#posterior probability in direction of posterior median
aggregate(.~ element, data = v_b, FUN = function(x) 
  sum(sign(median(x))==sign(x))/length(x))

#posterior 90% CIs
aggregate(.~ element, data = v_b, FUN = function(x) 
  round(quantile(x, c(0.05,0.95)),2), simplify = F)

#summarize beta coefficients for genetic correlations####

#generate beta coefficients for genetic correlations
cor_b = cor_beta_f(x = stan.df$X, traits = traits, cpc_b = post$B_cpc)

#posterior median
aggregate(.~ element, data = cor_b, FUN = function(x) 
  round(median(x),2))

#posterior probability in direction of posterior median
aggregate(.~ element, data = cor_b, FUN = function(x) 
  sum(sign(median(x))==sign(x))/length(x))

#posterior 90% CIs
aggregate(.~ element, data = cor_b, FUN = function(x) 
  round(quantile(x, c(0.05,0.95)),2), simplify = F)


#plot posteriors of beta coefficients####################
v_bl = melt(v_b)
v_bl2 = v_bl[v_bl$variable!="X.Intercept.", ] #excl intercepts
cor_bl = melt(cor_b)
cor_bl2 = cor_bl[cor_bl$variable!="X.Intercept.", ] #excl intercepts
cor_bl2$element = factor(cor_bl2$element, levels = c("FD_BS", "FD_GD", "BS_GD"))

v_bp = 
  ggplot(v_bl2, aes(x = value, color = variable, fill = variable))+
  geom_density(aes(y = after_stat(scaled)), alpha = 0.30)+
  geom_vline(xintercept = 0, lty = "dotted")+
  facet_wrap(.~element)+
  labs(x=bquote(bold(paste("log effect size ", sigma[a]^2))))+
  theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
    legend.position = "top",
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_text(size=12,face="bold"),
    axis.title.y = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(size = 1),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.background= element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
    guides()

cor_bp = 
  ggplot(cor_bl2, aes(x = value, color = variable, fill = variable))+
  geom_density(aes(y = after_stat(scaled)), alpha = 0.30)+
  geom_vline(xintercept = 0, lty = "dotted")+
  facet_wrap(.~element)+
  labs(x=bquote(bold(paste("atanh effect size ", r[a]))))+
  theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
    legend.position = "top",
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_text(size=12,face="bold"),
    axis.title.y = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(size = 1),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.background= element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
    guides()

library(cowplot)
cbp = plot_grid(v_bp, cor_bp, ncol = 1, rel_heights = c(1,1))
save_plot("fig 3a.png", cbp, base_height = 5, base_width = 8)


##############################################################
#plot covariance function for group size

#predictor matrix for generating predictions from model
#X_pred columns must follow order of original design matrix
#colnames(stan.df$X), but the colnames of X_pred are arbitrary

#sex x dom across ages
{
seqp = seq(-1,1, by = 0.2) #standardized values from -2 to +2
#agev = rep(seqp, 4)
#sexv = rep(c(0,0,1,1),each=length(seqp))
#domv = rep(c(0,1,0,1),each=length(seqp))
agev = rep(seqp, 4)
sexv = rep(c(0,0,1,1),each=length(seqp))
domv = rep(c(0,1,0,1),each=length(seqp))

X_pred = data.frame(int = 1, age = agev,
                    sex = sexv, dominant = domv,
                    groupsize = 0, I.groupsize.2. = 0,
                    age.sex.m = agev * sexv,
                    age.dominant = agev * domv, 
                    sex.m.dominant = sexv * domv,
                    age.sex.m.dominant = agev * sexv * domv)

mv_var = v_f(x = X_pred, traits = traits, v_b = post$B_v)
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits,  cpc = mv_cpc)
mv_cov = cov_f(x = X_pred, traits = traits, v_b = post$B_v, cpc = mv_cpc)

mv_cor$element = factor(mv_cor$element, levels = c("FD_BS", "FD_GD", "BS_GD"))
mv_cov$element = factor(mv_cov$element, levels = c("FD_BS", "FD_GD", "BS_GD"))

levels=6
gs_cov1 = 
    ggplot(mv_cov, 
      aes(x = age, y = value, color = as.factor(sex), fill = as.factor(sex)))+ 
      coord_cartesian(ylim = c(-1,1))+
  scale_x_continuous(expand=c(0,0), labels = c("-1","","0","","1"))+
  scale_y_continuous(expand=c(0,0), labels = c("-1","","0","","1"))+
      geom_hline(yintercept = 0, lty = "dashed")+
      geom_vline(xintercept = 0, lty = "dashed")+
      stat_lineribbon(size=2, .width = ppoints(levels), alpha=0.5/levels)+
      scale_color_manual(values = c("orange2","blue4"))+
      scale_fill_manual(values = c("orange2","blue4"))+
      labs(y=bquote(bold(paste(sigma[a]))), x = "age (standardized)")+
      facet_wrap(.~ dominant * element)+
      theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.line = element_line(size = 1),
        panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                  linetype="solid"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.spacing.x = unit(2, "lines"))+
        guides(fill="none", color ="none")
}

#group size
{
  seqp = seq(-1,1, by = 0.2) #standardized values from -2 to +2
  #gsv = rep(seqp, 4)
  #sexv = rep(c(0,1,0,1),each=length(seqp))
  #agev = rep(c(0,0,1,1),each=length(seqp))
  
  #gsv = rep(seqp, 4*2)
  #agev = rep(c(-1,1,-1,1,-1,1,-1,1),each=length(seqp))
  #sexv = rep(c(0,0,1,1),each=length(2*seqp))
  #domv = rep(c(0,1,0,1),each=length(2*seqp))
  
  #X_pred = data.frame(int = 1, age = agev, sex = sexv, dominant = 0, 
  #                    groupsize = gsv, I.groupsize.2. = gsv*gsv,
  #                    age.sex.m = agev*sexv,
  #                    age.dominant = agev*domv, 
  #                    sex.m.dominant = sexv*domv,
  #                    age.sex.m.dominant = agev*sexv*domv)
  
  X_pred = data.frame(int = 0, age = 0, sex = 0, dominant = 0, 
                      groupsize = seqp, I.groupsize.2. = seqp*seqp,
                      age.sex.m = 0,
                      age.dominant = 0, 
                      sex.m.dominant = 0,
                      age.sex.m.dominant = 0)
  
  mv_var = v_f(x = X_pred, traits = traits, v_b = post$B_v)
  mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
  mv_cor = cor_f(x = X_pred, traits = traits,  cpc = mv_cpc)
  mv_cov = cov_f(x = X_pred, traits = traits, v_b = post$B_v, cpc = mv_cpc)
  
  mv_cor$element = factor(mv_cor$element, levels = c("FD_BS", "FD_GD", "BS_GD"))
  mv_cov$element = factor(mv_cov$element, levels = c("FD_BS", "FD_GD", "BS_GD"))
  
  
  levels=6
  color = fill = "orchid"
  
  gs_cov2 = 
    ggplot(mv_cov, 
           aes(x = groupsize, y = value))+ 
    coord_cartesian(ylim = c(-0.3,0.3))+
    scale_x_continuous(expand=c(0,0), labels = c("-1","","0","","1"))+
    #scale_y_continuous(expand=c(0,0), labels = c("-0.2","","0","","0.5"))+
    geom_hline(yintercept = 0, lty = "dashed")+
    geom_vline(xintercept = 0, lty = "dashed")+
    stat_lineribbon(size=2, .width = ppoints(levels), alpha=0.5/levels,
                    fill = fill, color = color)+
    labs(y=bquote(bold(paste(sigma[a]))), x = "group size (standardized)")+
    facet_wrap(.~ element)+
    theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_text(size=12,face="bold"),
          axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.line = element_line(size = 1),
          panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
          panel.spacing.x = unit(2, "lines"))+
          guides(fill="none", color ="none")
}
    
cpp = plot_grid(gs_cov1, gs_cov2, ncol = 1, rel_heights = c(0.66,0.34))
save_plot("fig 4bc.png", cpp, base_height = 7, base_width = 6)


