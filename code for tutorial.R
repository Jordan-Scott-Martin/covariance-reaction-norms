#code for running analysis in tutorial

setwd("...")

source("CRN functions.R")

#simulate data
stan.dl = sim_CRN_QG(N = 600, Nc = 30, npred = 3, ntrait = 3, l_es = 0.5, u_es = 0.9)

library(shinystan)
library(cmdstanr)
library(rstan)

#directory for cmdstan installation
set_cmdstan_path("...")

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
  adapt_delta = 0.9, 
  max_treedepth = 10,
  refresh = 10)

saveRDS(est, "CRN_fit.RDS")
CRN_fit = readRDS("CRN_fit.RDS")
launch_shinystan(est)

#name of traits in order of stan.dl$z columns (1, 2, ..., D)
#(e.g. "count_egg", "limb.length", "BD1", etc.)
traits = paste0("Z",seq(1:stan.dl$D)) #z1-z3

#get posterior MCMC samples
post = extract(CRN_fit)

#get posteriors of RN parameters in long format
v_b = v_beta_f(x = stan.dl$X, traits = traits, v_b = post$B_v)
head(v_b, 3)

#posterior median
aggregate(.~ element, data = v_b, FUN = function(x) 
  round(median(x),2))

#posterior probability in direction of posterior median
aggregate(.~ element, data = v_b, FUN = function(x) 
  sum(sign(median(x))==sign(x))/length(x))

#posterior 90% CIs
aggregate(.~ element, data = v_b, FUN = function(x) 
  round(quantile(x, c(0.05,0.95)),2), simplify = F)

#generate beta coefficients for genetic correlations
cor_b = cor_beta_f(x = stan.dl$X, traits = traits, cpc_b = post$B_cpc)
head(cor_b, 3)

#posterior median
aggregate(.~ element, data = cor_b, FUN = function(x) 
  round(median(x),2))

#posterior probability in direction of posterior median
aggregate(.~ element, data = cor_b, FUN = function(x) 
  sum(sign(median(x))==sign(x))/length(x))

#posterior 90% CIs
aggregate(.~ element, data = cor_b, FUN = function(x) 
  round(quantile(x, c(0.05,0.95)),2), simplify = F)

#long format
v_bl = melt(v_b)
v_bl2 = v_bl[v_bl$variable!="X.Intercept.", ] #excl intercepts

#plot posteriors of variance RNs
ggplot(v_bl2, aes(x = value, color = variable, fill = variable))+
  geom_density(aes(y = after_stat(scaled)), alpha = 0.30)+
  geom_vline(xintercept = 0, lty = "dotted")+
  facet_wrap(.~element)+
  labs(x=bquote(bold(paste("log effect size ", sigma[a]^2))))+
  theme(
    legend.position = "top",
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size=12,face="bold"),
    axis.title.y = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_blank(),
    axis.line.x = element_line(linewidth = 1),
    axis.line.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.background= element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides()


cor_bl = melt(cor_b)
cor_bl2 = cor_bl[cor_bl$variable!="X.Intercept.", ] #excl intercepts

#plot posteriors of variance RNs
ggplot(cor_bl2, aes(x = value, color = variable, fill = variable))+
  geom_density(aes(y = after_stat(scaled)), alpha = 0.30)+
  geom_vline(xintercept = 0, lty = "dotted")+
  facet_wrap(.~element)+
  labs(x=bquote(bold(paste("atanh effect size ", r[a]))))+
  theme(
    legend.position = "top",
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size=12,face="bold"),
    axis.title.y = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_blank(),
    axis.line.x = element_line(linewidth = 1),
    axis.line.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.background= element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides()

#predictor matrix for generating predictions from model
#X_pred columns must follow order of original design matrix
#colnames(stan.df$X), but the colnames of X_pred are arbitrary
seq = seq(-1,1, by = 0.3) #standardized values from -1 to +1
X_pred = data.frame(int = 1, X1 = seq, X2 = 0, X3 = 0)

#functions for predicting genetic correlations and covariances
mv_cpc = cpc_f(x = X_pred, cpc_b = post$B_cpc)
mv_cor = cor_f(x = X_pred, traits = traits,  cpc = mv_cpc)
mv_cov = cov_f(x = X_pred, traits = traits, v_b = post$B_v, cpc = mv_cpc)

#visualize change in covariance
levels = 6 #number of bands for Bayesian CIs, try ppoints(levels)
ggplot(mv_cov, 
       aes(x = X1, y = value))+ 
  geom_hline(yintercept = 0, lty = "dashed")+
  geom_vline(xintercept = 0, lty = "dashed")+
  stat_lineribbon(size=2, .width = ppoints(levels), alpha=0.8/levels)+
  labs(x=bquote(bold(paste(X[1]))),
       y=bquote(bold(paste(sigma[a]))))+
  facet_wrap(.~ element)+
  theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.line = element_line(linewidth = 1),
        
        panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                  linetype="solid"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides(fill="none", color ="none")


