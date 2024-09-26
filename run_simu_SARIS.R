# ====================================================
# Script Name: utils_SARIS.R
# Description: Simulations of the article https://arxiv.org/abs/2408.13022
# Author: Tom Guédon tom.guedondurieu@gmail.com
# Date: 2024-09
# Version: 1.0
#
# R Version: R version 4.3.3 (2024-02-29 ucrt)
# Required Packages: (Liste des bibliothèques nécessaires)
#   - dplyr (1.1.4)
#   - ggplot2 (3.5.0)
#   - reshape2 (1.4.4)
#   - gridExtra (2.3)
#
# Usage:
# This script shall be executed linearly
# ====================================================





rm(list=ls())
# set working directory to source file location
source(file = "utils_SARIS.R")
library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(dplyr)
library(reshape2)

set.seed(0)

path_init = getwd()
newpath = "figures_new"
dir.create(newpath)
setwd(newpath)

##############################################-
##############################################-


# DEBUT CALCUL mu=1 sigma=1

##############################################-
##############################################-



rstar = 1

mu = 1
sigma = 1

f0 = function(x){return(rstar*dnorm(x))}
f1 = function(x){return(dnorm(x,mean = mu ,sd = sigma))}

r0 =0
K = 5000
Kheat = 300
M = 2*Kheat +2*K 

stepsize = function(x){if (x<K%/%2){return(0.1)}else{return(1/(1+x**0.66))}}
nexp = 50

epsilon_bridge = 10**(-8)

bridge_res = rep(0,nexp)

ris_bridge_res = rep(0,nexp)

saris_bridge_res = rep(0,nexp)

saris_bridgerm_res = rep(0,nexp)

saris_opt_res = rep(0,nexp)


for (n in 1:nexp){
  
  sample0 = metropolis(f0, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
  sample1 = metropolis(f1, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
  
  
  sample_mixte = rep(0,2*K)
  p0 = 1
  p1 = 1
  for (k in 1:(2*K)){
    
    if (runif(1)<0.5){
      
      sample_mixte[k]=sample0$samples[(Kheat+1):(Kheat+2*K)][p0]
      p0 = p0+1
      
    }
    else{
      
      sample_mixte[k]=sample1$samples[(Kheat+1):(Kheat+2*K)][p1]
      p1=p1+1
      
    }
  }
  
  rbridge=r0
  
  f0_sample1 = f0(sample1$samples[(Kheat+1):(Kheat+K)])
  f0_sample0 = f0(sample0$samples[(Kheat+1):(Kheat+K)])
  f1_sample0 = f1(sample0$samples[(Kheat+1):(Kheat+K)])
  f1_sample1 = f1(sample1$samples[(Kheat+1):(Kheat+K)])
  
  
delta = 10
while(delta>epsilon_bridge){
    delta = rbridge
    numerateur = f0_sample1/(f0_sample1 + rbridge * f1_sample1)
    denominateur = f1_sample0/(f0_sample0 + rbridge * f1_sample0)
    rbridge= sum(numerateur)/sum(denominateur)
    delta = abs(delta-rbridge)
  }
  
  bridge_res[n] = log(rbridge)
  
  f0_sample = f0(sample_mixte)
  f1_sample = f1(sample_mixte)
  
  S = function(r){return(sum(f0_sample/(f0_sample + r*f1_sample)) - sum((r*f1_sample/(f0_sample + r*f1_sample)) ))}
  
  ris_bridge = bissection(f = S, a = 0.0001 , b=100, tol = epsilon_bridge)
  
  ris_bridge_res[n] = log(ris_bridge)
  
  SARISlog_bridge = BridgeSARISlog(sample_mixte,f0,f1,r0 = 0,stepsize = stepsize,epsilon = 0)
  
  saris_bridgelog = SARISlog_bridge$r
  
  L = length(SARISlog_bridge$r)
  
  saris_bridge_res[n] = mean(saris_bridgelog[(L-K):(L)])
  
  SARIS_optlog = SARIS_RMlog(M, f0,f1, r0=0,  alpha = function(x,r){return(abs(f0(x)-r*f1(x)))},stepsize = stepsize,acc = 0.44)
  
  saris_optlog = SARIS_optlog$r
  
  L2 = length(saris_optlog)
  
  saris_opt_res[n] = mean(saris_optlog[(L2-K):L2])
  
  SARIS_bridgerm = SARIS_RMlog(M, f0,f1, r0=0,  alpha = function(x,r){return(f0(x)+r*f1(x))},stepsize = stepsize,acc = 0.44)
  
  saris_bridgerm = SARIS_bridgerm$r
  
  saris_bridgerm_res[n] = mean(saris_optlog[(L2-K):L2])
  
  print(n)
}



##############################################-
##############################################-

# FIN CALCUL mu=1 sigma=1


# DEBUT CALCUL mu=5 sigma=1

##############################################-
##############################################-
set.seed(0)
rstar = 1

mu = 5
sigma = 1

f0 = function(x){return(rstar*dnorm(x))}
f1 = function(x){return(dnorm(x,mean = mu ,sd = sigma))}

r0 =0
K = 5000
Kheat = 300
M = 2*Kheat +2*K 

stepsize = function(x){if (x<K%/%2){return(0.1)}else{return(1/(1+x**0.66))}}
nexp = 50

epsilon_bridge = 10**(-8)

bridge_res5 = rep(0,nexp)

ris_bridge_res5 = rep(0,nexp)

saris_bridge_res5 = rep(0,nexp)

saris_bridgerm_res5 = rep(0,nexp)

saris_opt_res5 = rep(0,nexp)


for (n in 1:nexp){
  
  sample0 = metropolis(f0, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
  sample1 = metropolis(f1, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
  
  
  sample_mixte = rep(0,2*K)
  p0 = 1
  p1 = 1
  for (k in 1:(2*K)){
    
    if (runif(1)<0.5){
      
      sample_mixte[k]=sample0$samples[(Kheat+1):(Kheat+2*K)][p0]
      p0 = p0+1
      
    }
    else{
      
      sample_mixte[k]=sample1$samples[(Kheat+1):(Kheat+2*K)][p1]
      p1=p1+1
      
    }
  }
  
  rbridge=r0
  
  f0_sample1 = f0(sample1$samples[(Kheat+1):(Kheat+K)])
  f0_sample0 = f0(sample0$samples[(Kheat+1):(Kheat+K)])
  f1_sample0 = f1(sample0$samples[(Kheat+1):(Kheat+K)])
  f1_sample1 = f1(sample1$samples[(Kheat+1):(Kheat+K)])
  
  
  delta = 10
  while(delta>epsilon_bridge){
    delta = rbridge
    numerateur = f0_sample1/(f0_sample1 + rbridge * f1_sample1)
    denominateur = f1_sample0/(f0_sample0 + rbridge * f1_sample0)
    rbridge= sum(numerateur)/sum(denominateur)
    delta = abs(delta-rbridge)
  }
  
  bridge_res5[n] = log(rbridge)
  
  f0_sample = f0(sample_mixte)
  f1_sample = f1(sample_mixte)
  
  S = function(r){return(sum(f0_sample/(f0_sample + r*f1_sample)) - sum((r*f1_sample/(f0_sample + r*f1_sample)) ))}
  
  ris_bridge = bissection(f = S, a = 0.0001 , b=100, tol = epsilon_bridge)
  
  ris_bridge_res5[n] = log(ris_bridge)
  
  SARISlog_bridge = BridgeSARISlog(sample_mixte,f0,f1,r0 = 0,stepsize = stepsize,epsilon = 0)
  
  saris_bridgelog = SARISlog_bridge$r
  
  L = length(SARISlog_bridge$r)
  
  saris_bridge_res5[n] = mean(saris_bridgelog[(L-K):(L)])
  
  SARIS_optlog = SARIS_RMlog(M, f0,f1, r0=0,  alpha = function(x,r){return(abs(f0(x)-r*f1(x)))},stepsize = stepsize,acc = 0.44)
  
  saris_optlog = SARIS_optlog$r
  
  L2 = length(saris_optlog)
  
  saris_opt_res5[n] = mean(saris_optlog[(L2-K):L2])
  
  SARIS_bridgerm = SARIS_RMlog(M, f0,f1, r0=0,  alpha = function(x,r){return(f0(x)+r*f1(x))},stepsize = stepsize,acc = 0.44)
  
  saris_bridgerm = SARIS_bridgerm$r
  
  saris_bridgerm_res5[n] = mean(saris_optlog[(L2-K):L2])
  
  print(n)
}


gros_plot2 = as.data.frame(list("r" = c(ris_bridge_res5,bridge_res5,saris_bridge_res5,saris_opt_res5, saris_bridgerm_res5),"method" = rep(c("RIS-MIXT","BRIDGE-OPT", "SARIS-MIXT", "SARIS-EXT-opt", "SARIS-EXT-mixt"),  rep(nexp,5))))
pl2 = ggplot(gros_plot2, aes(x=method, y=r))+
  geom_boxplot()+
  geom_hline(yintercept  = log(rstar))+
  scale_color_grey()+
  labs( x = "", y = "log(r)") +  # Ajouter des étiquettes d'axes
  
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 10,angle=0),
        axis.text.y = element_text(face = "bold", 
                                   size = 12),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=7)

pl2
limits<- ggplot_build(pl2)$layout$panel_params[[1]]$y.range
limits


gros_plot1 = as.data.frame(list("r" = c(ris_bridge_res,bridge_res,saris_bridge_res,saris_opt_res, saris_bridgerm_res),"method" = rep(c("RIS-MIXT","BRIDGE-OPT", "SARIS-MIXT", "SARIS-EXT-opt", "SARIS-EXT-mixt"),  rep(nexp,5))))

pl1 = ggplot(gros_plot1, aes(x=method, y=r))+
  geom_boxplot()+
  geom_hline(yintercept  = log(rstar))+
  scale_color_grey()+
  labs( x= "", y = "log(r)") +  # Ajouter des étiquettes d'axes
  
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 10,angle=0),
        axis.text.y = element_text(face = "bold", 
                                   size = 12),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"))+
  stat_summary(fun.y=mean, geom="point", shape=20, size=7) + ylim(limits) 

pl1


jpeg("p2.jpg")
print(pl2)
dev.off()

saveRDS(pl2,"pl2_plot.rds")


jpeg("pl1.jpg")
print(pl1)
dev.off()

saveRDS(pl1,"pl1_plot.rds")

##############################################-
##############################################-


# FIN CALCUL mu=5 sigma=1


# DEBUT CALCUL mu variable

##############################################-
##############################################-

epsilon_bridge = 10**(-8)
set.seed(0)
mu_vec = 1:10
rstar = 1
r0 =0
K = 5000
Kheat = 300
M = 2*Kheat +2*K

stepsize = function(x){if (x<K%/%2){return(0.1)}else{return(1/(1+x**0.66))}}
nexp = 50

sigma = 1


bridge_res = matrix(0,nrow = nexp,ncol = length(mu_vec))
saris_opt_res = matrix(0,nrow = nexp,ncol = length(mu_vec))

for (j in mu_vec){
  
  mu = mu_vec[j]
  
  f0 = function(x){return(rstar*dnorm(x))}
  f1 = function(x){return(dnorm(x,mean = mu ,sd = sigma))}
  
  
  
  
  
  for (n in 1:nexp){
    
    sample0 = metropolis(f0, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
    sample1 = metropolis(f1, sample_size = Kheat+2*K, burn_in = Kheat, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25)
    
    rbridge=r0
    
    f0_sample1 = f0(sample1$samples[(Kheat+1):(Kheat+K)])
    f0_sample0 = f0(sample0$samples[(Kheat+1):(Kheat+K)])
    f1_sample0 = f1(sample0$samples[(Kheat+1):(Kheat+K)])
    f1_sample1 = f1(sample1$samples[(Kheat+1):(Kheat+K)])
    

    delta = 10
    t=0
    while(delta>epsilon_bridge & t<200){
      delta = rbridge
      numerateur = f0_sample1/(f0_sample1 + rbridge * f1_sample1)
      denominateur = f1_sample0/(f0_sample0 + rbridge * f1_sample0)
      rbridge= sum(numerateur)/sum(denominateur)
      delta = abs(delta-rbridge)
      t=t+1
    }
    
    
    bridge_res[n,j] = log(rbridge)
    
    
    
    
    SARIS_geomlog = SARIS_RMlog(M, f0,f1, r0=r0,  alpha = function(x,r){return((sqrt(abs(f0(x)))+sqrt(abs(r*f1(x))))**2)},stepsize = stepsize,acc = 0.44)
    saris_geomlog = SARIS_geomlog$r
    len = length(saris_geomlog)
    
    saris_opt_res[n,j] = mean(saris_geomlog[(len-K):len])
    
    print(c(n,j))
    
    
  }
}
mean_geom = apply(X=saris_opt_res, FUN = mean,MARGIN =2)
sd_geom = apply(X=saris_opt_res, FUN = sd,MARGIN =2)

mean_bridge = apply(X=bridge_res, FUN = mean,MARGIN =2)
sd_bridge = apply(X=bridge_res, FUN = sd,MARGIN =2)

mean_tot = c(mean_bridge,mean_geom)
sd_tot = c(sd_bridge,sd_geom)
meth_tot = c(rep("bridge",length(mu_vec)),rep("saris", length(mu_vec)))

data_plot = data.frame(mean = mean_tot, sd = sd_tot, method = meth_tot, mu = as.factor(rep(mu_vec,2)))


pl3 <- ggplot(data_plot, aes(x = mu, y = mean, color = method)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5, position = position_dodge(width = 0.9), size = 1) +  # Ajouter les barres d'erreur avec un décalage
  geom_point(size = 3, position = position_dodge(width = 0.9)) +  # Ajouter les points avec un décalage
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "mu", y = "log(r)", color = "Estimator") +  # Ajouter des étiquettes d'axes
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))

# Afficher le graphique
print(pl3)
jpeg("pl3.jpg",width = 900,height=450)
print(pl3)
dev.off()
saveRDS(data_plot, "pl3_plot.rds")

##############################################-
##############################################-


# FIN CALCUL mu variable

# JOINT PROCEDURE

##############################################-
##############################################-

quantiles <- function(x) {
  c(q5 = quantile(x, 0.05), q95 = quantile(x, 0.95))
}

set.seed(0)


# joint procedure estimation et calculr lrt 

# vraies valeurs loi de X
mu_star = c(1,1)

gamma1_star = 1**2
gamma2_star = 1**2

Gamma_star = diag(c(gamma1_star,gamma2_star))

# bruit
sigma2_star = 2 # connu
sigma2 = 1
# paramètres du modèle

beta0_star = 0.1
beta_star = c(1,1)

# nombre d'obs

n = 200
n_mis = 20

# simulation des données

X_tot  = rmvnorm(n, mean = mu_star, sigma = Gamma_star)

residus = rnorm(n,mean = 0 , sd = sqrt(sigma2_star))

Y = beta0_star+ X_tot%*%beta_star + residus


X = X_tot

X[1:n_mis,2] = NA



####
# STOCKAGE
####

K=250

nexp = 20

# 
beta0_liste = rep(0,K)
beta_liste = matrix(0,nrow = K,ncol = 2)
gamma_liste = matrix(0,nrow =K,ncol = 2)
mu_liste = matrix(0,nrow = K, ncol =2)

beta0 = 0
beta = c(0,0)

mu = c(0,0)
gamma_vec = c(0.5,0.5)

beta0_liste[1] = beta0 
beta_liste[1,] = beta
gamma_liste[1,] = gamma_vec
mu_liste[1,] = mu
sample_liste = matrix(0,nrow =K-1,ncol = n_mis)
marglik_liste = rep(0,K-1)


# 
beta0_listeH0 = rep(0,K)
beta_listeH0 = matrix(0,nrow = K,ncol = 2)
gamma_listeH0 = matrix(0,nrow =K,ncol = 2)
mu_listeH0 = matrix(0,nrow = K, ncol =2)

beta0H0 = 0
betaH0 = c(0,0)

muH0 = c(0,0)
gamma_vecH0 = c(0.5,0.5)

beta0_listeH0[1] = beta0H0 
beta_listeH0[1,] = betaH0
gamma_listeH0[1,] = gamma_vecH0
mu_listeH0[1,] = muH0
sample_listeH0 = matrix(0,nrow =K-1,ncol = n_mis)
marglik_listeH0 = rep(0,K-1)

step.size = function(x){return(0.1)}
# 



lr_liste = matrix(0,nrow = K-1, ncol = nexp)
lr_exact = matrix(0,nrow = K-1, ncol = nexp)


for (t in 1:nexp){
  for (k in 1:(K-1)){
    
    X2_nouv = rep(0,n_mis)
    
    for (i in 1:n_mis){
      
      X2_nouv[i] = rnorm(1, mean =moy_post(Y[i],X[i,1],mu = mu[2],beta0 = beta0,beta1 = beta[1],beta2 = beta[2],gamma2 = gamma_vec[2],sigma2 = 1), sd = sqrt(var_post(beta = beta,gamma2 = gamma_vec[2],sigma2=1)) )
      
    }
    sample_liste[k,] = X2_nouv
    
    Xk = X
    Xk[1:n_mis,2] = X2_nouv
    
    g_beta0 = 0
    g_beta = 0
    g_mu = 0
    g_gamma=0
    
    for (i in 1:n){
      g_beta0 = g_beta0 + grad_beta0(Y[i], beta0, beta, Xk[i,], sigma2_star)/n
      g_beta = g_beta +grad_beta(Y[i], beta0, beta, Xk[i,], sigma2_star)/n
      
      g_mu = g_mu + grad_mu(Xk[i,],mu,gamma_vec)/n
      g_gamma = g_gamma + grad_gamma(Xk[i,],mu,gamma_vec)/n
      
    }
    
    beta0 = beta0 + step.size(k) * g_beta0
    beta = beta + step.size(k) * g_beta
    gamma_vec = gamma_vec + step.size(k) * g_gamma
    mu = mu + step.size(k) * g_mu
    
    
    marglik_liste[k] = marg_lik(Y, beta0, beta, X, sigma2, mu ,gamma_vec,n_mis)
    
    
    beta0_liste[k+1] = beta0 
    beta_liste[k+1,] = beta
    gamma_liste[k+1,] = gamma_vec
    mu_liste[k+1,] = mu
    
    ###
    # on refait pareil pour H0
    ###
    
    X2_nouv = rep(0,n_mis)
    
    # faire un vectoeur moyenne et une matrice variance et faire rmnvnorm
    
    for (i in 1:n_mis){
      
      X2_nouv[i] = rnorm(1, mean =moy_post(Y[i],X[i,1],mu = muH0[2],beta0 = beta0H0,beta1 = betaH0[1],beta2 = betaH0[2],gamma2 = gamma_vecH0[2],sigma2 = 1), sd = sqrt(var_post(beta = betaH0,gamma2 = gamma_vecH0[2],sigma2=1)) )
      
    }
    sample_listeH0[k,] = X2_nouv
    
    Xk = X
    Xk[1:n_mis,2] = X2_nouv
    
    g_beta0 = 0
    g_beta = 0
    g_mu = 0
    g_gamma=0
    
    for (i in 1:n){
      #g_beta0 = g_beta0 + grad_beta0(Y[i], beta0H0, betaH0, Xk[i,], sigma2_star)/n
      g_beta = g_beta +grad_beta(Y[i], beta0H0, betaH0, Xk[i,], sigma2_star)/n
      
      g_mu = g_mu + grad_mu(Xk[i,],muH0,gamma_vecH0)/n
      g_gamma = g_gamma + grad_gamma(Xk[i,],muH0,gamma_vecH0)/n
      
    }
    
    #beta0H0 = beta0H0 + step.size(k) * g_beta0
    betaH0 = betaH0 + step.size(k) * g_beta
    gamma_vecH0 = gamma_vecH0 + step.size(k) * g_gamma
    muH0 = muH0 + step.size(k) * g_mu
    
    marglik_listeH0[k] = marg_lik(Y, beta0=0, betaH0, X, sigma2, muH0 ,gamma_vecH0,n_mis)
    
    beta0_listeH0[k+1] = beta0H0 
    beta_listeH0[k+1,] = betaH0
    gamma_listeH0[k+1,] = gamma_vecH0
    mu_listeH0[k+1,] = muH0
    
    
    
    print(paste("t=",t,"k=",k))
  }
  
  
  
  
  # 
  
  

  sample_mixte = matrix(0,nrow = K-1,ncol = n_mis)
  
  for (k in 1:(K-1)){
    u = runif(1)
    if(u<0.5){sample_mixte[k,] = sample_listeH0[k,]
    }else{sample_mixte[k,] = sample_liste[k,]
    
    }
  }
  
  r_liste = matrix(0,nrow =K,ncol =n)
  
  
  logf0_sample = matrix(0,nrow =K-1,ncol = n)
  logf1_sample = matrix(0,nrow =K-1,ncol = n)
  
  
  
  for (k in 1:(K-1)){
    
    Xk = X
    Xk[1:n_mis,2] = sample_mixte[k,]
    
    for (i in 1:n){
      logf0_sample[k,i] = complete_llik(Y[i], beta0_listeH0[k], beta_listeH0[k,], Xk[i,], sigma2 = 1, mu_listeH0[k,] ,gamma_listeH0[k,] )
      logf1_sample[k,i] = complete_llik(Y[i], beta0_liste[k], beta_liste[k,], Xk[i,], sigma2 = 1, mu_liste[k,] ,gamma_liste[k,] )
    }
    
  }
  
  
  # une colonne c'est pour un ri 
  r_liste = matrix(0,nrow =K,ncol =n)
  step_size = function(x){return(0.1)}
  
  for (k in 1:(K-1)){
    for(i in 1:n){
      
      r_liste[k+1,i]= r_liste[k,i] + step_size(k)* (exp(logf0_sample[k,i]) - exp(logf1_sample[k,i]+r_liste[k,i]))/(exp(logf0_sample[k,i]) + exp(logf1_sample[k,i]+r_liste[k,i]))
      
    }
  }
  
  lr_liste[,t] = apply(-2*r_liste,FUN=sum,MARGIN=1)[2:K]
  lr_exact[,t] = -2*(marglik_listeH0- marglik_liste)
  
  
}

lr_liste_mean =colMeans(t(lr_liste))
lr_exact_mean = colMeans(t(lr_exact))
ite = rep(1:(K-1),nexp)



lr_sd = apply(lr_liste,FUN=sd,MARGIN=1)
lr_exact_sd = apply(lr_exact,FUN=sd,MARGIN=1)

#




diff_mse = (lr_liste-lr_exact)**2
diff_mse_sd = apply(diff_mse,FUN=sd,MARGIN=1)
diff_mse_mean = apply(diff_mse,FUN=mean,MARGIN=1)


quantiles_mse = t(apply(diff_mse, 1, quantiles))
quantiles_mse = as.data.frame(quantiles_mse)
colnames(quantiles_mse) = c("q5", "q95")


# Calculer les moyennes pour diff_mse
diff_mse_mean = apply(diff_mse, FUN = mean, MARGIN = 1)


ite_vec = c(50,100,150,200,249)

data_approx = data.frame("lr" = lr_liste_mean, "sd" =lr_sd, "method" = rep("saris", length(lr_liste_mean)),"iterations" = 1:length(lr_liste_mean))[ite_vec,]
data_exact = data.frame("lr" = lr_exact_mean, "sd" =lr_exact_sd, "method" = rep("exact", length(lr_exact_mean)), "iterations" = 1:length(lr_exact_mean))[ite_vec,]

data_joint = rbind(data_approx,data_exact)
data_mse = cbind(data.frame("mse" = diff_mse_mean,"iterations"= 1:length(diff_mse_mean)), quantiles_mse )[ite_vec,]

plot_joint <- ggplot(data_joint, aes(x = iterations, y = lr, color = method)) +
  geom_point(size = 5, position = position_dodge(10)) +  # Ajouter les points avec un d?calage
  geom_errorbar(aes(ymin = lr- sd, ymax = lr + sd),size=1, width = 10, position = position_dodge(width = 10)) +  # Ajouter les barres d'erreur avec le m?me d?calage
  labs(x = "iterations", y = "LRT") +  # Ajouter des ?tiquettes d'axes
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),) + ylim(-2, 2)
  plot_joint
  


plot_mse <- ggplot(data_mse, aes(x = iterations, y = mse)) +
  geom_point(size = 5, position = position_dodge(width = 0.9)) +  # Ajouter les points avec un d?calage
  geom_errorbar(aes(ymin = q5, ymax = q95),size=1, width = 10, position = position_dodge(width = 0.9)) +  # Ajouter les barres d'erreur avec le m?me d?calage
  labs(x = "iterations", y = "MSE") +  # Ajouter des ?tiquettes d'axes
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept=0) + ylim(-0.1, 0.5)



plot_mse

x11()
grid.arrange(plot_joint, plot_mse, ncol=2)
dev.off()

jpeg("lr_n10.jpg",width = 900,height=450)
plot_joint
dev.off()


jpeg("mse_n10.jpg",width = 900,height=450)
plot_mse
dev.off()

write.csv2(data.frame('lr' = lr_liste, "lr_exact" = lr_exact), file = 'lr_n10.csv')

##############################################-
##############################################-

# JOINT PROCEDURE

# missing 25 % of values

##############################################-
##############################################-


set.seed(0)


# joint procedure estimation et calculr lrt 

# vraies valeurs loi de X
mu_star = c(1,1)

gamma1_star = 1**2
gamma2_star = 1**2

Gamma_star = diag(c(gamma1_star,gamma2_star))

# bruit
sigma2_star = 2 # connu
sigma2 = 1
# paramètres du modèle

beta0_star = 0.1
beta_star = c(1,1)

# nombre d'obs

n = 200
n_mis = 50

# simulation des données

X_tot  = rmvnorm(n, mean = mu_star, sigma = Gamma_star)

residus = rnorm(n,mean = 0 , sd = sqrt(sigma2_star))

Y = beta0_star+ X_tot%*%beta_star + residus


X = X_tot

X[1:n_mis,2] = NA



####
# STOCKAGE
####

K=250

nexp = 20

# 
beta0_liste = rep(0,K)
beta_liste = matrix(0,nrow = K,ncol = 2)
gamma_liste = matrix(0,nrow =K,ncol = 2)
mu_liste = matrix(0,nrow = K, ncol =2)

beta0 = 0
beta = c(0,0)

mu = c(0,0)
gamma_vec = c(0.5,0.5)

beta0_liste[1] = beta0 
beta_liste[1,] = beta
gamma_liste[1,] = gamma_vec
mu_liste[1,] = mu
sample_liste = matrix(0,nrow =K-1,ncol = n_mis)
marglik_liste = rep(0,K-1)


# 
beta0_listeH0 = rep(0,K)
beta_listeH0 = matrix(0,nrow = K,ncol = 2)
gamma_listeH0 = matrix(0,nrow =K,ncol = 2)
mu_listeH0 = matrix(0,nrow = K, ncol =2)

beta0H0 = 0
betaH0 = c(0,0)

muH0 = c(0,0)
gamma_vecH0 = c(0.5,0.5)

beta0_listeH0[1] = beta0H0 
beta_listeH0[1,] = betaH0
gamma_listeH0[1,] = gamma_vecH0
mu_listeH0[1,] = muH0
sample_listeH0 = matrix(0,nrow =K-1,ncol = n_mis)
marglik_listeH0 = rep(0,K-1)

step.size = function(x){return(0.1)}
# 



lr_liste2 = matrix(0,nrow = K-1, ncol = nexp)
lr_exact2 = matrix(0,nrow = K-1, ncol = nexp)


for (t in 1:nexp){
  for (k in 1:(K-1)){
    
    X2_nouv = rep(0,n_mis)
    
    for (i in 1:n_mis){
      
      X2_nouv[i] = rnorm(1, mean =moy_post(Y[i],X[i,1],mu = mu[2],beta0 = beta0,beta1 = beta[1],beta2 = beta[2],gamma2 = gamma_vec[2],sigma2 = 1), sd = sqrt(var_post(beta = beta,gamma2 = gamma_vec[2],sigma2=1)) )
      
    }
    sample_liste[k,] = X2_nouv
    
    Xk = X
    Xk[1:n_mis,2] = X2_nouv
    
    g_beta0 = 0
    g_beta = 0
    g_mu = 0
    g_gamma=0
    
    for (i in 1:n){
      g_beta0 = g_beta0 + grad_beta0(Y[i], beta0, beta, Xk[i,], sigma2_star)/n
      g_beta = g_beta +grad_beta(Y[i], beta0, beta, Xk[i,], sigma2_star)/n
      
      g_mu = g_mu + grad_mu(Xk[i,],mu,gamma_vec)/n
      g_gamma = g_gamma + grad_gamma(Xk[i,],mu,gamma_vec)/n
      
    }
    
    beta0 = beta0 + step.size(k) * g_beta0
    beta = beta + step.size(k) * g_beta
    gamma_vec = gamma_vec + step.size(k) * g_gamma
    mu = mu + step.size(k) * g_mu
    
    
    marglik_liste[k] = marg_lik(Y, beta0, beta, X, sigma2, mu ,gamma_vec,n_mis)
    
    
    beta0_liste[k+1] = beta0 
    beta_liste[k+1,] = beta
    gamma_liste[k+1,] = gamma_vec
    mu_liste[k+1,] = mu
    
    ###
    # on refait pareil pour H0
    ###
    
    X2_nouv = rep(0,n_mis)
    
    # faire un vectoeur moyenne et une matrice variance et faire rmnvnorm
    
    for (i in 1:n_mis){
      
      X2_nouv[i] = rnorm(1, mean =moy_post(Y[i],X[i,1],mu = muH0[2],beta0 = beta0H0,beta1 = betaH0[1],beta2 = betaH0[2],gamma2 = gamma_vecH0[2],sigma2 = 1), sd = sqrt(var_post(beta = betaH0,gamma2 = gamma_vecH0[2],sigma2=1)) )
      
    }
    sample_listeH0[k,] = X2_nouv
    
    Xk = X
    Xk[1:n_mis,2] = X2_nouv
    
    g_beta0 = 0
    g_beta = 0
    g_mu = 0
    g_gamma=0
    
    for (i in 1:n){
      #g_beta0 = g_beta0 + grad_beta0(Y[i], beta0H0, betaH0, Xk[i,], sigma2_star)/n
      g_beta = g_beta +grad_beta(Y[i], beta0H0, betaH0, Xk[i,], sigma2_star)/n
      
      g_mu = g_mu + grad_mu(Xk[i,],muH0,gamma_vecH0)/n
      g_gamma = g_gamma + grad_gamma(Xk[i,],muH0,gamma_vecH0)/n
      
    }
    
    #beta0H0 = beta0H0 + step.size(k) * g_beta0
    betaH0 = betaH0 + step.size(k) * g_beta
    gamma_vecH0 = gamma_vecH0 + step.size(k) * g_gamma
    muH0 = muH0 + step.size(k) * g_mu
    
    marglik_listeH0[k] = marg_lik(Y, beta0=0, betaH0, X, sigma2, muH0 ,gamma_vecH0,n_mis)
    
    beta0_listeH0[k+1] = beta0H0 
    beta_listeH0[k+1,] = betaH0
    gamma_listeH0[k+1,] = gamma_vecH0
    mu_listeH0[k+1,] = muH0
    
    
    
    print(paste("t=",t,"k=",k))
  }
  
  
  
  
  # 
  
  
  
  sample_mixte = matrix(0,nrow = K-1,ncol = n_mis)
  
  for (k in 1:(K-1)){
    u = runif(1)
    if(u<0.5){sample_mixte[k,] = sample_listeH0[k,]
    }else{sample_mixte[k,] = sample_liste[k,]
    
    }
  }
  
  
  logf0_sample = matrix(0,nrow =K-1,ncol = n)
  logf1_sample = matrix(0,nrow =K-1,ncol = n)
  
  
  
  for (k in 1:(K-1)){
    
    Xk = X
    Xk[1:n_mis,2] = sample_mixte[k,]
    
    for (i in 1:n){
      logf0_sample[k,i] = complete_llik(Y[i], beta0_listeH0[k], beta_listeH0[k,], Xk[i,], sigma2 = 1, mu_listeH0[k,] ,gamma_listeH0[k,] )
      logf1_sample[k,i] = complete_llik(Y[i], beta0_liste[k], beta_liste[k,], Xk[i,], sigma2 = 1, mu_liste[k,] ,gamma_liste[k,] )
    }
    
  }
  
  
  # une colonne c'est pour un ri 
  r_liste2 = matrix(0,nrow =K,ncol =n)
  step_size = function(x){return(0.1)}
  
  for (k in 1:(K-1)){
    for(i in 1:n){
      
      r_liste2[k+1,i]= r_liste2[k,i] + step_size(k)* (exp(logf0_sample[k,i]) - exp(logf1_sample[k,i]+r_liste2[k,i]))/(exp(logf0_sample[k,i]) + exp(logf1_sample[k,i]+r_liste2[k,i]))
      
    }
  }
  
  lr_liste2[,t] = apply(-2*r_liste2,FUN=sum,MARGIN=1)[2:K]
  lr_exact2[,t] = -2*(marglik_listeH0- marglik_liste)
  
  
}

lr_liste_mean =colMeans(t(lr_liste2))
lr_exact_mean = colMeans(t(lr_exact2))
ite = rep(1:(K-1),nexp)



lr_sd = apply(lr_liste2,FUN=sd,MARGIN=1)
lr_exact_sd = apply(lr_exact2,FUN=sd,MARGIN=1)

#




diff_mse = (lr_liste2-lr_exact2)**2
diff_mse_mean = apply(diff_mse,FUN=mean,MARGIN=1)


quantiles_mse = t(apply(diff_mse, 1, quantiles))
quantiles_mse = as.data.frame(quantiles_mse)
colnames(quantiles_mse) = c("q5", "q95")


# Calculer les moyennes pour diff_mse
diff_mse_mean = apply(diff_mse, FUN = mean, MARGIN = 1)


ite_vec = c(50,100,150,200,249)

data_approx = data.frame("lr" = lr_liste_mean, "sd" =lr_sd, "method" = rep("saris", length(lr_liste_mean)),"iterations" = 1:length(lr_liste_mean))[ite_vec,]
data_exact = data.frame("lr" = lr_exact_mean, "sd" =lr_exact_sd, "method" = rep("exact", length(lr_exact_mean)), "iterations" = 1:length(lr_exact_mean))[ite_vec,]

data_joint = rbind(data_approx,data_exact)
data_mse = cbind(data.frame("mse" = diff_mse_mean,"iterations"= 1:length(diff_mse_mean)), quantiles_mse )[ite_vec,]

plot_joint <- ggplot(data_joint, aes(x = iterations, y = lr, color = method)) +
  geom_point(size = 5, position = position_dodge(width = 10)) +  # Ajouter les points avec un d?calage
  geom_errorbar(aes(ymin = lr- sd, ymax = lr + sd),size=1, width = 10, position = position_dodge(width = 10)) +  # Ajouter les barres d'erreur avec le m?me d?calage
  labs(x = "iterations", y = "LRT") +  # Ajouter des ?tiquettes d'axes
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) + ylim(-2, 2)

plot_joint

plot_mse <- ggplot(data_mse, aes(x = iterations, y = mse)) +
  geom_point(size = 5, position = position_dodge(width = 0.9)) +  # Ajouter les points avec un d?calage
  geom_errorbar(aes(ymin = q5, ymax = q95),size=1, width = 10, position = position_dodge(width = 0.9)) +  # Ajouter les barres d'erreur avec le m?me d?calage
  labs(x = "iterations", y = "MSE") +  # Ajouter des ?tiquettes d'axes
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept=0) + ylim(-0.1, 0.5)

plot_mse

x11()
grid.arrange(plot_joint, plot_mse, ncol=2)
dev.off()

jpeg("lr_n25.jpg",width = 900,height=450)
plot_joint
dev.off()


jpeg("mse_n25.jpg",width = 900,height=450)
plot_mse
dev.off()

write.csv2(data.frame('lr' = lr_liste2, "lr_exact" = lr_exact2), file = 'lr_n25.csv')
