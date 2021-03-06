# ---------------------------------------------- #
# Example code for sLDA-M
# 
# Code adapted from Glynn et al. (2019):
# https://github.com/G-Lynn/DLTM
# Reference:
# Glynn, C., Tokdar, S. T., Howard, B., & Banks, D. L. (2019). 
# Bayesian analysis of dynamic linear topic models.Bayesian Analysis, 14(1), 53-80.
# ---------------------------------------------- #

library(MASS)  
library(gtools)
library(Rcpp)
library(RcppArmadillo)
library(clue)

set.seed(51)
# --------------- Data generation ---------------------- #
source("sub_data_generation.R")

#---------------- Load utility functions --------------- #
sourceCpp("Cpp/z_step_slda_M.cpp")


K=4  # number of topics
burnIn = 1000 
nSamples = 2000 
thin =  1 
nIter = thin*nSamples # total number of iterations


#----Specify hyper-parameters -------------------
## Beta
beta_p = array(1, dim=c(K,V))

## Mean intercept and slope
beta_0 = rep(0,K)
Sigma_0 = 100^2*diag(K)

## Theta 
theta_p = array(1, dim=c(nDocs,K))

# ------------- MCMC initialization ------------------- #

## Topic (word distribution)
gam.hist = array(NA, dim=c(K,V,nSamples))   # storage for gamma
gam0 = matrix(NA, nrow=K, ncol=V)
for(k in 1:K) gam0[k,] = rdirichlet(1,alpha=beta_p[k,])


## Theta
theta.hist = array(NA, dim=c(nDocs,K,nSamples))   # storage for theta
theta0 = matrix(NA, nrow=nDocs, ncol=K)
for(d in 1:nDocs) theta0[d,] = rdirichlet(1,alpha=theta_p[d,])

## Z 
Z0 = sample(c(1:K), sum(N.D), replace=T )

## y_kv
y_kv = array(0, dim=c(K,V))
wz = as.data.frame(cbind(W[,1], Z0))
y.tmp = table(wz$V1, wz$Z0)
for(k in 1:K){
  y_kv[k,] = y.tmp[,k]
}

## x_dk
zw = as.data.frame(cbind(W[,2], Z0))
x.tmp = aggregate(zw, by=list(zw$V1), function(x) table(factor(x, levels=1:K)))
x_dk = x.tmp[,3]

## regression Coefficients
b1.hist = array(NA, dim=c(K,nSamples))
b1 = mvrnorm(1,mu=beta_0, Sigma=0.1*Sigma_0)
b2.hist = array(NA, dim=c(K,nSamples))
b2 = mvrnorm(1,mu=beta_0, Sigma=0.1*Sigma_0)

## error variance
ve = 0.01 #3
se = 0.01 #9

sig1 = ceiling(abs(rnorm(1,0,1000)))
sig1.hist = c()
sig2 = ceiling(abs(rnorm(1,0,1000)))
sig2.hist = c()

#----MCMC-------------------------------------------------------


b=0
for(m in 1:nIter){
  
  #-----gamma Step (topic) -------------------------------------"
  post.beta = y_kv + beta_p
  for(k in 1:K) gam0[k,] = rdirichlet(1,alpha=post.beta[k,])
  
  #-----Theta Step-------------------------------------
  post.theta = x_dk + theta_p
  for(d in 1:nDocs) theta0[d,] = rdirichlet(1,alpha=post.theta[d,])
  
  #---Z Step -----------------------------------------
  zstep = Zstep(gam0, W, Z0, theta0,  y_kv, x_dk, b1, b2, sig1, sig2, mem, N.D, resp )
  Z0 = zstep[[1]]
  y_kv = zstep[[2]]
  x_dk = zstep[[3]]
  
  #--- Beta Step 1 -----------------------------------------
  z.bar1 = (x_dk[mem==1,]/rowSums(x_dk[mem==1,]))
  Vb1 = solve(solve(Sigma_0) + 1/sig1 * t(z.bar1)%*%z.bar1)
  mu_b1 = Vb1%*%(solve(Sigma_0)%*%beta_0 + (1/sig1)*t(z.bar1)%*%resp[mem==1])
  b1 = mvrnorm(1, mu=mu_b1, Sigma=Vb1)
  
  #--- Beta Step 2 -----------------------------------------
  z.bar2 = (x_dk[mem==2,]/rowSums(x_dk[mem==2,]))
  Vb2 = solve(solve(Sigma_0) + 1/sig2 * t(z.bar2)%*%z.bar2)
  mu_b2 = Vb2%*%(solve(Sigma_0)%*%beta_0 + (1/sig2)*t(z.bar2)%*%resp[mem==2])
  b2 = mvrnorm(1, mu=mu_b2, Sigma=Vb2)
  
  #-----sig1 Step----------------------------------------
  ae1 = 0.5*sum(mem==1) + ve
  be1 = se + 0.5*t(resp[mem==1] - z.bar1%*%b1)%*%(resp[mem==1] - z.bar1%*%b1)
  sig1 = 1/rgamma(1,shape=ae1, rate=be1)

  #-----sig2 Step----------------------------------------
  ae2 = 0.5*sum(mem==2) + ve
  be2 = se + 0.5*t(resp[mem==2] - z.bar2%*%b2)%*%(resp[mem==2] - z.bar2%*%b2)
  sig2 = 1/rgamma(1,shape=ae2, rate=be2)

  if(m %% thin ==0){
    b = b+1
    gam.hist[,,b] = gam0
    theta.hist[,,b] = theta0
    b1.hist[,b] = b1
    b2.hist[,b] = b2
    sig1.hist[b] = sig1
    sig2.hist[b] = sig2
  }
  if(m%%200==0) message(paste("Sample: ",m,sep=""))
}




#---- Result summary --------------------------------------------------

# Topics
gam00 = gam.hist[,,(burnIn+1):nSamples]
gam_s =apply(gam00,c(1,2),mean)
gams = t(gam_s)

par(mfrow=c(1,4))
plot(gams[,1], type = 'l')
plot(gams[,2], type = 'l')
plot(gams[,3], type = 'l')
plot(gams[,4], type = 'l')


gam_dot = t.gam %*% gams
new.order = as.numeric(solve_LSAP(gam_dot, maximum = T))
gam.new = gams[,new.order]

par(mfrow=c(1,4))
plot(gam.new[,1], type = 'l')
plot(gam.new[,2], type = 'l')
plot(gam.new[,3], type = 'l')
plot(gam.new[,4], type = 'l')



# ----------- Trace plots -------------- #
## Regression coefficients
par(mfrow=c(2,2),   mar = c(3,3,1,1) + 1)
plot(b1.hist[1,],type='l', xlab='Iteration', ylab='b11')
plot(b1.hist[2,],type='l', xlab='Iteration', ylab='b12')
plot(b1.hist[3,],type='l', xlab='Iteration', ylab='b13')
plot(b1.hist[4,],type='l', xlab='Iteration', ylab='b14')

par(mfrow=c(2,2),   mar = c(3,3,1,1) + 1)
plot(b2.hist[1,],type='l', xlab='Iteration', ylab='b21')
plot(b2.hist[2,],type='l', xlab='Iteration', ylab='b22')
plot(b2.hist[3,],type='l', xlab='Iteration', ylab='b23')
plot(b2.hist[4,],type='l', xlab='Iteration', ylab='b24')

par(mfrow=c(2,2),   mar = c(3,3,1,1) + 1)
plot(density(b1.hist[1,(burnIn+1):nSamples]), ylab='b11', main='')
plot(density(b1.hist[2,(burnIn+1):nSamples]), ylab='b12', main='')
plot(density(b1.hist[3,(burnIn+1):nSamples]), ylab='b13', main='')
plot(density(b1.hist[4,(burnIn+1):nSamples]), ylab='b14', main='')

par(mfrow=c(2,2),   mar = c(3,3,1,1) + 1)
plot(density(b2.hist[1,(burnIn+1):nSamples]), ylab='b21', main='')
plot(density(b2.hist[2,(burnIn+1):nSamples]), ylab='b22', main='')
plot(density(b2.hist[3,(burnIn+1):nSamples]), ylab='b23', main='')
plot(density(b2.hist[4,(burnIn+1):nSamples]), ylab='b24', main='')

## Error variance
par(mfrow=c(1,2),   mar = c(3,3,1,1) + 1)
plot(sig1.hist,type='l', xlab='Iteration', ylab='sig1^2')
plot(sig2.hist,type='l', xlab='Iteration', ylab='sig2^2')

par(mfrow=c(1,2))
plot(density(sig1.hist[(burnIn+1):nSamples]), ylab='sig1^2', main='')
plot(density(sig2.hist[(burnIn+1):nSamples]), ylab='sig2^2', main='')

# ------- Parameter estimates ------------ #
## Regression coefficients
round(apply(b1.hist[new.order,(burnIn+1):nSamples],1,mean),2)
round(apply(b2.hist[new.order,(burnIn+1):nSamples],1,mean),2)

round(apply(b1.hist[new.order,(burnIn+1):nSamples],1,function(x) quantile(x, probs=c(0.025, 0.975))),2)
round(apply(b2.hist[new.order,(burnIn+1):nSamples],1,function(x) quantile(x, probs=c(0.025, 0.975))),2)

## Error variance
round(mean(sig1.hist[(burnIn+1):nSamples]),2)
round(mean(sig2.hist[(burnIn+1):nSamples]),2)

round(quantile(sig1.hist[(burnIn+1):nSamples], probs=c(0.025, 0.975)),2)
round(quantile(sig2.hist[(burnIn+1):nSamples], probs=c(0.025, 0.975)),2)
