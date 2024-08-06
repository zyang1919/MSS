#!/usr/bin/env Rscript

## clean history 
rm(list = ls())
library(MASS)
library(LaplacesDemon)


##########################
####### read data ########
##########################
sur_sim <- "/nfs/ihfs/home_metis/fluan"
df <- get(load(paste(sur_sim, "Simulated 10,000.RData", sep = '/')))

ftrue <- get(load(paste(sur_sim, "True F.RData", sep = '/')))
f1true <- ftrue[[1]]
f2true <- ftrue[[2]]
f3true <- ftrue[[3]]
getK <- function(x){
  # This function returns K matrix, and requires unique design points x
  
  n <- length(x)
  dx <- diff(x)
  dx2 <- diff(x,lag=2)
  Qt <- matrix(rep(0,(n-2)*n),n-2,n)
  for (k in 1:(n-2)){
    Qt[k,k]   <- 1/dx[k]
    Qt[k,k+1] <- -1/dx[k+1]-1/dx[k]
    Qt[k,k+2] <- 1/dx[k+1]
  }
  R  <- matrix(rep(0,(n-2)*(n-2)),n-2,n-2)
  R[1,1] <- dx[1]
  R[1,2] <- dx[2]
  for (k in 2:(n-3)){
    R[k,k-1] <- dx[k]
    R[k,k]   <- 2*dx2[k]
    R[k,k+1] <- dx[k+1]
  }
  R[n-2,n-3] <- dx[n-2]
  R[n-2,n-2] <- dx[n-1]
  R <- R/6
  invR <- solve(R)
  K <- t(Qt)%*%invR%*%Qt
  K
  
  
}
##########################
##########################

##########################
###### run simulation ####
##########################
df_final <- df[8001:10000]
nsim <- length(df_final)

sim_tau1 <- rep(0, nsim); sim_tau2 <- rep(0, nsim); sim_tau3 <- rep(0, nsim)
sim_sig1 <- rep(0, nsim); sim_sig2 <- rep(0, nsim); sim_sig3 <- rep(0, nsim)
sim_alp1 <- rep(0, nsim); sim_alp2 <- rep(0, nsim); sim_alp3 <- rep(0, nsim)
sim_rho12 <- rep(0, nsim); sim_rho13 <- rep(0, nsim); sim_rho23 <- rep(0, nsim)
sim_f1 <- list(); sim_f2 <- list(); sim_f3 <- list()
imse_all <- rep(0, nsim)

for(s in 1:nsim){
  
  ## get data
  Y <- df_final[[s]]
  y1 <- as.matrix(Y[,1], nrow = 100)
  y2 <- as.matrix(Y[,2], nrow = 100)
  y3 <- as.matrix(Y[,3], nrow = 100)
  Tm <- seq(-2, 2, length = 100)
  m <- 100
  Im <- diag(1, m)
  K <- getK(Tm)
  
  sig1t <- var(c(y1))
  sig2t <- var(c(y2))
  sig3t <- var(c(y3))
  
  r12 <- cor(y1,y2); r13 <- cor(y1,y3); r23 <- cor(y2,y3)
  zt12t <- 0.5*log((1+r12)/(1-r12))
  zt13t <- 0.5*log((1+r13)/(1-r13))
  zt23t <- 0.5*log((1+r23)/(1-r23))
  
  zt12 <-0.1; zt13 <- 0.1; zt23 <- 0.1
  
  tau1t <- tau2t <- tau3t <- 3
  
  rho12t <- 0.3
  rho13t <- 0.1
  rho23t <- 0.2
  
  alp1t <- tau1t*sig1t
  alp2t <- tau2t*sig2t
  alp3t <- tau3t*sig3t
  
  f1t <- rep(1, length(Tm)) 
  f2t <- rep(1, length(Tm))
  f3t <- rep(1, length(Tm))
  
  ## hyper paramter 
  At1 <- 1; At2 <- 1; At3 <- 1;
  Bt1 <- 1; Bt2 <- 1; Bt3 <- 1;
  
  As1 <- 1; As2 <- 1; As3 <- 1;
  Bs1 <- 1; Bs2 <- 1; Bs3 <- 1;
  
  
  niter <- 2000 # samples 
  nburn <- 2000 # burning 
  aiter <- niter + nburn # all 
  
  alp1tt <- rep(0, niter)
  alp2tt <- rep(0, niter)
  alp3tt <- rep(0, niter)
  
  f1tt <- matrix(0, niter, length(Tm)) 
  f2tt <- matrix(0, niter, length(Tm))
  f3tt <- matrix(0, niter, length(Tm))
  
  sig1tt <- rep(0, niter)
  sig2tt <- rep(0, niter)
  sig3tt <- rep(0, niter)
  
  rho12tt <- rep(0, niter)
  rho13tt <- rep(0, niter)
  rho23tt <- rep(0, niter)
  
  tau1tt <- rep(0, niter)
  tau2tt <- rep(0, niter)
  tau3tt <- rep(0, niter)
  
  for(k in 1:aiter){
    
    #sig1t <- sig2t <- sig3t <- 1
    
    ## get the inverse matrix for each block 
    deno <- (1-rho13t^2)*(1-rho12t^2) - (rho23t - rho12t*rho13t)^2
    iSigma11t <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13t)^2)*(1-rho12t^2))/(deno)))
    iSigma11t <- iSigma11t*Im
    
    iSigma22t <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13t)^2)/deno)
    iSigma22t <- iSigma22t*Im
    
    iSigma33t <- (1/sig3t)*((1-rho12t^2)/deno)
    iSigma33t <- iSigma33t*Im
    
    iSigma12t <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)*rho13t/deno))
    iSigma12t <- iSigma12t*Im
    
    iSigma13t <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13t <- iSigma13t*Im
    
    iSigma23t <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23t <- iSigma23t*Im
    ## end of inverse of block 
    
    ## identity matrix 
    Im <- diag(1, 100)
    
    ## update f
    ## update f1t
    ## get inverse of Sigma1t
    #iSig <- getiSig(sig1t,sig2t,sig3t,rho12t,rho13t,rho23t,100)
    #iSigma11t <- iSig$invSig11
    #iSigma12t <- iSig$invSig12
    #iSigma13t <- iSig$invSig13
    
    
    iA1 <- iSigma11t + c(alp1t)*K
    s1 <- svd(iA1)
    
    D1 <- solve(diag(s1$d))
    Aalp1 <- s1$v%*%D1%*%t(s1$u)
    Halp1 <- s1$v%*%sqrt(D1)
    
    
    #  A1 <- solve(iA1)
    
    
    B1 <- iSigma11t%*%y1 + iSigma12t%*%(y2 - f2t) + iSigma13t%*%(y3 - f3t)
    
    Z1 <- rnorm(m, 0, 1)
    f1t <- (Aalp1%*%B1) + c(sqrt(sig1t))*(Halp1%*%Z1)
    
    
    #  Mu1 <- A1%*%B1
    #  Var1 <- A1
    
    
    
    #  f1t <- matrix(mvrnorm(1, Mu1, Var1), nrow = length(Tm))
    
    ## update f2t
    ## get inverse of Sigma2t
    #iSigma22t <- iSig$invSig22
    #iSigma12t <- iSig$invSig12
    #iSigma23t <- iSig$invSig23
    
    
    iA2 <- iSigma22t + c(alp2t)*K
    s2 <- svd(iA2)
    
    D2 <- solve(diag(s2$d))
    Aalp2 <- s2$v%*%D2%*%t(s2$u)
    Halp2 <- s2$v%*%sqrt(D2)
    #  A2 <- solve(iA2)
    
    
    
    B2 <- iSigma22t%*%y2 + iSigma12t%*%(y1 - f1t) + iSigma23t%*%(y3 - f3t)
    
    Z2 <- rnorm(m, 0, 1)
    f2t <- (Aalp2%*%B2) + c(sqrt(sig2t))*(Halp2%*%Z2)
    
    #  Mu2 <- A2%*%B2
    #  Var2 <- A2
    #  Var2 <- nearPD(Var2)$mat
    
    #  f2t <- matrix(mvrnorm(1, Mu2, Var2), nrow = length(Tm))
    
    ## update f3t
    ## get invers of Sigma3t
    #iSigma33t <- iSig$invSig33
    #iSigma13t <- iSig$invSig13
    #iSigma23t <- iSig$invSig23
    
    
    iA3 <- iSigma33t + c(alp3t)*K
    s3 <- svd(iA3)
    
    D3 <- solve(diag(s3$d))
    Aalp3 <- s3$v%*%D3%*%t(s3$u)
    Halp3 <- s3$v%*%sqrt(D3)
    #  A3 <- solve(iA3)
    
    
    
    B3 <- iSigma33t%*%y3 + iSigma13t%*%(y1 - f1t) + iSigma23t%*%(y2 - f2t)
    
    Z3 <- rnorm(m, 0, 1)
    f3t <- (Aalp3%*%B3) + c(sqrt(sig3t))*(Halp3%*%Z3)
    
    #  Mu3 <- A3%*%B3
    #  Var3 <- A3
    #  Var3 <- nearPD(Var3)$mat
    
    #  f3t <- matrix(mvrnorm(1, Mu3, Var3), nrow = length(Tm))
    ### end of updating f ######
    
    
    ## update tau_i
    ## update tau1
    At1star <- At1 + (m-2)/2
    Bt1star <- Bt1 + 0.5*(t(f1t)%*%K%*%f1t)
    
    tau1t <- rgamma(1, At1star, Bt1star)
    
    ## update tau2
    At2star <- At2 + (m-2)/2
    Bt2star <- Bt2 + 0.5*(t(f2t)%*%K%*%f2t)
    
    tau2t <- rgamma(1, At2star, Bt2star)
    
    ## update tau3
    At3star <- At3 + (m-2)/2
    Bt3star <- Bt3 + 0.5*(t(f3t)%*%K%*%f3t)
    
    tau3t <- rgamma(1, At3star, Bt3star)
    ## end of updating tau_i
    
    
    ## update sigma_i
    ## update sigma 1
    ## proposal IG
    sig1c <- rinvgamma(1, shape = (var(y1) + 2), scale = (var(y1) + 1)*var(y1))
    #  sig1c <- rgamma(1, sig1t^2, rate = sig1t)
    ## get inverse for all iSigma1jt
    iSigma11c <- (1/sig1c)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13t)^2)*(1-rho12t^2))/(deno)))
    iSigma11c <- iSigma11c*Im
    
    iSigma12c <- (-1/(sqrt(sig1c)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)*rho13t/deno))
    iSigma12c <- iSigma12c*Im
    
    iSigma13c <- (-1/(sqrt(sig1c)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13c <- iSigma13c*Im
    
    RSSs1_sig1c <- t(y1 - f1t)%*%iSigma11c%*%(y1 - f1t)
    RSSs2_sig1c <- t(y1 - f1t)%*%iSigma12c%*%(y2 - f2t)
    RSSs3_sig1c <- t(y1 - f1t)%*%iSigma13c%*%(y3 - f3t)
    
    numers1 <- -(m/2+As1+1)*log(sig1c)+(-0.5*(RSSs1_sig1c+2*RSSs2_sig1c+2*RSSs3_sig1c)+Bs1/sig1c)
    
    
    ## target 
    RSSs1_sig1t <- t(y1 - f1t)%*%iSigma11t%*%(y1 - f1t)
    RSSs2_sig1t <- t(y1 - f1t)%*%iSigma12t%*%(y2 - f2t)
    RSSs3_sig1t <- t(y1 - f1t)%*%iSigma13t%*%(y3 - f3t)
    
    denoms1 <- -(m/2+As1+1)*log(sig1t)+(-0.5*(RSSs1_sig1t+2*RSSs2_sig1t+2*RSSs3_sig1t)+Bs1/sig1t)
    
    
    ## pdf of proposal 
    logsig1c <- dinvgamma(sig1t, shape = (sig1c^2 + 2), scale = ((sig1c^2 + 1)*sig1c), log = TRUE)
    #  logsig1c <- dgamma(sig1t, shape = sig1c^2, rate = sig1c, log = TRUE)
    logsig1t <- dinvgamma(sig1c, shape = (sig1t^2 + 2), scale = ((sig1t^2 + 1)*sig1t), log = TRUE)
    #  logsig1t <- dgamma(sig1c, shape = sig1t^2, rate = sig1t, log = TRUE)
    
    ## transition 
    logalpha_sig1 <- numers1 - denoms1 + logsig1c - logsig1t
    
    ## move or not 
    logu_sig1 <- log(runif(1, 0, 1))
    if(logu_sig1<logalpha_sig1){
      sig1t <- sig1c
    }
    
    ## update sigma 2
    ## proposal IG
    sig2c <- rinvgamma(1, shape = (var(y2) + 2), scale = ((var(y2) + 1)*var(y2)))
    #  sig2c <- rgamma(1, sig2t^2, rate = sig2t)
    ## get inverse for all iSigma2jc
    iSigma22c <- (1/(sig2c*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13t)^2)/deno)
    iSigma22c <- iSigma22c*Im
    
    iSigma12c <- (-1/(sqrt(sig1t)*sqrt(sig2c)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)/deno))
    iSigma12c <- iSigma12c*Im
    iSigma12t <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)/deno))
    iSigma12t <- iSigma12t*Im
    
    iSigma23c <- (-1/(sqrt(sig2c)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23c <- iSigma23c*Im
    
    
    RSSs1_sig2c <- t(y2 - f2t)%*%iSigma12c%*%(y1 - f1t)
    RSSs2_sig2c <- t(y2 - f2t)%*%iSigma22c%*%(y2 - f2t)
    RSSs3_sig2c <- t(y2 - f2t)%*%iSigma23c%*%(y3 - f3t)
    
    numers2 <- -(m/2+As2+1)*log(sig2c)+(-0.5*(2*RSSs1_sig2c+RSSs2_sig2c+2*RSSs3_sig2c)+Bs2/sig2c)
    
    
    ## target 
    RSSs1_sig2t <- t(y2 - f2t)%*%iSigma12t%*%(y1 - f1t)
    RSSs2_sig2t <- t(y2 - f2t)%*%iSigma22t%*%(y2 - f2t)
    RSSs3_sig2t <- t(y2 - f2t)%*%iSigma23t%*%(y3 - f3t)
    
    denoms2 <- -(m/2+As2+1)*log(sig2t)+(-0.5*(2*RSSs1_sig2t+RSSs2_sig2t+2*RSSs3_sig2t)+Bs2/sig2t)
    
    
    ## pdf of proposal 
    logsig2c <- dinvgamma(sig2t, shape = (sig2c^2 + 2), scale = ((sig2c^2 + 1)*sig2c), log = TRUE)
    #  logsig2c <- dgamma(sig2t, shape = sig2c^2, rate = sig2c, log = TRUE)
    logsig2t <- dinvgamma(sig2c, shape = (sig2t^2 + 2), scale = ((sig2t^2 + 1)*sig2t), log = TRUE)
    #  logsig2t <- dgamma(sig2c, shape = sig2t^2, rate = sig2t, log = TRUE)
    
    ## transition 
    logalpha_sig2 <- numers2 - denoms2 + logsig2c - logsig2t
    
    ## move or not 
    logu_sig2 <- log(runif(1, 0, 1))
    if(logu_sig2<logalpha_sig2){
      sig2t <- sig2c
    }
    
    
    ## update sigma 3
    ## proposal IG
    sig3c <- rinvgamma(1, shape = (var(y3) + 2), scale = ((var(y3) + 1)*var(y3)))
    #  sig3c <- rgamma(1, sig3t^2, rate = sig3t)
    ## get inverse for all iSigma3jc
    iSigma33c <- (1/sig3c)*((1-rho12t^2)/deno)
    iSigma33c <- iSigma33c*Im
    
    iSigma13c <- (-1/(sqrt(sig1t)*sqrt(sig3c)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13c <- iSigma13c*Im
    iSigma13t <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13t <- iSigma13t*Im
    
    iSigma23c <- (-1/(sqrt(sig2t)*sqrt(sig3c)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23c <- iSigma23c*Im
    iSigma23t <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23t <- iSigma23t*Im
    
    
    RSSs1_sig3c <- t(y3 - f3t)%*%iSigma13c%*%(y1 - f1t)
    RSSs2_sig3c <- t(y3 - f3t)%*%iSigma23c%*%(y2 - f2t)
    RSSs3_sig3c <- t(y3 - f3t)%*%iSigma33c%*%(y3 - f3t)
    
    numers3 <- -(m/2+As3+1)*log(sig3c)+(-0.5*(2*RSSs1_sig3c+2*RSSs2_sig3c+RSSs3_sig3c)+Bs3/sig3c)
    
    
    ## target 
    
    RSSs1_sig3t <- t(y3 - f3t)%*%iSigma13t%*%(y1 - f1t)
    RSSs2_sig3t <- t(y3 - f3t)%*%iSigma23t%*%(y2 - f2t)
    RSSs3_sig3t <- t(y3 - f3t)%*%iSigma33t%*%(y3 - f3t)
    
    denoms3 <- -(m/2+As3+1)*log(sig3t)+(-0.5*(2*RSSs1_sig3t+2*RSSs2_sig3t+RSSs3_sig3t)+Bs3/sig3t)
    
    
    ## pdf of proposal 
    logsig3c <- dinvgamma(sig3t, shape = (sig3c^2 + 2), scale = ((sig3c^2 + 1)*sig3c), log = TRUE)
    #  logsig3c <- dgamma(sig3t, shape = sig3c^2, rate = sig3c, log = TRUE)
    logsig3t <- dinvgamma(sig3c, shape = (sig3t^2 + 2), scale = ((sig3t^2 + 1)*sig3t), log = TRUE)
    #  logsig3t <- dgamma(sig3c, shape = sig3t^2, rate = sig3t, log = TRUE)
    
    ## transition 
    logalpha_sig3 <- numers3 - denoms3 + logsig3c - logsig3t
    
    ## move or not 
    logu_sig3 <- log(runif(1, 0, 1))
    if(logu_sig3<logalpha_sig3){
      sig3t <- sig3c
    }
    ## end of updating sigma_i
    
    
    ## update rho_ij
    ## update rho12
    #  rho12c <- runif(1, -1, 1)
    zc12 <- rnorm(1, zt12t,sd =.001)
    rho12c <- (exp(2*zc12)-1)/(exp(2*zc12)+1)
    rho12t <- (exp(2*zt12)-1)/(exp(2*zt12)+1)
    ## get inverse of iSigma12
    denor12c <- (1-rho13t^2)*(1-rho12c^2) - (rho23t - rho12c*rho13t)^2
    iSigma12c <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12c^2)))*(rho12c-((rho13t-rho12c*rho23t)*(rho23t-rho13t*rho12c)/denor12c))
    iSigma12c <- iSigma12t*Im
    iSigma12t <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)/deno))
    iSigma12t <- iSigma12t*Im
    
    iSigma11c <- (1/sig1t)*(1 + (((rho12c^2)/(1-rho12c^2))+(((rho23t-rho12c*rho13t)^2)*(1-rho12c^2))/(denor12c)))
    iSigma11c <- iSigma11c*Im
    iSigma11t <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13t)^2)*(1-rho12t^2))/(deno)))
    iSigma11t <- iSigma11t*Im
    
    iSigma22c <- (1/(sig2t*(1-rho12c^2)))*(1+((rho23t-rho12c*rho13t)^2)/denor12c)
    iSigma22c <- iSigma22c*Im
    iSigma22t <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13t)^2)/deno)
    iSigma22t <- iSigma22t*Im
    
    iSigma33c <- (1/sig3t)*((1-rho12c^2)/denor12c)
    iSigma33c <- iSigma33c*Im
    iSigma33t <- (1/sig3t)*((1-rho12t^2)/deno)
    iSigma33t <- iSigma33t*Im
    
    iSigma13c <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12c*rho23t)/denor12c)
    iSigma13c <- iSigma13c*Im
    iSigma13t <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13t <- iSigma13t*Im
    
    iSigma23c <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12c*rho13t)/denor12c)
    iSigma23c <- iSigma23c*Im
    iSigma23t <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23t <- iSigma23t*Im
    
    ## RSS ij
    RSS11_r12c <- t(y1-f1t)%*%iSigma11c%*%(y1-f1t)
    RSS11_r12t <- t(y1-f1t)%*%iSigma11t%*%(y1-f1t)
    
    RSS22_r12c <- t(y2-f2t)%*%iSigma22c%*%(y2-f2t)
    RSS22_r12t <- t(y2-f2t)%*%iSigma22t%*%(y2-f2t)
    
    RSS33_r12c <- t(y3-f3t)%*%iSigma33c%*%(y3-f3t)
    RSS33_r12t <- t(y3-f3t)%*%iSigma33t%*%(y3-f3t)
    
    RSS12_r12c <- t(y1-f1t)%*%iSigma12c%*%(y2-f2t)
    RSS12_r12t <- t(y1-f1t)%*%iSigma12t%*%(y2-f2t)
    
    RSS13_r12c <- t(y1-f1t)%*%iSigma13c%*%(y3-f3t)
    RSS13_r12t <- t(y1-f1t)%*%iSigma13t%*%(y3-f3t)
    
    RSS23_r12c <- t(y2-f2t)%*%iSigma23c%*%(y3-f3t)
    RSS23_r12t <- t(y2-f2t)%*%iSigma23t%*%(y3-f3t)
    
    RSS_r12c <- RSS11_r12c+RSS22_r12c+RSS33_r12c+2*RSS12_r12c+2*RSS13_r12c+2*RSS23_r12c
    RSS_r12t <- RSS11_r12t+RSS22_r12t+RSS33_r12t+2*RSS12_r12t+2*RSS13_r12t+2*RSS23_r12t
    
    # transition
    logzc12 <- 2*zc12 - 2*log(exp(2*zc12)+1)
    logzt12 <- 2*zt12 - 2*log(exp(2*zt12)+1)
    
    
    
    numerr12 <- (-m/2)*log((1-rho13t^2)*(1-rho23t^2)-(rho12c-rho13t*rho23t)^2)+(-0.5)*(RSS_r12c)
    denomr12 <- (-m/2)*log((1-rho13t^2)*(1-rho23t^2)-(rho12t-rho13t*rho23t)^2)+(-0.5)*(RSS_r12t)
    
    
    logalpharh12 <- numerr12-denomr12+logzc12 - logzt12
    logu12 <- log(runif(1, 0, 1))
    if(logu12<logalpharh12){
      zt12 <- zc12
    }
    rho12t <- (exp(2*zt12)-1)/(exp(2*zt12)+1)
    
    
    ## update rho13
    zc13 <- rnorm(1, zt13t,sd =.001)
    rho13c <- (exp(2*zc13)-1)/(exp(2*zc13)+1)
    rho13t <- (exp(2*zt13)-1)/(exp(2*zt13)+1)
    ## get inverse of iSigma13c
    deno <- (1-rho13t^2)*(1-rho12t^2) - (rho23t - rho12t*rho13t)^2
    denor13c <- (1-rho13c^2)*(1-rho12t^2) - (rho23t - rho12t*rho13c)^2
    
    iSigma12c <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13c-rho12t*rho23t)*(rho23t-rho13c*rho12t)/denor13c))
    iSigma12c <- iSigma12t*Im
    iSigma12t <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)/deno))
    iSigma12t <- iSigma12t*Im
    
    iSigma11c <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13c)^2)*(1-rho12t^2))/(denor13c)))
    iSigma11c <- iSigma11c*Im
    iSigma11t <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13t)^2)*(1-rho12t^2))/(deno)))
    iSigma11t <- iSigma11t*Im
    
    iSigma22c <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13c)^2)/denor13c)
    iSigma22c <- iSigma22c*Im
    iSigma22t <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13t)^2)/deno)
    iSigma22t <- iSigma22t*Im
    
    iSigma33c <- (1/sig3t)*((1-rho12t^2)/denor13c)
    iSigma33c <- iSigma33c*Im
    iSigma33t <- (1/sig3t)*((1-rho12t^2)/deno)
    iSigma33t <- iSigma33t*Im
    
    iSigma13c <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13c-rho12t*rho23t)/denor13c)
    iSigma13c <- iSigma13c*Im
    iSigma13t <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13t <- iSigma13t*Im
    
    iSigma23c <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13c)/denor13c)
    iSigma23c <- iSigma23c*Im
    iSigma23t <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23t <- iSigma23t*Im
    
    ## RSS ij
    RSS11_r13c <- t(y1-f1t)%*%iSigma11c%*%(y1-f1t)
    RSS11_r13t <- t(y1-f1t)%*%iSigma11t%*%(y1-f1t)
    
    RSS22_r13c <- t(y2-f2t)%*%iSigma22c%*%(y2-f2t)
    RSS22_r13t <- t(y2-f2t)%*%iSigma22t%*%(y2-f2t)
    
    RSS33_r13c <- t(y3-f3t)%*%iSigma33c%*%(y3-f3t)
    RSS33_r13t <- t(y3-f3t)%*%iSigma33t%*%(y3-f3t)
    
    RSS12_r13c <- t(y1-f1t)%*%iSigma12c%*%(y2-f2t)
    RSS12_r13t <- t(y1-f1t)%*%iSigma12t%*%(y2-f2t)
    
    RSS13_r13c <- t(y1-f1t)%*%iSigma13c%*%(y3-f3t)
    RSS13_r13t <- t(y1-f1t)%*%iSigma13t%*%(y3-f3t)
    
    RSS23_r13c <- t(y2-f2t)%*%iSigma23c%*%(y3-f3t)
    RSS23_r13t <- t(y2-f2t)%*%iSigma23t%*%(y3-f3t)
    
    RSS_r13c <- RSS11_r13c+RSS22_r13c+RSS33_r13c+2*RSS12_r13c+2*RSS13_r13c+2*RSS23_r13c
    RSS_r13t <- RSS11_r13t+RSS22_r13t+RSS33_r13t+2*RSS12_r13t+2*RSS13_r13t+2*RSS23_r13t
    
    # transition
    logzc13 <- 2*zc13 - 2*log(exp(2*zc13)+1)
    logzt13 <- 2*zt13 - 2*log(exp(2*zt13)+1)
    
    numerr13 <- (-m/2)*log(1-rho12t^2-rho13c^2-rho23t^2+2*rho12t*rho13c*rho23t)+(-0.5)*(RSS_r13c)
    denomr13 <- (-m/2)*log(1-rho12t^2-rho13t^2-rho23t^2+2*rho12t*rho13t*rho23t)+(-0.5)*(RSS_r13t)
    
    logalpharh13 <- numerr13-denomr13+logzc13 - logzt13
    logu13 <- log(runif(1, 0, 1))
    if(logu13<logalpharh13){
      zt13 <- zc13
    }
    rho13t <- (exp(2*zt13)-1)/(exp(2*zt13)+1)
    
    
    ## update rho23
    zc23 <- rnorm(1, zt23t, sd=.001)
    rho23c <- (exp(2*zc23)-1)/(exp(2*zc23)+1)
    #  while(rho23c>0.5){
    #    zc23 <- rnorm(1, zt23t/2, sd=0.1)
    #    rho23c <- (exp(2*zc23)-1)/(exp(2*zc23)+1)
    #  }
    rho23t <- (exp(2*zt23)-1)/(exp(2*zt23)+1)
    ## get inverse of iSigma12c
    deno <- (1-rho13t^2)*(1-rho12t^2) - (rho23t - rho12t*rho13t)^2
    denor23c <- (1-rho13t^2)*(1-rho12t^2) - (rho23c - rho12t*rho13t)^2
    
    iSigma12c <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23c)*(rho23c-rho13c*rho12t)/denor23c))
    iSigma12c <- iSigma12t*Im
    iSigma12t <- (-1/(sqrt(sig1t)*sqrt(sig2t)*(1-rho12t^2)))*(rho12t-((rho13t-rho12t*rho23t)*(rho23t-rho13t*rho12t)/deno))
    iSigma12t <- iSigma12t*Im
    
    iSigma11c <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23c-rho12t*rho13t)^2)*(1-rho12t^2))/(denor23c)))
    iSigma11c <- iSigma11c*Im
    iSigma11t <- (1/sig1t)*(1 + (((rho12t^2)/(1-rho12t^2))+(((rho23t-rho12t*rho13t)^2)*(1-rho12t^2))/(deno)))
    iSigma11t <- iSigma11t*Im
    
    iSigma22c <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23c-rho12t*rho13t)^2)/denor23c)
    iSigma22c <- iSigma22c*Im
    iSigma22t <- (1/(sig2t*(1-rho12t^2)))*(1+((rho23t-rho12t*rho13t)^2)/deno)
    iSigma22t <- iSigma22t*Im
    
    iSigma33c <- (1/sig3t)*((1-rho12t^2)/denor23c)
    iSigma33c <- iSigma33c*Im
    iSigma33t <- (1/sig3t)*((1-rho12t^2)/deno)
    iSigma33t <- iSigma33t*Im
    
    iSigma13c <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23c)/denor23c)
    iSigma13c <- iSigma13c*Im
    iSigma13t <- (-1/(sqrt(sig1t)*sqrt(sig3t)))*((rho13t-rho12t*rho23t)/deno)
    iSigma13t <- iSigma13t*Im
    
    iSigma23c <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23c-rho12t*rho13t)/denor23c)
    iSigma23c <- iSigma23c*Im
    iSigma23t <- (-1/(sqrt(sig2t)*sqrt(sig3t)))*((rho23t-rho12t*rho13t)/deno)
    iSigma23t <- iSigma23t*Im
    
    ## RSS ij
    RSS11_r23c <- t(y1-f1t)%*%iSigma11c%*%(y1-f1t)
    RSS11_r23t <- t(y1-f1t)%*%iSigma11t%*%(y1-f1t)
    
    RSS22_r23c <- t(y2-f2t)%*%iSigma22c%*%(y2-f2t)
    RSS22_r23t <- t(y2-f2t)%*%iSigma22t%*%(y2-f2t)
    
    RSS33_r23c <- t(y3-f3t)%*%iSigma33c%*%(y3-f3t)
    RSS33_r23t <- t(y3-f3t)%*%iSigma33t%*%(y3-f3t)
    
    RSS12_r23c <- t(y1-f1t)%*%iSigma12c%*%(y2-f2t)
    RSS12_r23t <- t(y1-f1t)%*%iSigma12t%*%(y2-f2t)
    
    RSS13_r23c <- t(y1-f1t)%*%iSigma13c%*%(y3-f3t)
    RSS13_r23t <- t(y1-f1t)%*%iSigma13t%*%(y3-f3t)
    
    RSS23_r23c <- t(y2-f2t)%*%iSigma23c%*%(y3-f3t)
    RSS23_r23t <- t(y2-f2t)%*%iSigma23t%*%(y3-f3t)
    
    RSS_r23c <- RSS11_r23c+RSS22_r23c+RSS33_r23c+2*RSS12_r23c+2*RSS13_r23c+2*RSS23_r23c
    RSS_r23t <- RSS11_r23t+RSS22_r23t+RSS33_r23t+2*RSS12_r23t+2*RSS13_r23t+2*RSS23_r23t
    
    numerr23 <- (-m/2)*log(1-rho12t^2-rho13t^2-rho23c^2+2*rho12t*rho13t*rho23c)+(-0.5)*(RSS_r23c)
    
    denomr23 <- (-m/2)*log(1-rho12t^2-rho13t^2-rho23t^2+2*rho12t*rho13t*rho23t)+(-0.5)*(RSS_r23t)
    
    # transition
    logzc23 <- 2*zc23 - 2*log(exp(2*zc23)+1)
    logzt23 <- 2*zt23 - 2*log(exp(2*zt23)+1)
    
    
    logalpharh23 <- numerr23-denomr23+logzc23 - logzt23
    logu23 <- log(runif(1, 0, 1))
    if(logu23<logalpharh23){
      zt23 <- zc23
    }
    rho23t <- (exp(2*zt23)-1)/(exp(2*zt23)+1)
    ### end of rho_ij
    
    ## alphat
    ## alph1
    alp1t <- tau1t*sig1t
    alp2t <- tau2t*sig2t
    alp3t <- tau3t*sig3t
    ### end of alphat
    
    ## collect samples
    if(k>nburn){
      alp1tt[k-nburn] <- alp1t
      alp2tt[k-nburn] <- alp2t
      alp3tt[k-nburn] <- alp3t
      
      f1tt[k-nburn,] <- f1t
      f2tt[k-nburn,] <- f2t
      f3tt[k-nburn,] <- f3t
      
      sig1tt[k-nburn] <- sig1t
      sig2tt[k-nburn] <- sig2t
      sig3tt[k-nburn] <- sig3t
      
      rho12tt[k-nburn] <- rho12t
      rho13tt[k-nburn] <- rho13t
      rho23tt[k-nburn] <- rho23t
      
      tau1tt[k-nburn] <- tau1t
      tau2tt[k-nburn] <- tau2t
      tau3tt[k-nburn] <- tau3t
      
    }
  }
  
  sim_tau1[s] <- mean(tau1tt); sim_tau2[s] <- mean(tau2tt); sim_tau3[s] <- mean(tau3tt)
  sim_sig1[s] <- mean(sig1tt); sim_sig2[s] <- mean(sig2tt); sim_sig3[s] <- mean(sig3tt)
  sim_alp1[s] <- mean(alp1tt); sim_alp2[s] <- mean(alp2tt); sim_alp3[s] <- mean(alp3tt)
  sim_rho12[s] <- mean(rho12tt); sim_rho13[s] <- mean(rho13tt); sim_rho23[s] <- mean(rho23tt)
  sim_f1[[s]] <- apply(f1tt, 2, mean)
  sim_f2[[s]] <- apply(f2tt, 2, mean) 
  sim_f3[[s]] <- apply(f3tt, 2, mean)

    
  ### calculate MISE 
  f1i <- apply(f1tt, 2, mean); f2i <- apply(f2tt, 2, mean); f3i <- apply(f3tt, 2, mean)
  X_design <- seq(-2, 2, length = 100)

  dx <- diff(X_design)
  nx <- length(X_design)-1

  Ftrue <- cbind(f1true, f2true, f3true)
  Fhatij <- cbind(f1i,f2i,f3i)
  
  FtF <- (Ftrue - Fhatij)%*%t(Ftrue - Fhatij)
  dFtF <- diag(FtF)
#  (diag((Ftrue - Fhatij)%*%t(Ftrue - Fhatij))[1:99])%*%dx[1:99]
  imse_all[s] <- dFtF[1:nx]%*%dx[1:nx]
  
}


Rst.Sim.30 <- list()
Rst.Sim.30[[1]] <- sim_tau1
Rst.Sim.30[[2]] <- sim_tau2
Rst.Sim.30[[3]] <- sim_tau3
Rst.Sim.30[[4]] <- sim_sig1
Rst.Sim.30[[5]] <- sim_sig2
Rst.Sim.30[[6]] <- sim_sig3
Rst.Sim.30[[7]] <- sim_alp1
Rst.Sim.30[[8]] <- sim_alp2
Rst.Sim.30[[9]] <- sim_alp3
Rst.Sim.30[[10]] <- sim_rho12
Rst.Sim.30[[11]] <- sim_rho13
Rst.Sim.30[[12]] <- sim_rho23
Rst.Sim.30[[13]] <- sim_f1
Rst.Sim.30[[14]] <- sim_f2
Rst.Sim.30[[15]] <- sim_f3
Rst.Sim.30[[16]] <- imse_all



names(Rst.Sim.30) <- c("tau1", "tau2", "tau3", 
                         "sig1", "sig2", "sig3", 
                         "alp1", "alp2", "alp3", 
                         "rho12", "rho13", "rho23", 
                         "f1", "f2", "f3", "imse")


save(Rst.Sim.30, file = paste(sur_sim, "Rst.Sim.10000.RData", sep = '/'))



