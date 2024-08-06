#!/usr/bin/env Rscript

## clean history 
rm(list = ls())
library(MASS)
library(LaplacesDemon)


## read data 
sur_sim <- "/nfs/ihfs/home_metis/fluan"
df <- get(load(paste(sur_sim, "Simulated 10,000.RData", sep = '/')))

ftrue <- get(load(paste(sur_sim, "True F.RData", sep = '/')))
f1true <- ftrue[[1]]
f2true <- ftrue[[2]]
f3true <- ftrue[[3]]


## data information
#n <- length(df) ## i = 50 subject in total
m <- 100 ## j = 100 time points for each subjects
q <- 3 ## q = 3 - variate multiresponse

## K function 
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
  sR <- eigen(R)
  sD <- sR$values
  sD <- sD[sD>10^(-4)]
  ni <- length(sD)
  sS <- sR$vectors
  invR <- sS[,1:ni]%*%solve(diag(sD))%*%t(sS[, 1:ni])
  
  K <- t(Qt)%*%invR%*%Qt
  
  K
}

Tm <- seq(-2, 2, length = 100)
K <- getK(Tm)
##############

##########################
###### run simulation ####
##########################
df_final <- df[8001:10000]
nsim <- length(df_final)



sim_alp1 <- rep(0, nsim); sim_alp2 <- rep(0, nsim); sim_alp3 <- rep(0, nsim)
sim_f1 <- list(); sim_f2 <- list(); sim_f3 <- list()
SIG_sim <- list()

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
  
  ## initial values for MCMC samples
  rho12t <- 0.3
  rho13t <- 0.1
  rho23t <- 0.2
  sig1 <- 1; sig2 <- 0.8;  sig3 <- 0.5;
  Sigmat <- rbind(cbind(sig1, sqrt(sig1)*sqrt(sig2)*rho12t, sqrt(sig1)*sqrt(sig3)*rho13t), 
                  cbind(sqrt(sig1)*sqrt(sig2)*rho12t, sig2, sqrt(sig2)*sqrt(sig3)*rho13t), 
                  cbind(sqrt(sig1)*sqrt(sig3)*rho13t, sqrt(sig1)*sqrt(sig2)*rho12t, sig3))
  
  sig1t <- sig1
  sig2t <- sig2
  sig3t <- sig3
  alp1t <- 0.01 
  alp2t <- 0.01
  alp3t <- 0.01
  tau1t <- alp1t/sig1t
  tau2t <- alp2t/sig2t
  tau3t <- alp3t/sig3t
  f1t <- rep(1, length(Tm)) 
  f2t <- rep(1, length(Tm))
  f3t <- rep(1, length(Tm))
  ## hyper paramter 
  A1 <- 1; A2 <- 1; A3 <- 1;
  B1 <- 1; B2 <- 1; B3 <- 10^3;
  
  
  niter <- 2000 # samples 
  nburn <- 2000 # burning 
  aiter <- niter + nburn # all
  alp1tt <- rep(0, niter)
  alp2tt <- rep(0, niter)
  alp3tt <- rep(0, niter)
  
  f1tt <- matrix(0, niter, length(Tm)) 
  f2tt <- matrix(0, niter, length(Tm))
  f3tt <- matrix(0, niter, length(Tm))
  
  SIG <- NULL
  
  for(k in 1:aiter){
    
    ## identity matrix 
    Im <- diag(1, 100)
    
    ## update f
    ## update f1t
    iV1 <- Im + c(alp1t)*K
    s1 <- svd(iV1)
    
    iV1 <- Im + alp1t*K
    V1 <- solve(iV1)
    
    Mu1 <- V1%*%y1
    Var1 <- c(sqrt(sig1t))*V1
    
    f1t <- matrix(mvrnorm(1, Mu1, Var1), nrow = length(Tm))
    
    ## update f2t
    iV2 <- Im + alp2t*K
    V2 <- solve(iV2)
    
    Mu2 <- V2%*%y2
    Var2 <- sqrt(sig2t)*V2
    
    f2t <- matrix(mvrnorm(1, Mu2, Var2), nrow = length(Tm))
    
    ## update f3t
    iV3 <- Im + alp3t*K
    V3 <- solve(iV3)
    
    Mu3 <- V3%*%y3
    Var3 <- sqrt(sig3t)*V3
    
    f3t <- matrix(mvrnorm(1, Mu3, Var3), nrow = length(Tm))
    ### end of updating f ######
    
    ## update Sigma
    nu <- 5
    S0 <- diag(1, q)
    
    nus <- m + nu
    R1 <- y1 - f1t
    R2 <- y2 - f2t
    R3 <- y3 - f3t
    R <- as.matrix(cbind(R1, R2, R3))
    
    rss <- t(R)%*%R
    Ss <- rss + S0
    
    Sigmac <- rinvwishart(nus, Ss)
    
    ## sig1t
    sig1t <- Sigmac[1,1]
    sig2t <- Sigmac[2,2]
    sig3t <- Sigmac[3,3]
    ### end of Sigma
    
    
    ## update tau(q)
    ## tau1t
    A1star <- 0.5*m + A1
    B1star <- B1 + 0.5*(t(f1t)%*%K%*%f1t)
    tau1t <- rgamma(1, A1star, B1star)
    
    ## tau2t
    A2star <- 0.5*m + A2
    B2star <- B2 + 0.5*(t(f2t)%*%K%*%f2t)
    tau2t <- rgamma(1, A2star, B2star)
    
    ## tau3t
    A3star <- 0.5*m + A3
    B3star <- B3 + 0.5*(t(f3t)%*%K%*%f3t)
    tau3t <- rgamma(1, A3star, B3star)
    ### end of taut
    
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
      
      SIG <- rbind(SIG, c(Sigmac))
      
    }
  }

  sim_alp1[s] <- mean(alp1tt); sim_alp2[s] <- mean(alp2tt); sim_alp3[s] <- mean(alp3tt)
  sim_f1[[s]] <- apply(f1tt, 2, mean)
  sim_f2[[s]] <- apply(f2tt, 2, mean) 
  sim_f3[[s]] <- apply(f3tt, 2, mean)
  SIG_sim[[s]] <- matrix(apply(SIG, 2, mean), q,q)
  
  ## calculate MISE 
  f1i <- apply(f1tt, 2, mean); f2i <- apply(f2tt, 2, mean); f3i <- apply(f3tt, 2, mean)
  X_design <- seq(-2, 2, length = 100)

  dx <- diff(X_design)
  nx <- length(X_design)-1

  Ftrue <- cbind(f1true, f2true, f3true)
  Fhatij <- cbind(f1i,f2i,f3i)
  
  FtF <- (Ftrue - Fhatij)%*%t(Ftrue - Fhatij)
  dFtF <- diag(FtF)
  
  imse_all[s] <- dFtF[1:nx]%*%dx[1:nx]
  
}

Rst.Sim.30 <- list()

Rst.Sim.30[[1]] <- sim_alp1
Rst.Sim.30[[2]] <- sim_alp2
Rst.Sim.30[[3]] <- sim_alp3

Rst.Sim.30[[4]] <- sim_f1
Rst.Sim.30[[5]] <- sim_f2
Rst.Sim.30[[6]] <- sim_f3

Rst.Sim.30[[7]] <- SIG_sim
Rst.Sim.30[[8]] <- imse_all



names(Rst.Sim.30) <- c("alp1", "alp2", "alp3", 
                       "f1", "f2", "f3", 
                       "SIG", "imse")


save(Rst.Sim.30, file = paste(sur_sim, "Rst.Sim.10000.SUR.RData", sep = '/'))



