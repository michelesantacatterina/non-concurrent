rm(list=ls())

library(geex)
library(sandwich)
library(ggplot2)
library(dplyr)

set.seed(1234)

inv_logit <- function(x){
  1/(1 + exp(-x))
}


data_genera <- function(n,
                        tr_effectA, 
                        beta_W_Y=0.25,
                        beta_W_E,
                        beta_E_Y, 
                        beta_W_A,
                        beta_Vk_E = 0,
                        beta_E_Vk,
                        heteroTE,
                        maxE,
                        thresholdA,
                        isDeterministic){
  
  #E          <- sort(runif(n,0,maxE))
  E          <- sort(rnorm(n,0,1))
  W          <- -mean(beta_W_E*E) + beta_W_E*E + rnorm(n,0,1)
  #thresholdA <- runif(1,0.2,0.8)
  thresholdA <- thresholdA*maxE #more or less (1 - proportion of concurrent)
  
  #V is a deterministic function of E
  if(isTRUE(isDeterministic)){
    V_trA      <- ifelse(E<thresholdA,0,1) #this mean that ~100*thresholdA % nonconcurrent vs 100*(1-thresholdA)% concurrent
  }
  
  #V is NOT a deterministic function of E
  if(!isTRUE(isDeterministic)){
    pr_Vk <- (1+exp( -(-thresholdA -mean(E)*beta_E_Vk + E*beta_E_Vk ) ))^-1
    V_trA <- rbinom(n,1,pr_Vk)
  }
  
  A          <- rep(0,n)
  #prA <- 1/2 #complete randomization
  WVtrA <- W[V_trA==1]
  prA        <- (1+exp( -(-mean(WVtrA)*beta_W_A + beta_W_A*WVtrA) ))^-1 #model (a) without enrichment
  #prA        <- (1+exp( -(beta_W_A*W) ))^-1 #model (a) without enrichment
  
  
  #Among those for which treatment A is NOT available;
  #receive control
  A[which(V_trA==0)] <- 0
  
  #Among those for which treatment A is available;
  #sample from binomial with probability prA
  A[which(V_trA==1)] <- rbinom(sum(V_trA),1,prob=prA)
  
  Y0 <- beta_W_Y*W + beta_E_Y*E + beta_Vk_E*V_trA*E + rnorm(n)
  # Y0 <- NULL
  # na0 <- sum(1-V_trA)
  # na1 <- sum(V_trA)
  # Y0[which(V_trA==0)] <- beta_W_Y*W[which(V_trA==0)] + beta_E_Y*E[which(V_trA==0)] + 
  #   beta_Vk_E*V_trA[which(V_trA==0)]*E[which(V_trA==0)] + rnorm(na0)
  # 
  # Y0[which(V_trA==1)] <- beta_W_Y*W[which(V_trA==1)] + beta_E_Y*E[which(V_trA==1)] + 
  #   beta_Vk_E*V_trA[which(V_trA==1)]*E[which(V_trA==1)] + rnorm(na1)
  
  
  if(isTRUE(heteroTE)){
    Y1 <- NULL
    Y1[which(V_trA==0)]   <- Y0[which(V_trA==0)] - tr_effectA #harmful among the non concurrent
    Y1[which(V_trA==1)]   <- Y0[which(V_trA==1)] + tr_effectA #beneficial among the concurrent
  }else{
    # Y1 <- NULL
    # Y1[which(V_trA==0)]  <- Y0[which(V_trA==0)] + tr_effectA # homogeneous treatment effect
    # Y1[which(V_trA==1)]  <- Y0[which(V_trA==1)] + tr_effectA # homogeneous treatment effect
    Y1  <- Y0 + tr_effectA # homogeneous treatment effect
  }
  
  Y               <- Y0
  Y[which(A==1)]  <- Y1[which(A==1)]
  
  data <- data.frame(Y,Y1,Y0,W,E,A,V_trA,thresholdA)
  
  # summary(glm(A~W,data=data,family = "binomial"))
  # data_Vk <- data[data$V_trA==1,]
  # summary(glm(A~W,data=data_Vk,family = "binomial"))
  # 
  # summary(glm(Y~W+E,data=data[data$V_trA==1,]))
  
  return(data)
  
}

################################################################################
#Notes on estimators

# This is based on Notability - Identification Summary notes on Nov 10, 2023
  # - XXX_OC_W:  onlyconcurrent using identification results for cCate(k,W)
  # - XXX_OC_WE: onlyconcurrent using identification results for cCate(k,W,E)
  #   Note: gcomp_OC_W and gcomp_OC_WE is the same estimator but just additionally conditioning on E
  # 
  # - XXX_All_W: allcontrols using identification results for cCate(k,W)
  # - XXX_All_WE: allcontrols using identification results for cCate(k,W,E)

#XXX is 
 # - plugin: same a gcomp but just using plugin and not using EE/sandwich estimator for SE
 # - gcomp: based on EE (parametric)
 # - aipw: based on semiparametric efficiency


################################################################################
#plug in

#Note: XX_OC_W and XX_OC_WE is the same estimator but just additionally conditioning on E
    #OC_W and OC_WE
    plugin_OC <- function(formula_m0,formula_m1,data){
      
      #----------controls
      
        #regress among controls
        m0_only_cont <- glm(formula_m0,data=data[data$A==0,])
        
        #predict among all data
        p0_all_m0_only_cont <- predict(m0_only_cont, newdata = data)
        
        e0 <- mean(p0_all_m0_only_cont) 
        
      #----------treated 
        
        #regress among treated
        m1_only_cont <- glm(formula_m1,data=data[data$A==1,])
        
        #predict among all data
        p1_all_m1_only_cont <- predict(m1_only_cont, newdata = data)
        
        e1 <- mean(p1_all_m1_only_cont) 
      
      return(list(mu0 = e0,
                  mu1 = e1,
                  mu  = e1-e0))
      
    }
    
    #All_W
    plugin_all_W <- function(formula_m0,formula_m1,data){
      
      #----------controls
      
        #regress among all controls 
        m0_all_cont <- glm(formula_m0,data=data[data$A==0,])
        X <- model.matrix(formula_m0, data = data,drop = FALSE)
        W <- model.matrix(Y~W, data = data,drop = FALSE)
        
        #predict among only concurrent
        #xbeta <- (X * data$V_trA) %*% coef(m0_all_cont)
        xbeta <- predict(m0_all_cont, newdata = data_Vk)
        
        #second expectation regressing xbeta onto W among only concurrent
        m2 <- glm(xbeta ~ W,data=data[data$V_trA==1,])
        
        walpha <- (W * data$V_trA) %*% coef(m2)
        
        e0 <- sum(walpha)/sum(data$V_trA)
      
      #----------treated 
      
        #regress among treated - only concurrent data$V_trA==1
        m1_only_cont <- glm(formula_m1,data=data[data$A==1 & data$V_trA==1,])
        
        #predict among all data  - only concurrent data$V_trA==1
        p1_all_m1_only_cont <- predict(m1_only_cont, newdata = data[data$V_trA==1,])
        
        e1 <- mean(p1_all_m1_only_cont) 
        
        return(list(mu0 = e0,
                    mu1 = e1,
                    mu  = e1-e0))
      
    }
    
    #All_WE
    plugin_all_WE <- function(formula_m0,formula_m1,data){
      
      #----------controls
      #regress among controls all controls 
      m0_all_cont <- glm(formula_m0,data=data[data$A==0,])
      
      Y <- data$Y
      X <- model.matrix(formula_m0, data = data,drop = FALSE)
      
      xbeta <- (X) %*% coef(m0_all_cont) # same as predict(m0_only_cont, newdata = data)
      e0 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
      
      #----------treated 
      
      
      #regress among controls only concurrent 
      m1_only_conc <- glm(formula_m1,data=data[data$A==1 & data$V_trA==1,])
      
      Y <- data$Y
      X <- model.matrix(formula_m1, data = data,drop = FALSE) 
      
      xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
      e1 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
      
      return(list(mu0 = e0,
                  mu1 = e1,
                  mu  = e1-e0))
      
    }
    
    
    
## Test simu
# formula_m0_W  <- Y ~ W
# formula_m1_W  <- Y ~ W
# formula_m0_WE <- Y ~ W + E
# formula_m1_WE <- Y ~ W + E
# 
# 
# data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                     beta_W_A,beta_Vk_E,beta_E_Vk,
#                     heteroTE,maxE,thresholdA,isDeterministic)
# data_Vk <- data[data$V_trA==1,]
# 
# #OC_W
# plugin_OC(formula_m0_W,formula_m1_W,data_Vk)
# #OC_WE
# plugin_OC(formula_m0_WE,formula_m1_WE,data_Vk)
# #All_W
# plugin_all_W(formula_m0_WE,formula_m1_WE,data)
# #All_WE
# plugin_all_WE(formula_m0_WE,formula_m1_WE,data)
# 
# 
# 
# point_oc_w <- point_oc_we <- point_all_w <- point_all_we <- true <- NULL
# itera <- 1000
# for(i in 1:itera){
#   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#   data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                       beta_W_A,beta_Vk_E,beta_E_Vk,
#                       heteroTE,maxE,thresholdA,isDeterministic)
# 
#   data$Z <- (data$W*data$E/5 + 0.5)^3
# 
#   data_Vk <- data[data$V_trA==1,]
# 
#   oc_w  <-  plugin_OC(formula_m0_W,formula_m1_W,data_Vk) #IMPORTANT data_VK
#   oc_we <-  plugin_OC(formula_m0_WE,formula_m1_WE,data_Vk) #IMPORTANT data_VK
# 
#   all_w  <- plugin_all_W(formula_m0_WE,formula_m1_WE,data)
#   all_we <- plugin_all_WE(formula_m0_WE,formula_m1_WE,data)
#   
#   point_oc_w[i]   <- oc_w$mu
#   point_oc_we[i]  <- oc_we$mu
#   point_all_w[i]  <- all_w$mu
#   point_all_we[i] <- all_we$mu
#   true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
# 
# }
# 
# mean(point_oc_w);sd(point_oc_w)
# mean(point_oc_we);sd(point_oc_we)
# 
# mean(point_all_w);sd(point_all_w)
# mean(point_all_we);sd(point_all_we)
# 
# 
# sd(point_oc_w)/sd(point_all_w)
# sd(point_oc_we)/sd(point_all_we)
# 
# sd(point_oc_w)/sd(point_all_we)
# sd(point_oc_we)/sd(point_all_w)
    
    
    
################################################################################
#plug in - IPW
    plugin_IPW_OC <- function(formula_me,data){
      
      data_Vk <- data[data$V_trA==1,]
      n <- nrow(data)
      
      #obtain propensities
      fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
      
      #obtain propensities in the Vk population
      eps <- predict(fitTr, type="response", newdata=data_Vk) 
      
      #obtain probability of Vk
      pv <- mean(data$V_trA)
      
      #indicator of Vk = 1
      IVk <- as.numeric(data$V_trA==1)
      
      #obtain phi1
      temp1                   <-  (data_Vk$A/eps)*(data_Vk$Y) 
      temp1_2                 <-  rep(0,n)
      temp1_2[data$V_trA==1]  <-  temp1
      phi1                    <-  (IVk/pv) * temp1_2
      e1 <- mean(phi1)
      
      #obtain phi0
      temp1                   <-  ((1-data_Vk$A)/(1-eps))*(data_Vk$Y) 
      temp1_2                 <-  rep(0,n)
      temp1_2[data$V_trA==1]  <-  temp1
      phi0                    <-  (IVk/pv) * temp1_2
      
      e0 <- mean(phi0)
      
      mu <- mean(phi1-phi0)
      
      return(list(mu0 = e0,
                  mu1 = e1,
                  mu  = e1-e0))
      
    }
    
    
    plugin_IPW_all <- function(formula_me,data){
      
      data_Vk <- data[data$V_trA==1,]
      n <- nrow(data)
      
      #obtain propensities
      fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
  
            
      #obtain propensities in the Vk population
      eps <- predict(fitTr, type="response", newdata=data_Vk) 
      
      #obtain probability of Vk
      pv <- mean(data$V_trA)
      
      #indicator of Vk = 1
      IVk <- as.numeric(data$V_trA==1)
      
      #obtain phi1
      temp1                   <-  (data_Vk$A/eps)*(data_Vk$Y) 
      temp1_2                 <-  rep(0,n)
      temp1_2[data$V_trA==1]  <-  temp1
      phi1                    <-  (IVk/pv) * temp1_2
      e1 <- mean(phi1)
      
      #obtain phi0
      eps0 <- rep(0,n)
      eps0[data$V_trA==1] <- eps
  
    
      temp0 <-  (IVk/pv) *( ((1-data$A)/(1-eps0) )*(data$Y) )

      
      phi0                    <-  temp0
      e0 <- mean(phi0)
      
      mu <- mean(phi1-phi0)
  
      
      return(list(mu0 = e0,
                  mu1 = e1,
                  mu  = e1-e0))
      
    }
    
    
    
# # Test simu
# formula_me  <- A ~ W + E
# 
# 
# data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                     beta_W_A,beta_Vk_E,beta_E_Vk,
#                     heteroTE,maxE,thresholdA,isDeterministic)
# data_Vk <- data[data$V_trA==1,]
# # 
# # #OC_W
# plugin_IPW_OC(A ~ W + E,data_Vk)
# # #OC_WE
# # plugin_IPW_OC(A ~ W ,data_Vk)
# # #all_W
# plugin_IPW_all(A ~ W + E,data)
# #all_WE
# plugin_IPW_all(A ~ W ,data)
# 
# 
# 
# point_oc_w <- point_oc_we <- point_all_w <- point_all_we <- true <- NULL
# itera <- 1000
# for(i in 1:itera){
#   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#   data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                       beta_W_A,beta_Vk_E,beta_E_Vk,
#                       heteroTE,maxE,thresholdA,isDeterministic)
# 
#   data$Z <- (data$W*data$E/5 + 0.5)^3
# 
#   data_Vk <- data[data$V_trA==1,]
# 
#   oc_w  <-  plugin_IPW_OC(A ~ W,data_Vk) #IMPORTANT data_VK
#   oc_we <-  plugin_IPW_OC(A ~ W + E,data_Vk) #IMPORTANT data_VK
# 
#   all_w  <- plugin_IPW_all(A ~ W,data)
#   all_we <- plugin_IPW_all(A ~ W + E ,data)
# 
#   point_oc_w[i]   <- oc_w$mu
#   point_oc_we[i]  <- oc_we$mu
#   point_all_w[i]  <- all_w$mu
#   point_all_we[i] <- all_we$mu
#   true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
# 
# }
# 
# mean(point_oc_w);sd(point_oc_w)
# mean(point_oc_we);sd(point_oc_we)
# 
# mean(point_all_w);sd(point_all_w)
# mean(point_all_we);sd(point_all_we)
# 
# 
# sd(point_oc_w)/sd(point_all_w)
# sd(point_oc_we)/sd(point_all_we)
# 
# sd(point_oc_w)/sd(point_all_we)
# sd(point_oc_we)/sd(point_all_w)
    


################################################################################
#Gcomp - EIF

#--------------------Concurrent only 
    gcomp_OC <- function(formula_m0,formula_m1,data){
      
      #----------controls
      
          #regress among controls only concurrent 
          m0_only_conc <- glm(formula_m0,data=data[data$A==0 & data$V_trA==1,])
          
          Y <- data$Y
          X <- model.matrix(formula_m0, data = data,drop = FALSE)
          
          xbeta <- (X) %*% coef(m0_only_conc) # same as predict(m0_only_conc, newdata = data)
          e0 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
          
          h0 <- X * data$V_trA * (1-data$A) * as.numeric(t(Y -  xbeta) )
          h1 <- data$V_trA*(xbeta - e0)
          
          #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
          g11 <- ( t(X*data$V_trA*(1-data$A)) %*% X) / nrow(data)
          g11inv <- solve(g11)
          
          g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
          
          g22 <- sum(data$V_trA)
          g22inv <- 1/g22
          
          phi_beta_hat <- g11inv %*% t(h0)
          # B <- t(h0) %*% h0
          # A <- g11inv
          # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
          # same as 
          # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
          var_beta_hat0 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
          
          #phi_mu0 <- n*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
          phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
          #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
          se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
          
      
      #----------treated 
          
          
          #regress among controls only concurrent 
          m1_only_conc <- glm(formula_m1,data=data[data$A==1 & data$V_trA==1,])
          
          Y <- data$Y
          X <- model.matrix(formula_m1, data = data,drop = FALSE)
          
          xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
          e1 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
          
          h0 <- X * data$V_trA * (data$A) * as.numeric(t(Y -  xbeta) )
          h1 <- data$V_trA*(xbeta - e1)
          
          #g11 <- ( t(X0) %*% diag(data$A) %*% X0) / nrow(data)
          g11 <- ( t(X*data$V_trA*(data$A)) %*% X) / nrow(data)
          g11inv <- solve(g11)
          
          g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
          
          g22 <- sum(data$V_trA)
          g22inv <- 1/g22
          
          phi_beta_hat <- g11inv %*% t(h0)
          # B <- t(h0) %*% h0
          # A <- g11inv
          # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
          # same as 
          # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
          var_beta_hat1 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
          
          #phi_mu1 <- n*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
          phi_mu1 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
          #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
          se1 <- sqrt(var(phi_mu1)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
          
          
      #concATE
          
          mu    <- e1 - e0
          se_if <- sqrt(var(phi_mu1-phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      return(list(mu0 = e0,
                  se0_if = se0,
                  mu1 = e1,
                  se1_if = se1,
                  mu = mu,
                  se_if = se_if))
      
    }
    
    #--------------------All controls
    gcomp_all_W <- function(formula_m0,formula_m1,data){
      
      
      #----------controls
      
      #regress among controls all controls 
      m0_all_cont <- glm(formula_m0,data=data[data$A==0,])
      
      Y <- data$Y
      X <- model.matrix(formula_m0, data = data,drop = FALSE)
      W <- model.matrix(Y~W, data = data,drop = FALSE)
      
      #predict among all controls
      xbeta <- (X) %*% coef(m0_all_cont)
      #xbeta <- predict(m0_all_cont, newdata = data_Vk)
      
      #second expectation regressing xbeta onto W among only concurrent
      m2 <- glm(xbeta[data$V_trA==1,] ~ W,data=data[data$V_trA==1,])
      
      walpha <- (W) %*% coef(m2)
      
      e0 <- sum(walpha* data$V_trA)/sum(data$V_trA)
      
      h1 <- X * (1-data$A) * as.numeric(t(Y -  xbeta) )
      h2 <- W * data$V_trA * as.numeric(t(xbeta -  walpha) )
      h3 <- data$V_trA*(walpha - e0)
      
      #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
      g11    <- ( t(X*(1-data$A)) %*% X) / nrow(data)
      g11inv <- solve(g11)
      
      g21 <- - (t(W)%*% diag(data$V_trA) %*%X)/nrow(data)
      g22 <- (t(W)%*% diag(data$V_trA) %*%W)/nrow(data)
      g22inv <- solve(g22)
      
      g32 <- apply( - W*data$V_trA,2,sum)/nrow(data) #-data$V_trA%*%W/nrow(data)
      g33 <- sum(data$V_trA)
      g33inv <- 1/g33
      
      phi_beta_hat <- g11inv %*% t(h1)
      # B <- t(h1) %*% h1
      # A <- g11inv
      # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_all_cont))
      # #same as
      # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      varcov_beta_hat <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      
      phi_alpha_hat <- g22inv %*% t((h2 - t( (g21) %*% (phi_beta_hat) ) ))
      
      varcov_alpha_hat <- phi_alpha_hat %*% t(phi_alpha_hat)/nrow(data)^2 
      
      #Sandwich estimator does not consider first step estimation
      # aa <- g22inv %*% t(h2)
      # aa %*% t(aa)/nrow(data)^2
      # diag(sandwich(m2))
      
      phi_mu0 <- as.numeric(n*g33inv %*% t((h3 - t( (g32) %*% (phi_alpha_hat) ) )))
      #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
      se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      
      #----------treated 
      
      
      #regress among controls only concurrent 
      m1_only_conc <- glm(formula_m1,data=data[data$A==1 & data$V_trA==1,])
      
      Y <- data$Y
      X <- model.matrix(formula_m1, data = data,drop = FALSE)
      
      xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
      e1 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
      
      h0 <- X * data$V_trA * (data$A) * as.numeric(t(Y -  xbeta) )
      h1 <- data$V_trA*(xbeta - e1)
      
      #g11 <- ( t(X0) %*% diag(data$A) %*% X0) / nrow(data)
      g11 <- ( t(X*data$V_trA*(data$A)) %*% X) / nrow(data)
      g11inv <- solve(g11)
      
      g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
      
      g22 <- sum(data$V_trA)
      g22inv <- 1/g22
      
      phi_beta_hat <- g11inv %*% t(h0)
      # B <- t(h0) %*% h0
      # A <- g11inv
      # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
      # same as 
      # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      var_beta_hat1 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      
      phi_mu1 <- n*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
      #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
      se1 <- sqrt(var(phi_mu1)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      
      #----------concATE
      
      mu    <- e1 - e0
      se_if <- sqrt(var(phi_mu1-phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      
      #----------return
      return(list(mu0 = e0,
                  se0_if = se0,
                  mu1 = e1,
                  se1_if = se1,
                  mu = mu,
                  se_if = se_if))
      
    }    

#--------------------All controls - without regressing walpha
    gcomp_all_WE <- function(formula_m0,formula_m1,data){
      
      #----------controls
      
      #regress among controls all controls 
      m0_all_cont <- glm(formula_m0,data=data[data$A==0,])
      
      Y <- data$Y
      X <- model.matrix(formula_m0, data = data,drop = FALSE)
      
      xbeta <- (X) %*% coef(m0_all_cont) # same as predict(m0_only_cont, newdata = data)
      e0 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
      
      h0 <- X * (1-data$A) * as.numeric(t(Y -  xbeta) ) 
      h1 <- data$V_trA*(xbeta - e0)
      
      #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
      g11 <- ( t(X*(1-data$A)) %*% X) / nrow(data)
      g11inv <- solve(g11)
      
      g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
      
      g22 <- sum(data$V_trA)
      g22inv <- 1/g22
      
      phi_beta_hat <- g11inv %*% t(h0)
      # B <- t(h0) %*% h0
      # A <- g11inv
      # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_all_cont))
      # # same as
      # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
      var_beta_hat0 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      
      phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
      #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
      se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      #----------treated 
      
      
      #regress among controls only concurrent 
      m1_only_conc <- glm(formula_m1,data=data[data$A==1 & data$V_trA==1,])
      
      Y <- data$Y
      X <- model.matrix(formula_m1, data = data,drop = FALSE)
      
      xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
      e1 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
      
      h0 <- X * data$V_trA * (data$A) * as.numeric(t(Y -  xbeta) )
      h1 <- data$V_trA*(xbeta - e1)
      
      #g11 <- ( t(X0) %*% diag(data$A) %*% X0) / nrow(data)
      g11 <- ( t(X*data$V_trA*(data$A)) %*% X) / nrow(data)
      g11inv <- solve(g11)
      
      g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
      
      g22 <- sum(data$V_trA)
      g22inv <- 1/g22
      
      phi_beta_hat <- g11inv %*% t(h0)
      # B <- t(h0) %*% h0
      # A <- g11inv
      # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
      # same as 
      # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      var_beta_hat1 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
      
      phi_mu1 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
      #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
      se1 <- sqrt(var(phi_mu1)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      #concATE
      
      mu    <- e1 - e0
      se_if <- sqrt(var(phi_mu1-phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
      
      return(list(mu0 = e0,
                  se0_if = se0,
                  mu1 = e1,
                  se1_if = se1,
                  mu = mu,
                  se_if = se_if))
      
    }
    


# ## Test simu
# formula_m0_W  <- Y ~ W
# formula_m1_W  <- Y ~ W
# formula_m0_WE <- Y ~ W + E
# formula_m1_WE <- Y ~ W + E
# 
# 
# data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                     beta_W_A,beta_Vk_E,beta_E_Vk,
#                     heteroTE,maxE,thresholdA,isDeterministic)
# data_Vk <- data[data$V_trA==1,]
# 
# #OC_W
# gcomp_OC(formula_m0_W,formula_m1_W,data_Vk)
# #OC_WE
# gcomp_OC(formula_m0_WE,formula_m1_WE,data_Vk)
# #All_W
# gcomp_all_W(formula_m0_WE,formula_m1_WE,data)
# #All_WE
# gcomp_all_WE(formula_m0_WE,formula_m1_WE,data)
# 
# 
# 
# point_oc_w <- point_oc_we <- point_all_w <- point_all_we <- true <- NULL
# se_oc_w <- se_oc_we <- se_all_w <- se_all_we <- NULL
# itera <- 1000
# for(i in 1:itera){
#   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#   data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                       beta_W_A,beta_Vk_E,beta_E_Vk,
#                       heteroTE,maxE,thresholdA,isDeterministic)
# 
#   data$Z <- (data$W*data$E/5 + 0.5)^3
# 
#   data_Vk <- data[data$V_trA==1,]
# 
#   oc_w  <-  gcomp_OC(formula_m0_W,formula_m1_W,data_Vk) #IMPORTANT data_VK
#   oc_we <-  gcomp_OC(formula_m0_WE,formula_m1_WE,data_Vk) #IMPORTANT data_VK
# 
#   all_w  <- gcomp_all_W(formula_m0_WE,formula_m1_WE,data)
#   all_we <- gcomp_all_WE(formula_m0_WE,formula_m1_WE,data)
# 
#   point_oc_w[i]   <- oc_w$mu
#   point_oc_we[i]  <- oc_we$mu
#   point_all_w[i]  <- all_w$mu
#   point_all_we[i] <- all_we$mu
#   
#   se_oc_w[i]   <- oc_w$se_if
#   se_oc_we[i]  <- oc_we$se_if
#   se_all_w[i]  <- all_w$se_if
#   se_all_we[i] <- all_we$se_if
#   true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
# 
# }
# 
# mean(point_oc_w);sd(point_oc_w);mean(se_oc_w)
# mean(point_oc_we);sd(point_oc_we);mean(se_oc_we)
# 
# mean(point_all_w);sd(point_all_w);mean(se_all_w)
# mean(point_all_we);sd(point_all_we);mean(se_all_we)
# 
# 
# sd(point_oc_w)/sd(point_all_w)
# sd(point_oc_we)/sd(point_all_we)
# 
# sd(point_oc_w)/sd(point_all_we)
# sd(point_oc_we)/sd(point_all_w)


    
    
    
################################################################################
#IPW - EIF
    
#--------------------Concurrent only 
IPW_OC <- function(formula_me,data){
  
  #regress among concurrent
  m0_only_conc <- glm(formula_me,data=data[data$V_trA==1,],family="binomial")
  X <- model.matrix(formula_me, data = data,drop = FALSE)
  
  #propensity score among concurrent
  xbeta <- (X) %*% coef(m0_only_conc) 
  exp_xbeta <- as.numeric(exp(xbeta))
  exp_m_xbeta <- exp(-xbeta)
  ps <- as.numeric(exp_xbeta/(1+exp_xbeta)) #same as predict(m0_only_conc, newdata = data, type="response")

  pv <- mean(data$V_trA)
  #----------Controls
  
    w0 <- (1-data$V_trA*data$A)/(1-ps)
    #e0 <- sum(w0*data$V_trA*data$Y)/sum(w0*data$V_trA) #mean(data_Vk$Y0)
    e0 <- mean( (data$V_trA/pv) *( ( w0 )*(data$Y) ) )
    
    h0 <- X * data$V_trA * ( (1-data$A) - (1-as.numeric(t(ps))) )
    
    g11 <- ( t(X * data$V_trA * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
    g11inv <- solve(g11)  
      
    phi_beta_hat <- g11inv %*% t(h0)
    # B <- t(h0) %*% h0
    # A <- g11inv
    # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
    # same as 
    # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
    
    
    h1 <- data$V_trA*(w0*data$Y - e0)
    
    a0 <- 1-data$A
    #g21 <- as.numeric((- t( a0 * data$V_trA * data$Y *exp_xbeta ) %*% X)/nrow(data))
    g21 <- as.numeric((t( a0 * data$V_trA * data$Y *exp_xbeta ) %*% X)/nrow(data))
    
    #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
    
    g22 <- sum(data$V_trA)
    g22inv <- 1/g22
    
    
    phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
    se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  
    
  #----------Treated
    
    w1 <- (data$V_trA*data$A)/(ps)
    #e1 <- sum(w1*data$V_trA*data$Y)/sum(w1*data$V_trA) #mean(data_Vk$Y1)
    e1 <- mean( (data$V_trA/pv) *( ( w1 )*(data$Y) ) )
    
    h0 <- X * data$V_trA * ( data$A - as.numeric(t(ps)) )
    
    g11 <- ( t(X * data$V_trA * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
    g11inv <- solve(g11)  
    
    phi_beta_hat <- g11inv %*% t(h0)
    # B <- t(h0) %*% h0
    # A <- g11inv
    # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
    # same as 
    # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
    
    
    h1 <- data$V_trA*(w1*data$Y - e1)
    
    a1 <- data$A
    g21 <- as.numeric((t( a1 * data$V_trA * data$Y *exp_m_xbeta ) %*% X)/nrow(data))
    
    #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
    
    g22 <- sum(data$V_trA)
    g22inv <- 1/g22
    
    
    phi_mu1 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
    se1 <- sqrt(var(phi_mu1)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  #concATE
  
  mu    <- e1 - e0
  se_if <- sqrt(var(phi_mu1-phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  return(list(mu0 = e0,
              se0_if = se0,
              mu1 = e1,
              se1_if = se1,
              mu = mu,
              se_if = se_if))
  
}
    
    
    
#--------------------All controls
IPW_all <- function(formula_me,data){
  
  #regress among concurrent
  m0_only_conc <- glm(formula_me,data=data[data$V_trA==1,],family="binomial")
  X <- model.matrix(formula_me, data = data,drop = FALSE)
  
  #propensity score among concurrent
  xbeta <- (X) %*% coef(m0_only_conc) 
  exp_xbeta <- as.numeric(exp(xbeta))
  exp_m_xbeta <- exp(-xbeta)
  ps <- as.numeric(exp_xbeta/(1+exp_xbeta)) #same as predict(m0_only_conc, newdata = data, type="response")
  ps0 <- rep(0,nrow(data))
  ps0[data$V_trA==1] <- ps[data$V_trA==1]
  
  pv <- mean(data$V_trA)
  
  #----------Controls
  
  #w0 <- (1-data$V_trA*data$A)/(1-ps)
  w0 <- (1-data$A)/(1-ps0)
  #e0 <- sum(w0*data$V_trA*data$Y)/sum(w0*data$V_trA) #mean(data_Vk$Y0)
  e0 <- mean( (data$V_trA/pv) *( ( w0 )*(data$Y) ) )
  
  #h0 <- X * data$V_trA * ( (1-data$A) - (1-as.numeric(t(ps))) )
  h0 <- X * data$V_trA * ( (1-data$A) - (1-as.numeric(t(ps0))) )
  
  g11 <- ( t(X * data$V_trA * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
  g11inv <- solve(g11)  
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
  # # same as
  # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
  
  
  h1 <- (data$V_trA)*((w0)*data$Y - e0)
  
  a0 <- 1-data$A
  #g21 <- as.numeric((- t( a0 * data$V_trA * data$Y *exp_xbeta ) %*% X)/nrow(data))
  g21 <- as.numeric((t( a0 * data$V_trA * data$Y *exp_xbeta ) %*% X)/nrow(data))
  
  #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
  
  g22 <- sum(data$V_trA)
  g22inv <- 1/g22
  
  
  phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
  se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  
  
  #----------Treated
  
  w1 <- (data$V_trA*data$A)/(ps)
  #e1 <- sum(w1*data$V_trA*data$Y)/sum(w1*data$V_trA) #mean(data_Vk$Y1)
  e1 <- mean( (data$V_trA/pv) *( ( w1 )*(data$Y) ) )
  
  h0 <- X * data$V_trA * ( data$A - as.numeric(t(ps)) )
  
  g11 <- ( t(X * data$V_trA * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
  g11inv <- solve(g11)  
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
  # same as 
  # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
  
  
  h1 <- data$V_trA*(w1*data$Y - e1)
  
  a1 <- data$A
  g21 <- as.numeric((t( a1 * data$V_trA * data$Y *exp_m_xbeta ) %*% X)/nrow(data))
  
  #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
  
  g22 <- sum(data$V_trA)
  g22inv <- 1/g22
  
  
  phi_mu1 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
  se1 <- sqrt(var(phi_mu1)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  #concATE
  
  mu    <- e1 - e0
  se_if <- sqrt(var(phi_mu1-phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  return(list(mu0 = e0,
              se0_if = se0,
              mu1 = e1,
              se1_if = se1,
              mu = mu,
              se_if = se_if))
  
}
    
    
    ## Test simu
    formula_me  <- A ~ W + E


    # data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
    #                     beta_W_A,beta_Vk_E,beta_E_Vk,
    #                     heteroTE,maxE,thresholdA,isDeterministic)
    # data_Vk <- data[data$V_trA==1,]
    # 
    # #OC_W
    # IPW_OC(A ~ W + E,data_Vk)
    # #OC_WE
    # IPW_all(A ~ W + E,data)
    # 
    # 
    # 
    # 
    # point_oc_w <- point_oc_we <- point_all_w <- point_all_we <- true <- NULL
    # se_oc_w <- se0_oc_w <- se1_oc_w <- se_oc_we <- se_all_we <- se0_oc_we <- se1_oc_we <- NULL
    # point0_oc_w <- point0_oc_we <- point1_oc_w <- point1_oc_we <- point1_all_we <- NULL
    # itera <- 1000
    # for(i in 1:itera){
    #   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
    #   data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
    #                       beta_W_A,beta_Vk_E,beta_E_Vk,
    #                       heteroTE,maxE,thresholdA,isDeterministic)
    # 
    #   data$Z <- (data$W*data$E/5 + 0.5)^3
    # 
    #   data_Vk <- data[data$V_trA==1,]
    # 
    #   oc_w  <-  IPW_OC(A ~ W,data_Vk) #IMPORTANT data_VK
    #   oc_we <-  IPW_OC(A ~ W + E,data_Vk) #IMPORTANT data_VK
    #   
    #   all_we  <-  IPW_all(A ~ W + E,data) #IMPORTANT data_VK
    #   
    # 
    #   point_oc_w[i]   <- oc_w$mu
    #   point_oc_we[i]  <- oc_we$mu
    #   point_all_we[i]  <- all_we$mu
    # 
    #   point0_oc_w[i]   <- oc_w$mu0
    #   point0_oc_we[i]  <- oc_we$mu0
    # 
    #   point1_oc_w[i]   <- oc_w$mu1
    #   point1_oc_we[i]  <- oc_we$mu1
    # 
    #   se_oc_w[i] <- oc_w$se_if
    #   se_oc_we[i] <- oc_we$se_if
    #   se_all_we[i] <- all_we$se_if
    # 
    #   se0_oc_w[i] <- oc_w$se0_if
    #   se0_oc_we[i] <- oc_we$se0_if
    # 
    #   se1_oc_w[i] <- oc_w$se1_if
    #   se1_oc_we[i] <- oc_we$se1_if
    # 
    #   true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
    # 
    # }
    # 
    # mean(point_oc_we);sd(point_oc_we);mean(se_oc_we);
    # mean(point_all_we);sd(point_all_we);mean(se_all_we);
    
    # mean(point_oc_w);sd(point_oc_w);mean(se_oc_w)
    # mean(point0_oc_w);sd(point0_oc_w);mean(se0_oc_w);
    # mean(point0_oc_we);sd(point0_oc_we);mean(se0_oc_we);
    # mean(point1_oc_w);sd(point1_oc_w);mean(se1_oc_w);
    # mean(point1_oc_we);sd(point1_oc_we);mean(se1_oc_we);
    # 
    # 
    # sd(point_oc_w)/sd(point_all_w)
    # sd(point_oc_we)/sd(point_all_we)
    # 
    # sd(point_oc_w)/sd(point_all_we)
    # sd(point_oc_we)/sd(point_all_w)
    

################################################################################


#-----------AIPW conc
aipw_OC <- function(formula_m0,formula_m1,formula_me,data){
  
  data_Vk <- data[data$V_trA==1,]
  n <- nrow(data)
  
  #regress among only concurrent controls
  m0_only_concC <- glm(formula_m0,data=data_Vk[data_Vk$A==0,])
  
  #predict among only concurrent controls
  p0_only_concC_m0_only_concC <- predict(m0_only_concC, newdata = data_Vk)
  
  #regress among only concurrent treated
  mk_only_concC <- glm(formula_m1,data=data_Vk[data_Vk$A==1,])
  
  #predict among only concurrent treated
  pk_only_concC_mk_only_concC <- predict(mk_only_concC, newdata = data_Vk)
  
  #obtain propensities
  fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
  
  #obtain propensities in the Vk population
  eps <- predict(fitTr, type="response", newdata=data_Vk) 
  
  #obtain probability of Vk
  pv <- mean(data$V_trA)
  
  #indicator of Vk = 1
  IVk <- as.numeric(data$V_trA==1)
  
  
  #obtain phi1
  temp1                   <- (data_Vk$A/eps)*(data_Vk$Y-pk_only_concC_mk_only_concC) + pk_only_concC_mk_only_concC
  temp1_2                 <- rep(0,n)
  temp1_2[data$V_trA==1]  <-  temp1
  phi1                    <-  (IVk/pv) * temp1_2
  se1 <- (sqrt(var(phi1)/n)) * sqrt((n-1)/n)
  e1 <- mean(phi1)
  
  #obtain phi0
  temp0                   <- ( (1-data_Vk$A)/(1-eps) )*(data_Vk$Y-p0_only_concC_m0_only_concC) + p0_only_concC_m0_only_concC
  temp0_2                 <- rep(0,n)
  temp0_2[data$V_trA==1]  <-  temp0
  phi0                    <- (IVk/pv) * temp0_2
  se0 <- (sqrt(var(phi0)/n)) * sqrt((n-1)/n)
  e0 <- mean(phi0)
  
  mu <- mean(phi1-phi0)
  se_if <- (sqrt(var(phi1-phi0)/n)) * sqrt((n-1)/n)
  
  #----------return
  return(list(mu0 = e0,
              se0_if = se0,
              mu1 = e1,
              se1_if = se1,
              mu = mu,
              se_if = se_if))
  
  
}


#-----------AIPW all
aipw_all <- function(formula_m0,formula_m1,formula_me,data,nonpara){
  
  data_Vk <- data[data$V_trA==1,]
  n <- nrow(data)
  
  #regress among all controls
  m0_allC <- glm(formula_m0,data=data[data$A==0,])
  Y <- data$Y
  X <- model.matrix(formula_m0, data = data,drop = FALSE)
  # 
  #predict among only concurrent controls
  #p0_only_concC_m0_allC <- predict(m0_allC, newdata = data_Vk)
  p0_only_concC_m0_allC <- (X) %*% coef(m0_allC)
  
  #regress among only concurrent treated
  mk_only_concC <- glm(formula_m1,data=data_Vk[data_Vk$A==1,])
  
  #predict among only concurrent treated
  pk_only_concC_mk_only_concC <- predict(mk_only_concC, newdata = data_Vk)
  
  #obtain propensities
  fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
  
  #obtain propensities in the Vk population
  eps <- predict(fitTr, type="response", newdata=data_Vk) 
  
  #obtain probability of Vk
  pv <- mean(data$V_trA)
  
  #indicator of Vk = 1
  IVk <- as.numeric(data$V_trA==1)
  
  #obtain phi1
  temp1                   <-  (data_Vk$A/eps)*(data_Vk$Y-pk_only_concC_mk_only_concC) + pk_only_concC_mk_only_concC
  temp1_2                 <-  rep(0,n)
  temp1_2[data$V_trA==1]  <-  temp1
  phi1                    <-  (IVk/pv) * temp1_2
  se1 <- (sqrt(var(phi1)/n)) * sqrt((n-1)/n)
  e1 <- mean(phi1)
  
  #obtain phi0
  # temp0                   <-  ( (1-data_Vk$A)/(1-eps) )*(data_Vk$Y-p0_only_concC_m0_allC) + p0_only_concC_m0_allC
  eps0 <- rep(0,n)
  eps0[data$V_trA==1] <- eps
  #temp0                   <-  ( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
  #temp0                   <-  (IVk/pv) *( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + p0_only_concC_m0_allC )
  # temp0_2                 <-  rep(0,n)
  # temp0_2[data$V_trA==1]  <-  temp0
  
  if(isTRUE(nonpara)){
    temp0 <-  (IVk/pv) *( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + p0_only_concC_m0_allC )
  }else{
    temp0 <- ( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
  }
  
  phi0                    <-  temp0
  se0 <- (sqrt(var(phi0)/n)) * sqrt((n-1)/n)
  e0 <- mean(phi0)
  
  mu <- mean(phi1-phi0)
  se_if <- (sqrt(var(phi1-phi0)/n)) * sqrt((n-1)/n)
  
  #----------return
  return(list(mu0 = e0,
              se0_if = se0,
              mu1 = e1,
              se1_if = se1,
              mu = mu,
              se_if = se_if))
  
  
}



# ## Test simu
# formula_m0_W  <- Y ~ W
# formula_m1_W  <- Y ~ W
# formula_m0_WE <- Y ~ W + E
# formula_m1_WE <- Y ~ W + E
# formula_me <- A ~ W 
# 
# data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                     beta_W_A,beta_Vk_E,beta_E_Vk,
#                     heteroTE,maxE,thresholdA,isDeterministic)
# data_Vk <- data[data$V_trA==1,]
# 
# #OC_W
# aipw_OC(formula_m0_W,formula_m1_W,formula_me,data_Vk)
# #OC_WE
# aipw_OC(formula_m0_WE,formula_m1_WE,formula_me,data_Vk)
# #All_WE
# aipw_all(formula_m0_WE,formula_m1_WE,formula_me,data,nonpara=T)
# 
#
# 
# 
# point_oc_w <- point_oc_we <- point_all_w <- point_all_we <- true <- NULL
# se_oc_w <- se_oc_we <- se_all_w <- se_all_we <- NULL
# itera <- 1000
# for(i in 1:itera){
#   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#   data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                       beta_W_A,beta_Vk_E,beta_E_Vk,
#                       heteroTE,maxE,thresholdA,isDeterministic)
#   
#   data$Z <- (data$W*data$E/5 + 0.5)^3
#   
#   data_Vk <- data[data$V_trA==1,]
#   
#   oc_w  <-  aipw_OC(formula_m0_W,formula_m1_W,formula_me,data_Vk) #IMPORTANT data_VK
#   oc_we <-  aipw_OC(formula_m0_WE,formula_m1_WE,formula_me,data_Vk) #IMPORTANT data_VK
#   
#   all_w  <- aipw_all(formula_m0_WE,formula_m1_WE,formula_me,data,nonpara=T)
#   
#   point_oc_w[i]   <- oc_w$mu
#   point_oc_we[i]  <- oc_we$mu
#   point_all_w[i]  <- all_w$mu
#   
#   se_oc_w[i]   <- oc_w$se_if
#   se_oc_we[i]  <- oc_we$se_if
#   se_all_w[i]  <- all_w$se_if
#   true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
#   
# }
# 
# mean(point_oc_w);sd(point_oc_w);mean(se_oc_w)
# mean(point_oc_we);sd(point_oc_we);mean(se_oc_we)
# 
# mean(point_all_w);sd(point_all_w);mean(se_all_w)
# 
# 
# sd(point_oc_w)/sd(point_all_w)
# sd(point_oc_we)/sd(point_all_w)






################################################################################
################################################################################
# Simulations
################################################################################
################################################################################





#COMMENT: We neet to find a way to generate data for which the effect of W on Y is still there in Vk=1 population

# to evaluate estimators: , where b is bias^2 and m is mse
# bias/se
ratio = function(b,m) { sqrt(b/(m-b)) }


# Generate data
tr_effectA      <- 0.8    #large treatment effect 
beta_E_Y        <- 0.5    #effect size E->Y

beta_Vk_E       <- 0      #no interaction Vk/E 
maxE            <- 1      #max entry time (standardized between 0 and maxE)
heteroTE        <- FALSE  #homogeneous treatment effect
thresholdA      <- 0.6   #40% concurrent ::: 100*(1-thresholdA)% concurrent
beta_W_A        <- 0.8 + 0*(thresholdA)*maxE   #effect size W->A 
beta_W_Y        <- 0.8 + 0*(thresholdA)*maxE  #effect size W->Y
beta_W_E        <- 0.8 + 0*(thresholdA)*maxE

isDeterministic <- TRUE

n <- 100
itera <- 1000



# #Correct-Outcome Correct-Treatement
# formula_m0_W  <- Y ~ W
# formula_m1_W  <- Y ~ W
# formula_m0_WE <- Y ~ W + E
# formula_m1_WE <- Y ~ W + E
# formula_me    <- A ~ W + E
# formula_me_W  <- A ~ W
# formula_me_WE <- A ~ W + E
# formula_mv    <- V_trA ~ W + E




# #Miss-Outcome Correct-Treatement
formula_m0_W  <- Y ~ 1
formula_m1_W  <- Y ~ 1
formula_m0_WE <- Y ~ 1
formula_m1_WE <- Y ~ 1
formula_me    <- A ~ W + E
formula_me_W  <- A ~ W
formula_me_WE <- A ~ W + E
# 
# 
# #Correct-Outcome Miss-Treatement
# formula_m0_W  <- Y ~ W
# formula_m1_W  <- Y ~ W
# formula_m0_WE <- Y ~ W + E
# formula_m1_WE <- Y ~ W + E
# formula_me    <- A ~ 1
# formula_me_W  <- A ~ 1
# formula_me_WE <- A ~ 1
# 
# 
# #Miss-Outcome Miss-Treatement
# formula_m0_W  <- Y ~ 1
# formula_m1_W  <- Y ~ 1
# formula_m0_WE <- Y ~ 1
# formula_m1_WE <- Y ~ 1
# formula_me    <- A ~ 1
# formula_me_W  <- A ~ 1
# formula_me_WE <- A ~ 1


# data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
#                     beta_W_A,beta_Vk_E,beta_E_Vk,
#                     heteroTE,maxE,thresholdA,isDeterministic)
# data_Vk <- data[data$V_trA==1,]


gcomp_point_oc_w <- gcomp_point_oc_we <- gcomp_point_all_w <- gcomp_point_all_we <- true <- NULL
aipw_point_oc_w <- aipw_point_oc_we <- aipw_point_all <- NULL
gcomp_se_oc_w <- gcomp_se_oc_we <- gcomp_se_all_w <- gcomp_se_all_we <- NULL
aipw_se_oc_w <- aipw_se_oc_we <- aipw_se_all <- NULL

IPW_point_oc_w <- IPW_point_oc_we <- ipw_se_oc_w <- ipw_se_oc_we <- NULL

se_oc_w <- se_oc_we <- se_all_w <- se_all_we <- NULL


mse_gcomp_point_oc_w <- mse_gcomp_point_oc_we <- mse_gcomp_point_all_we <- 
  mse_aipw_point_oc_w <- mse_aipw_point_oc_we <- mse_aipw_point_all_we <- 
  mse_ipw_point_oc_w <- mse_ipw_point_oc_we <- NULL

bias_gcomp_point_oc_w <- bias_gcomp_point_oc_we <- bias_gcomp_point_all_we <- 
  bias_aipw_point_oc_w <- bias_aipw_point_oc_we <- bias_aipw_point_all_we <-     
  bias_ipw_point_oc_w <- bias_ipw_point_oc_we <- NULL

se_ipw_point_oc_w_sv <-  se_ipw_point_oc_we_sv <- 
  se_ipw_point_oc_w_if <-  se_ipw_point_oc_we_if <- NULL

se_gcomp_point_oc_w_sv <- se_gcomp_point_oc_we_sv <- se_gcomp_point_all_we_sv <- 
  se_gcomp_point_oc_w_if <- se_gcomp_point_oc_we_if <- se_gcomp_point_all_we_if <- NULL

se_aipw_point_oc_w_sv <- se_aipw_point_oc_we_sv <- se_aipw_point_all_we_sv <- 
  se_aipw_point_oc_w_if <- se_aipw_point_oc_we_if <- se_aipw_point_all_we_if <- NULL

ratio_gcomp_oc_w_vs_all <- ratio_gcomp_oc_we_vs_all <- ratio_aipw_oc_w_vs_all <- 
  ratio_aipw_oc_we_vs_all <- NULL


M_mse_gcomp_point_oc_w <- M_mse_gcomp_point_oc_we <- M_mse_gcomp_point_all_we <-
  M_mse_aipw_point_oc_w <- M_mse_aipw_point_oc_we <- M_mse_aipw_point_all_we <- 
  M_mse_ipw_point_oc_w <- M_mse_ipw_point_oc_we <- NULL



coverage_gcompOR_w <- coverage_gcompOR_we <- coverage_gcompall_w <- 
  coverage_gcompall_we <- coverage_aipw_point_oc_w <- coverage_aipw_point_oc_we <-
  coverage_aipw_point_all_w <- coverage_ipw_point_oc_w <- coverage_ipw_point_oc_we <- NULL


M_coverage_gcompOR_w <- M_coverage_gcompOR_we <- M_coverage_gcomp_all_w <-
  M_coverage_gcomall_we <- M_coverage_aipw_oc_w <- M_coverage_aipw_oc_we <- 
  M_coverage_aipw_all <- M_coverage_ipw_oc_w <- M_coverage_ipw_oc_we <- NULL

#seq_th <- seq(0.1,0.9,by=0.1)
seq_th <- seq(0.1,0.5,by=0.1) #use this otherwise too few controls
#seq_th <- seq(0.1,0.2,by=0.1)



coverage <- function(point,se,true){
  
  k <- qnorm(0.975) 
  
  L <- point - k*se
  U <- point + k*se
  
  return(ifelse((L < true) & (true < U) ,1,0))
  
}


itera <- 10000
k <- 0
for(th in seq_th){
  print(paste("########## % of concurrent,", 100*(1-th)))
  thresholdA <- th
  k <- k+1 
    for(i in 1:itera){
          if( (i == 1) | (i %% 500)==0){print(paste("Iteration,", i, "out of", itera)) }
          data <- data_genera(n,tr_effectA,beta_W_Y,beta_W_E,beta_E_Y,
                              beta_W_A,beta_Vk_E,beta_E_Vk,
                              heteroTE,maxE,thresholdA,isDeterministic)
          
          data$Z <- (data$W*data$E/5 + 0.5)^3
          
          data_Vk <- data[data$V_trA==1,]
          
          gcomp_oc_w  <-  gcomp_OC(formula_m0_W,formula_m1_W,data_Vk) #IMPORTANT data_VK
          gcomp_oc_we <-  gcomp_OC(formula_m0_WE,formula_m1_WE,data_Vk) #IMPORTANT data_VK
          gcomp_all_w  <- gcomp_all_W(formula_m0_WE,formula_m1_WE,data)
          gcomp_all_we <- gcomp_all_WE(formula_m0_WE,formula_m1_WE,data)
          
          ipw_oc_w  <-  IPW_OC(formula_me_W,data) 
          ipw_oc_we <-  IPW_OC(formula_me_WE,data) 
          
          aipw_oc_w  <-  aipw_OC(formula_m0_W,formula_m1_W,formula_me,data_Vk) #IMPORTANT data_VK
          aipw_oc_we <-  aipw_OC(formula_m0_WE,formula_m1_WE,formula_me,data_Vk) #IMPORTANT data_VK
          
          aipw_all_w  <- aipw_all(formula_m0_WE,formula_m1_WE,formula_me,data,nonpara=T)
          
          gcomp_point_oc_w[i]   <- gcomp_oc_w$mu
          gcomp_point_oc_we[i]  <- gcomp_oc_we$mu
          gcomp_point_all_w[i]  <- gcomp_all_w$mu
          gcomp_point_all_we[i] <- gcomp_all_we$mu
          
          IPW_point_oc_w[i]   <- ipw_oc_w$mu
          IPW_point_oc_we[i]  <- ipw_oc_we$mu
          
          aipw_point_oc_w[i]   <- aipw_oc_w$mu
          aipw_point_oc_we[i]  <- aipw_oc_we$mu
          aipw_point_all[i]  <- aipw_all_w$mu
          
          se_oc_w[i]   <- gcomp_oc_w$se_if
          se_oc_we[i]  <- gcomp_oc_we$se_if
          se_all_w[i]  <- gcomp_all_w$se_if
          se_all_we[i]  <- gcomp_all_we$se_if
          
          ipw_se_oc_w[i] <- ipw_oc_w$se_if
          ipw_se_oc_we[i] <- ipw_oc_we$se_if
          
          aipw_se_oc_w[i]   <- aipw_oc_w$se_if
          aipw_se_oc_we[i]  <- aipw_oc_we$se_if
          aipw_se_all[i]  <- aipw_all_w$se_if
          true[i] <- mean(data_Vk$Y1) - mean(data_Vk$Y0)
          
          mse_gcomp_point_oc_w[i] <- (gcomp_point_oc_w[i] - true[i])^2
          mse_gcomp_point_oc_we[i]  <- (gcomp_point_oc_we[i] - true[i])^2
          mse_gcomp_point_all_we[i]  <- (gcomp_point_all_we[i] - true[i])^2
          
          mse_aipw_point_oc_w[i]  <- (aipw_point_oc_w[i] - true[i])^2
          mse_aipw_point_oc_we[i]  <- (aipw_point_oc_we[i] - true[i])^2
          mse_aipw_point_all_we[i]  <- (aipw_point_all[i] - true[i])^2
          
          mse_ipw_point_oc_w[i]  <- (IPW_point_oc_w[i] - true[i])^2
          mse_ipw_point_oc_we[i]  <- (IPW_point_oc_we[i] - true[i])^2
          
          coverage_gcompOR_w[i]  <- coverage(gcomp_point_oc_w[i],
                                             se_oc_w[i],true[i])
          
          coverage_gcompOR_we[i] <- coverage(gcomp_point_oc_we[i],
                                             se_oc_we[i],true[i])
          
          coverage_gcompall_w[i] <- coverage(gcomp_point_all_w[i],
                                             se_all_w[i],true[i])
          
          coverage_gcompall_we[i] <- coverage(gcomp_point_all_we[i],
                                              se_all_we[i],true[i])
          
          
          coverage_ipw_point_oc_w[i] <- coverage(IPW_point_oc_w[i],
                                                 ipw_se_oc_w[i],true[i])
          
          coverage_ipw_point_oc_we[i] <- coverage(IPW_point_oc_we[i],
                                                  ipw_se_oc_we[i],true[i])
          
          
          coverage_aipw_point_oc_w[i] <- coverage(aipw_point_oc_w[i],
                                                  aipw_se_oc_w[i],true[i])
          
          coverage_aipw_point_oc_we[i] <- coverage(aipw_point_oc_we[i],
                                                   aipw_se_oc_we[i],true[i])
          
          coverage_aipw_point_all_w[i] <- coverage(aipw_point_all[i],
                                                   aipw_se_all[i],true[i])
      
          
      
    }#end for (i)
  
    
    #bias
    bias_gcomp_point_oc_w[k] <- (mean(gcomp_point_oc_w,na.rm = TRUE) - mean(true))
    bias_gcomp_point_oc_we[k] <- (mean(gcomp_point_oc_we,na.rm = TRUE) - mean(true))
    bias_gcomp_point_all_we[k] <- (mean(gcomp_point_all_we,na.rm = TRUE) - mean(true))
    
    bias_ipw_point_oc_w[k] <- (mean(IPW_point_oc_w,na.rm = TRUE) - mean(true))
    bias_ipw_point_oc_we[k] <- (mean(IPW_point_oc_we,na.rm = TRUE) - mean(true))
    
    bias_aipw_point_oc_w[k] <- (mean(aipw_point_oc_w,na.rm = TRUE) - mean(true))
    bias_aipw_point_oc_we[k] <- (mean(aipw_point_oc_we,na.rm = TRUE) - mean(true))
    bias_aipw_point_all_we[k] <- (mean(aipw_point_all,na.rm = TRUE) - mean(true))
    
    #Standard error
    se_gcomp_point_oc_w_sv[k] <- sd(gcomp_point_oc_w,na.rm = TRUE)
    se_gcomp_point_oc_we_sv[k] <- sd(gcomp_point_oc_we,na.rm = TRUE)
    se_gcomp_point_all_we_sv[k] <- sd(gcomp_point_all_we,na.rm = TRUE)
    
    se_gcomp_point_oc_w_if[k] <- mean(se_oc_w,na.rm = TRUE)
    se_gcomp_point_oc_we_if[k] <- mean(se_oc_we,na.rm = TRUE)
    se_gcomp_point_all_we_if[k] <- mean(se_all_we,na.rm = TRUE)
    
    
    se_ipw_point_oc_w_sv[k] <- sd(IPW_point_oc_w,na.rm = TRUE)
    se_ipw_point_oc_we_sv[k] <- sd(IPW_point_oc_we,na.rm = TRUE)
    se_ipw_point_oc_w_if[k] <- mean(ipw_se_oc_w,na.rm = TRUE)
    se_ipw_point_oc_we_if[k] <- mean(ipw_se_oc_we,na.rm = TRUE)
    
    
    se_aipw_point_oc_w_sv[k] <- sd(aipw_point_oc_w,na.rm = TRUE)
    se_aipw_point_oc_we_sv[k] <- sd(aipw_point_oc_we,na.rm = TRUE)
    se_aipw_point_all_we_sv[k] <- sd(aipw_point_all,na.rm = TRUE)
    
    se_aipw_point_oc_w_if[k] <- mean(aipw_se_oc_w,na.rm = TRUE)
    se_aipw_point_oc_we_if[k] <- mean(aipw_se_oc_we,na.rm = TRUE)
    se_aipw_point_all_we_if[k] <- mean(aipw_se_all,na.rm = TRUE)
    

    
    #efficiency gains:
    # sd(gcomp_point_oc_w)/sd(gcomp_point_oc_we)
    ratio_gcomp_oc_w_vs_all[k] <- sd(gcomp_point_oc_w,na.rm = TRUE)/sd(gcomp_point_all_w,na.rm = TRUE)
    ratio_gcomp_oc_we_vs_all[k] <- sd(gcomp_point_oc_we,na.rm = TRUE)/sd(gcomp_point_all_we,na.rm = TRUE)
    
    # sd(aipw_point_oc_w)/sd(aipw_point_oc_we)
    ratio_aipw_oc_w_vs_all[k] <- sd(aipw_point_oc_w,na.rm = TRUE)/sd(aipw_point_all,na.rm = TRUE)
    ratio_aipw_oc_we_vs_all[k] <- sd(aipw_point_oc_we,na.rm = TRUE)/sd(aipw_point_all,na.rm = TRUE)
    
    M_mse_gcomp_point_oc_w[k] <- mean(mse_gcomp_point_oc_w,na.rm = TRUE)
    M_mse_gcomp_point_oc_we[k] <- mean(mse_gcomp_point_oc_we,na.rm = TRUE)
    M_mse_gcomp_point_all_we[k] <- mean(mse_gcomp_point_all_we,na.rm = TRUE)
    
    M_mse_aipw_point_oc_w[k] <- mean(mse_aipw_point_oc_w,na.rm = TRUE)
    M_mse_aipw_point_oc_we[k] <- mean(mse_aipw_point_oc_we,na.rm = TRUE)
    M_mse_aipw_point_all_we[k] <- mean(mse_aipw_point_all_we,na.rm = TRUE)
    
    
    M_mse_ipw_point_oc_w[k] <- mean(mse_ipw_point_oc_w,na.rm = TRUE)
    M_mse_ipw_point_oc_we[k] <- mean(mse_ipw_point_oc_we,na.rm = TRUE)
    
  
    M_coverage_gcompOR_w[k] <- mean(coverage_gcompOR_w,na.rm = TRUE)
    M_coverage_gcompOR_we[k] <- mean(coverage_gcompOR_we,na.rm = TRUE)
    M_coverage_gcomp_all_w[k] <- mean(coverage_gcompall_w,na.rm = TRUE)
    M_coverage_gcomall_we[k] <- mean(coverage_gcompall_we,na.rm = TRUE)

    M_coverage_ipw_oc_w[k] <- mean(coverage_ipw_point_oc_w,na.rm = TRUE)
    M_coverage_ipw_oc_we[k] <- mean(coverage_ipw_point_oc_we,na.rm = TRUE)
        
    M_coverage_aipw_oc_w[k] <- mean(coverage_aipw_point_oc_w,na.rm = TRUE)
    M_coverage_aipw_oc_we[k] <- mean(coverage_aipw_point_oc_we,na.rm = TRUE)
    M_coverage_aipw_all[k] <- mean(coverage_aipw_point_all_w,na.rm = TRUE)
    
    

}#end threshold



bias <- c(bias_gcomp_point_oc_w,
          bias_gcomp_point_oc_we,
          bias_gcomp_point_all_we,
          bias_ipw_point_oc_w,
          bias_ipw_point_oc_we,
          bias_aipw_point_oc_w,
          bias_aipw_point_oc_we,
          bias_aipw_point_all_we)


se_sv <- c(se_gcomp_point_oc_w_sv,
           se_gcomp_point_oc_we_sv,
           se_gcomp_point_all_we_sv,
           se_ipw_point_oc_w_sv,
           se_ipw_point_oc_we_sv,
           se_aipw_point_oc_w_sv,
           se_aipw_point_oc_we_sv,
           se_aipw_point_all_we_sv)

se_if <- c(se_gcomp_point_oc_w_if,
           se_gcomp_point_oc_we_if,
           se_gcomp_point_all_we_if,
           se_ipw_point_oc_w_if,
           se_ipw_point_oc_we_if,
           se_aipw_point_oc_w_if,
           se_aipw_point_oc_we_if,
           se_aipw_point_all_we_if)

mse <- c(M_mse_gcomp_point_oc_w,
         M_mse_gcomp_point_oc_we,
         M_mse_gcomp_point_all_we,
         M_mse_ipw_point_oc_w,
         M_mse_ipw_point_oc_we,
         M_mse_aipw_point_oc_w,
         M_mse_aipw_point_oc_we,
         M_mse_aipw_point_all_we)

coverage <- c(M_coverage_gcompOR_w,
              M_coverage_gcompOR_we,
              M_coverage_gcomp_all_w,
              M_coverage_ipw_oc_w,
              M_coverage_ipw_oc_we, 
              M_coverage_aipw_oc_w,
              M_coverage_aipw_oc_we, 
              M_coverage_aipw_all)


method <- c(rep("Gcomp-OC-W",length(seq_th)),
            rep("Gcomp-OC-WE",length(seq_th)),
            rep("Gcomp-All",length(seq_th)),
            rep("IPW-OC-W",length(seq_th)),
            rep("IPW-OC-WE",length(seq_th)),
            rep("AIPW-OC-W",length(seq_th)),
            rep("AIPW-OC-WE",length(seq_th)),
            rep("AIPW-All",length(seq_th)))

threshold <- rep(seq_th,8)


data_plot <- data.frame(bias,se_sv,se_if,mse,coverage,method,threshold)

data_plot

#save.image("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24_n100.RData")
save.image("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24_n100.RData")
#save.image("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24_n100.RData")
#save.image("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24_n100.RData")






