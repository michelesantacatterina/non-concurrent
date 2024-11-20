rm(list=ls())

library(geex)
library(sandwich)
library(ggplot2)
library(dplyr)

set.seed(1234)

inv_logit <- function(x){
  1/(1 + exp(-x))
}


softmax <- function(x) {
  exp(x) / sum(exp(x))
}

data_genera <- function(n,
                        tr_effectA, 
                        tr_effectB,
                        beta_W_Y=0.25,
                        beta_W_E,
                        beta_E_Y, 
                        beta_W_A,
                        beta_W_B,
                        beta_W_AB,
                        beta_Vk_E = 0,
                        beta_E_Vk,
                        heteroTE,
                        maxE,
                        thresholdA,
                        thresholdB,
                        isDeterministic){
  
  #E          <- sort(runif(n,0,maxE))
  E          <- sort(rnorm(n,0,1))
  W          <- -mean(beta_W_E*E) + beta_W_E*E + rnorm(n,0,1)
  #thresholdA <- runif(1,0.2,0.8)
  thresholdA <- thresholdA*maxE #more or less (1 - proportion of concurrent)
  
  #V is a deterministic function of E
  if(isTRUE(isDeterministic)){
    V_trA      <- ifelse(E>thresholdA,0,1) #this mean that ~100*thresholdA % nonconcurrent vs 100*(1-thresholdA)% concurrent
    V_trB      <- ifelse(E<thresholdB,0,1)
  }
  
  #V is NOT a deterministic function of E
  # if(!isTRUE(isDeterministic)){
  #   pr_Vk <- (1+exp( -(-thresholdA -mean(E)*beta_E_Vk + E*beta_E_Vk ) ))^-1
  #   V_trA <- rbinom(n,1,pr_Vk)
  # }
  
  
  # The setting is the following:
  # Control       XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  # Treatment A   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  # Treatment B                               XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  #               V_trtA == 1 & V_trtB == 0   V_trtA &    V_trtA == 0 & V_trtB == 1
  #                                           V_trtB == 1
  
  
  Trt          <- rep(0,n)
  
  #Control and A available (only)
  WVtrA <- W[which(V_trA==1 & V_trB == 0)]
  prA        <- (1+exp( -(-mean(WVtrA)*beta_W_A + beta_W_A*WVtrA) ))^-1 #model (a) without enrichment
  
  #Control and B available (only)
  WVtrB <- W[which(V_trA==0 & V_trB == 1)]
  prB        <- (1+exp( -(-mean(WVtrB)*beta_W_B + beta_W_B*WVtrB) ))^-1
  
  #Control A and B available 
  WVtrAB <- W[which(V_trA==1 & V_trB == 1)]
  
  # Calculate linear predictors
  z1 <- -mean(WVtrAB)*beta_W_AB + beta_W_AB*WVtrAB
  z2 <- -mean(WVtrAB)*beta_W_AB + beta_W_AB*WVtrAB
  
  # Calculate probabilities
  prAB <- cbind(exp(0), exp(z1), exp(z2))
  prAB <- prAB / rowSums(prAB)


  #Control and A available (only);
  #receive either control or treatment A
  Trt[which(V_trA==1 & V_trB == 0)] <- rbinom(sum(V_trA==1 & V_trB == 0),1,prob=prA)
  
  #Control and B available (only);
  #receive either control or treatment B
  Trt[which(V_trA==0 & V_trB == 1)] <- rbinom(sum(V_trA==0 & V_trB == 1),1,prob=prB) * 2 # I multiply by 2 because I want trt B to be 2
  
  #Among those for which treatment A and B are available;
  #sample from multinomial with probability prAB
  Trt[which(V_trA==1 & V_trB == 1)] <- 
    as.numeric( apply(prAB, 1, function(x) sample(0:2, size = 1, prob = x)) )
  
  
  Y0 <- beta_W_Y*W + beta_E_Y*E + beta_Vk_E*V_trA*E + rnorm(n)
  # Y0 <- NULL
  # na0 <- sum(1-V_trA)
  # na1 <- sum(V_trA)
  # Y0[which(V_trA==0)] <- beta_W_Y*W[which(V_trA==0)] + beta_E_Y*E[which(V_trA==0)] + 
  #   beta_Vk_E*V_trA[which(V_trA==0)]*E[which(V_trA==0)] + rnorm(na0)
  # 
  # Y0[which(V_trA==1)] <- beta_W_Y*W[which(V_trA==1)] + beta_E_Y*E[which(V_trA==1)] + 
  #   beta_Vk_E*V_trA[which(V_trA==1)]*E[which(V_trA==1)] + rnorm(na1)
  
  
  # if(isTRUE(heteroTE)){
  #   Y1 <- NULL
  #   Y1[which(V_trA==0)]   <- Y0[which(V_trA==0)] - tr_effectA #harmful among the non concurrent
  #   Y1[which(V_trA==1)]   <- Y0[which(V_trA==1)] + tr_effectA #beneficial among the concurrent
  # }else{
    Y1  <- Y0 + tr_effectA # homogeneous treatment effect
    Y2  <- Y0 + tr_effectB # homogeneous treatment effect
  # }
  
  Y               <- Y0
  Y[which(Trt==1)]  <- Y1[which(Trt==1)]
  Y[which(Trt==2)]  <- Y2[which(Trt==2)]
  
  
  data <- data.frame(Y,Y2,Y1,Y0,W,E,Trt,V_trA,V_trB,thresholdA,thresholdB)
  
  # summary(glm(A~W,data=data,family = "binomial"))
  # data_Vk <- data[data$V_trA==1,]
  # summary(glm(A~W,data=data_Vk,family = "binomial"))
  # 
  # summary(glm(Y~W+E,data=data[data$V_trA==1,]))
  
  return(data)
  
}


# Function to compute the target population

target_pop_fun <- function(pop_index,data){
  
  target_pop <- rep(0, nrow(data))
  
  if(pop_index == "A1"){ target_pop[data$V_trA == 1] <- 1 } 
  
  if(pop_index == "B1"){ target_pop[data$V_trB == 1] <- 1 } 
  
  if(pop_index == "A1B0"){ target_pop[data$V_trA==1 & data$V_trB == 0] <- 1 } 
  
  if(pop_index == "A1B1"){ target_pop[data$V_trA==1 & data$V_trB == 1] <- 1 }
  
  if(pop_index == "A0B1"){ target_pop[data$V_trA==0 & data$V_trB == 1] <- 1 }
  
  return(target_pop)
}


# Function to compute the true values of the targeted parameters

true_par_fun <- function(pop_index,trt_contrast,data){
  
  target_pop <- target_pop_fun(pop_index,data)
  
  if(trt_contrast == 1){
    true_par <- sum( target_pop*(data$Y1 - data$Y0) ) / sum(target_pop)
  }
  
  if(trt_contrast == 2){
    true_par <- sum( target_pop*(data$Y2 - data$Y0) ) / sum(target_pop)
  }
  
  return(true_par)
}





################################################################################
#Gcomp - EIF

#--------------------Concurrent only 
gcomp_OC <- function(formula_m0,formula_m1,trt_contrast,data,pop_index){
  
  #----get the target population
  target_pop <- target_pop_fun(pop_index,data)
  
  #----------controls
  
  #regress among controls only concurrent 
  m0_only_conc <- glm(formula_m0,data=data[data$Trt==0 & target_pop==1,])
  
  Y <- data$Y
  X <- model.matrix(formula_m0, data = data,drop = FALSE)
  
  xbeta <- (X) %*% coef(m0_only_conc) # same as predict(m0_only_conc, newdata = data)
  #e0 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
  e0 <- sum( target_pop*xbeta)/sum(target_pop) 
  
  h0 <- X * target_pop * (data$Trt==0) * as.numeric(t(Y -  xbeta) )
  h1 <- target_pop*(xbeta - e0)
  
  #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
  g11 <- ( t(X*target_pop*(data$Trt==0)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  g21 <- apply(-X*target_pop,2,sum)/nrow(data)
  
  g22 <- sum(target_pop)
  g22inv <- 1/g22
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
  # #same as
  # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
  var_beta_hat0 <- phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
  
  #phi_mu0 <- n*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
  phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
  #se <- sqrt( (t(phi_mu0) %*% phi_mu0 )/(nrow(data)-1) )/(nrow(data))
  se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
  
  
  #----------treated 
  
  
  #regress among controls only concurrent 
  m1_only_conc <- glm(formula_m1,data=data[data$Trt==trt_contrast & target_pop==1,])
  
  Y <- data$Y
  X <- model.matrix(formula_m1, data = data,drop = FALSE)
  
  xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
  e1 <- sum(target_pop*xbeta)/sum(target_pop) 
  
  h0 <- X * target_pop * (data$Trt==trt_contrast) * as.numeric(t(Y -  xbeta) )
  h1 <- target_pop*(xbeta - e1)
  
  #g11 <- ( t(X0) %*% diag(data$A) %*% X0) / nrow(data)
  g11 <- ( t(X*target_pop*(data$Trt==trt_contrast)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  g21 <- apply(-X*target_pop,2,sum)/nrow(data)
  
  g22 <- sum(target_pop)
  g22inv <- 1/g22
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m1_only_conc))
  # # same as
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
gcomp_all_WE <- function(formula_m0,formula_m1,trt_contrast,data,pop_index){
  
  #----get the target population
  target_pop <- target_pop_fun(pop_index,data)
  
  #----------controls
  
  #regress among controls all controls 
  m0_all_cont <- glm(formula_m0,data=data[data$Trt==0,])
  
  Y <- data$Y
  X <- model.matrix(formula_m0, data = data,drop = FALSE)
  
  xbeta <- (X) %*% coef(m0_all_cont) # same as predict(m0_only_cont, newdata = data)
  e0 <- sum(target_pop*xbeta)/sum(target_pop) 
  
  h0 <- X * (data$Trt==0) * as.numeric(t(Y -  xbeta) ) 
  h1 <- target_pop*(xbeta - e0)
  
  #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
  g11 <- ( t(X*(data$Trt==0)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  g21 <- apply(-X*target_pop,2,sum)/nrow(data)
  
  g22 <- sum(target_pop)
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
  m1_only_conc <- glm(formula_m1,data=data[data$Trt==trt_contrast & target_pop==1,])
  
  Y <- data$Y
  X <- model.matrix(formula_m1, data = data,drop = FALSE)
  
  xbeta <- (X) %*% coef(m1_only_conc) # same as predict(m1_only_cont, newdata = data)
  e1 <- sum(target_pop*xbeta)/sum(target_pop) 
  
  h0 <- X * target_pop * (data$Trt==trt_contrast) * as.numeric(t(Y -  xbeta) )
  h1 <- target_pop*(xbeta - e1)
  
  #g11 <- ( t(X0) %*% diag(data$A) %*% X0) / nrow(data)
  g11 <- ( t(X*target_pop*(data$Trt==trt_contrast)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  g21 <- apply(-X*target_pop,2,sum)/nrow(data)
  
  g22 <- sum(target_pop)
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


##******************************************************************************
## Test simu

      # ### Itera
      # 
      # point_oc_we <- point_all_we <- true <- NULL
      # se_oc_we <- se_all_we <- NULL
      # itera <- 1000
      # for(i in 1:itera){
      #   if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
      #   data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
      #                       beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
      #                       heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
      # 
      #   # pop_index <- "A1"
      #   # trt_contrast <- 1
      #   
      #   pop_index <- "B1"
      #   trt_contrast <- 2
      #   
      #   #----get the target population
      #   target_pop <- target_pop_fun(pop_index,data)
      #   data_Vk <- data[target_pop==1,]
      # 
      #   oc_we <-  gcomp_OC(formula_m0,formula_m1,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
      # 
      #   all_we <- gcomp_all_WE(formula_m0_WE,formula_m1_WE,trt_contrast,data,pop_index)
      # 
      #   point_oc_we[i]  <- oc_we$mu
      #   point_all_we[i] <- all_we$mu
      # 
      #   se_oc_we[i]  <- oc_we$se_if
      #   se_all_we[i] <- all_we$se_if
      #   true[i] <- true_par_fun(pop_index,trt_contrast,data)
      # 
      # }
      # 
      # mean(point_oc_we);sd(point_oc_we);mean(se_oc_we)
      # 
      # mean(point_all_we);sd(point_all_we);mean(se_all_we)
      # 
      # sd(point_oc_we)/sd(point_all_we)
      # 
      # mean(true)
      
      
      
      
      
################################################################################
#IPW - EIF
      
      #--------------------Concurrent only 
      IPW_OC <- function(formula_me,trt_contrast,data,pop_index){
        
        #----get the target population
        data <- data[(data$Trt==0 | data$Trt==trt_contrast),]
        target_pop <- target_pop_fun(pop_index,data=data)
        
        #regress among concurrent
        m0_only_conc <- glm(formula_me,data=data[target_pop==1,],family="binomial")
        X <- model.matrix(formula_me, data = data,drop = FALSE)
        
        #propensity score among concurrent
        xbeta <- (X) %*% coef(m0_only_conc) 
        exp_xbeta <- as.numeric(exp(xbeta))
        exp_m_xbeta <- exp(-xbeta)
        ps <- as.numeric(exp_xbeta/(1+exp_xbeta)) #same as predict(m0_only_conc, newdata = data, type="response")
        
        pv <- mean(target_pop)
        #----------Controls
        
        w0 <- ((data$Trt==0))/(1-ps)
        #e0 <- sum(w0*target_pop*data$Y)/sum(w0*target_pop) #mean(data_Vk$Y0)
        e0 <- mean( (target_pop/pv) *( ( w0 )*(data$Y) ) )
        
        h0 <- X * target_pop * ( data$Trt==0  - (1-as.numeric(t(ps))) )
        
        g11 <- ( t(X * target_pop * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
        g11inv <- solve(g11)  
        
        phi_beta_hat <- g11inv %*% t(h0)
        # B <- t(h0) %*% h0
        # A <- g11inv
        # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
        # # same as
        # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
        
        
        h1 <- target_pop*(w0*data$Y - e0)
        
        a0 <- data$Trt==0
        #g21 <- as.numeric((- t( a0 * data$V_trA * data$Y *exp_xbeta ) %*% X)/nrow(data))
        g21 <- as.numeric((t( a0 * target_pop * data$Y *exp_xbeta ) %*% X)/nrow(data))
        
        #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
        
        g22 <- sum(target_pop)
        g22inv <- 1/g22
        
        
        phi_mu0 <- nrow(data)*g22inv * ( h1 - (t(phi_beta_hat)%*%g21) ) 
        se0 <- sqrt(var(phi_mu0)/nrow(data)) * sqrt((nrow(data)-1)/nrow(data))
        
        
        
        #----------Treated
        
        w1 <- (target_pop*(data$Trt==trt_contrast))/(ps)
        #e1 <- sum(w1*data$V_trA*data$Y)/sum(w1*data$V_trA) #mean(data_Vk$Y1)
        e1 <- mean( (target_pop/pv) *( ( w1 )*(data$Y) ) )
        
        h0 <- X * target_pop * ( (data$Trt==trt_contrast) - as.numeric(t(ps)) )
        
        g11 <- ( t(X * target_pop * ( exp_xbeta/((exp_xbeta+1)^2) ) ) %*% X ) / nrow(data)
        g11inv <- solve(g11)  
        
        phi_beta_hat <- g11inv %*% t(h0)
        # B <- t(h0) %*% h0
        # A <- g11inv
        # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_only_conc))
        # same as 
        # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
        
        
        h1 <- target_pop*(w1*data$Y - e1)
        
        a1 <- data$Trt==trt_contrast
        g21 <- as.numeric((t( a1 * target_pop * data$Y *exp_m_xbeta ) %*% X)/nrow(data))
        
        #apply(-X*a0*data$V_trA*data$Y*exp_xbeta,2,sum)/nrow(data)
        
        g22 <- sum(target_pop)
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
      
      
# ##******************************************************************************
# ## Test simu
# 
#     ### Itera
# 
#     point_oc_we <- true <- NULL
#     se_oc_we <-  NULL
#     itera <- 1000
#     for(i in 1:itera){
#       if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#       data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
#                           beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
#                           heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
# 
#       # pop_index <- "A1"
#       # trt_contrast <- 1
# 
#       pop_index <- "B1"
#       trt_contrast <- 2
# 
#       #----get the target population
#       target_pop <- target_pop_fun(pop_index,data)
#       data_Vk <- data[target_pop==1,]
# 
#       oc_we <-  IPW_OC(formula_me,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
# 
#       point_oc_we[i]  <- oc_we$mu
# 
#       se_oc_we[i]  <- oc_we$se_if
#       true[i] <- true_par_fun(pop_index,trt_contrast,data)
# 
#     }
# 
#     mean(point_oc_we);sd(point_oc_we);mean(se_oc_we)
# 
#     mean(true)
      


################################################################################


    #-----------AIPW conc
    aipw_OC <- function(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index){
      
      #----get the target population
      data <- data[(data$Trt==0 | data$Trt==trt_contrast),]
      target_pop <- target_pop_fun(pop_index,data=data)
      
      data_Vk <- data[target_pop==1,]
      n <- nrow(data)
      
      #regress among only concurrent controls
      m0_only_concC <- glm(formula_m0,data=data_Vk[data_Vk$Trt==0,])
      
      #predict among only concurrent controls
      p0_only_concC_m0_only_concC <- predict(m0_only_concC, newdata = data_Vk)
      
      #regress among only concurrent treated
      mk_only_concC <- glm(formula_m1,data=data_Vk[data_Vk$Trt==trt_contrast,])
      
      #predict among only concurrent treated
      pk_only_concC_mk_only_concC <- predict(mk_only_concC, newdata = data_Vk)
      
      #obtain propensities
      fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
      
      #obtain propensities in the Vk population
      eps <- predict(fitTr, type="response", newdata=data_Vk) 
      
      #obtain probability of Vk
      pv <- mean(target_pop)
      
      #indicator of Vk = 1
      IVk <- as.numeric(target_pop==1)
      
      
      #obtain phi1
      temp1                   <- ((data_Vk$Trt==trt_contrast)/eps)*(data_Vk$Y-pk_only_concC_mk_only_concC) + pk_only_concC_mk_only_concC
      temp1_2                 <- rep(0,n)
      temp1_2[target_pop==1]  <-  temp1
      phi1                    <-  (IVk/pv) * temp1_2
      se1 <- (sqrt(var(phi1)/n)) * sqrt((n-1)/n)
      e1 <- mean(phi1)
      
      #obtain phi0
      temp0                   <- ( (data_Vk$Trt==0)/(1-eps) )*(data_Vk$Y-p0_only_concC_m0_only_concC) + p0_only_concC_m0_only_concC
      temp0_2                 <- rep(0,n)
      temp0_2[target_pop==1]  <-  temp0
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
    aipw_all <- function(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index,nonpara){
      
      
      #----get the target population
      data <- data[(data$Trt==0 | data$Trt==trt_contrast),]
      target_pop <- target_pop_fun(pop_index,data=data)
      
      data_Vk <- data[target_pop==1,]
      n <- nrow(data)
      
      #regress among all controls
      m0_allC <- glm(formula_m0,data=data[data$Trt==0,])
      Y <- data$Y
      X <- model.matrix(formula_m0, data = data,drop = FALSE)
      # 
      #predict among only concurrent controls
      #p0_only_concC_m0_allC <- predict(m0_allC, newdata = data_Vk)
      p0_only_concC_m0_allC <- (X) %*% coef(m0_allC)
      
      #regress among only concurrent treated
      mk_only_concC <- glm(formula_m1,data=data_Vk[data_Vk$Trt==trt_contrast,])
      
      #predict among only concurrent treated
      pk_only_concC_mk_only_concC <- predict(mk_only_concC, newdata = data_Vk)
      
      #obtain propensities
      fitTr <- glm(formula_me, data = data_Vk, family = binomial(link = "logit")) #V_trA is NA
      
      #obtain propensities in the Vk population
      eps <- predict(fitTr, type="response", newdata=data_Vk) 
      
      #obtain probability of Vk
      pv <- mean(target_pop)
      
      #indicator of Vk = 1
      IVk <- as.numeric(target_pop==1)
      
      #obtain phi1
      temp1                   <-  ((data_Vk$Trt==trt_contrast)/eps)*(data_Vk$Y-pk_only_concC_mk_only_concC) + pk_only_concC_mk_only_concC
      temp1_2                 <-  rep(0,n)
      temp1_2[target_pop==1]  <-  temp1
      phi1                    <-  (IVk/pv) * temp1_2
      se1 <- (sqrt(var(phi1)/n)) * sqrt((n-1)/n)
      e1 <- mean(phi1)
      
      #obtain phi0
      # temp0                   <-  ( (1-data_Vk$A)/(1-eps) )*(data_Vk$Y-p0_only_concC_m0_allC) + p0_only_concC_m0_allC
      eps0 <- rep(0,n)
      eps0[target_pop==1] <- eps
      #temp0                   <-  ( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
      #temp0                   <-  (IVk/pv) *( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + p0_only_concC_m0_allC )
      # temp0_2                 <-  rep(0,n)
      # temp0_2[target_pop==1]  <-  temp0
      
      if(isTRUE(nonpara)){
        temp0 <-  (IVk/pv) *( ((data$Trt==0)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + p0_only_concC_m0_allC )
      }else{
        temp0 <- ( ((data$Trt==0)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
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



# ##******************************************************************************
# ## Test simu
# 
#     ### Itera
# 
#     point_oc_we <- point_all_we <- true <- NULL
#     se_oc_we <- se_all_we <- NULL
#     itera <- 1000
#     for(i in 1:itera){
#       if( (i == 1) | (i %% 250)==0){print(paste("Iteration,", i, "out of", itera)) }
#       data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
#                           beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
#                           heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
# 
#       # pop_index <- "A1"
#       # trt_contrast <- 1
# 
#       pop_index <- "A1"
#       trt_contrast <- 1
# 
#       #----get the target population
#       target_pop <- target_pop_fun(pop_index,data)
#       data_Vk <- data[target_pop==1,]
# 
#       oc_we <-  aipw_OC(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index) #IMPORTANT data_VK
# 
#       all_we <- aipw_all(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index,nonpara=T)
# 
#       point_oc_we[i]  <- oc_we$mu
#       point_all_we[i] <- all_we$mu
# 
#       se_oc_we[i]  <- oc_we$se_if
#       se_all_we[i] <- all_we$se_if
#       true[i] <- true_par_fun(pop_index,trt_contrast,data)
# 
#     }
# 
#     mean(point_oc_we);sd(point_oc_we);mean(se_oc_we)
# 
#     mean(point_all_we);sd(point_all_we);mean(se_all_we)
# 
#     sd(point_oc_we)/sd(point_all_we)
# 
#     mean(true)
      
      
      
            
# Generate data
      tr_effectA      <- 0.5    #medium treatment effect 
      tr_effectB      <- 1      #large treatment effect 
      beta_E_Y        <- 0.5    #effect size E->Y
      
      beta_Vk_E       <- 0      #no interaction Vk/E 
      maxE            <- 1      #max entry time (standardized between 0 and maxE)
      heteroTE        <- FALSE  #homogeneous treatment effect
      thresholdA      <- 0.6   #40% concurrent ::: 100*(1-thresholdA)% concurrent
      thresholdB      <- 0.4   #smaller than thresholdA to allow for overlap between V_A and V_B
      beta_W_A        <- 0.8 + 0*(thresholdA)*maxE   #effect size W->A 
      beta_W_B        <- 0.8 + 0*(thresholdA)*maxE   #effect size W->A 
      beta_W_AB       <- 0.8 + 0*(thresholdA)*maxE   #effect size W->A 
      beta_W_Y        <- 0.8 + 0*(thresholdA)*maxE  #effect size W->Y
      beta_W_E        <- 0.8 + 0*(thresholdA)*maxE
      
      isDeterministic <- TRUE
      
      n <- 1000
      itera <- 1000
      
      
      formula_m0 <- Y ~ W + E
      formula_m1 <- Y ~ W + E
      formula_me <- factor(Trt) ~ W + E
      
      
      
      data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
                          beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
                          heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
      
      
      
      table(data$Trt[data$V_trA==1 & data$V_trB ==0])
      sum(data$V_trA==1 & data$V_trB ==0)
      
      table(data$Trt[data$V_trA==1 & data$V_trB ==1])
      sum(data$V_trA==1 & data$V_trB ==1)
      
      table(data$Trt[data$V_trA==0 & data$V_trB ==1])
      sum(data$V_trA==0 & data$V_trB ==1)
      
      
      
      
      
      
      
      
      
      
################################################################################
################################################################################
################################################################################


gcomp_point_oc_we <- gcomp_point_all_we <- true <- NULL
aipw_point_oc_we <- aipw_point_all <- NULL
gcomp_se_oc_we <- gcomp_se_all_we <- NULL
aipw_se_oc_we <- aipw_se_all <- NULL

IPW_point_oc_we <- ipw_se_oc_w <- ipw_se_oc_we <- NULL

se_oc_we <- se_all_w <- se_all_we <- NULL


mse_gcomp_point_oc_we <- mse_gcomp_point_all_we <- 
mse_aipw_point_oc_we <- mse_aipw_point_all_we <- 
mse_ipw_point_oc_we <- NULL

bias_gcomp_point_oc_we <- bias_gcomp_point_all_we <- 
bias_aipw_point_oc_we <- bias_aipw_point_all_we <-     
bias_ipw_point_oc_we <- NULL

se_ipw_point_oc_we_sv <- se_ipw_point_oc_we_if <- NULL

se_gcomp_point_oc_we_sv <- se_gcomp_point_all_we_sv <- 
se_gcomp_point_oc_we_if <- se_gcomp_point_all_we_if <- NULL

se_aipw_point_oc_we_sv <- se_aipw_point_all_we_sv <- 
se_aipw_point_oc_we_if <- se_aipw_point_all_we_if <- NULL


M_mse_gcomp_point_oc_we <- M_mse_gcomp_point_all_we <-
M_mse_aipw_point_oc_we <- M_mse_aipw_point_all_we <- 
M_mse_ipw_point_oc_we <- NULL



coverage_gcompOR_we <- coverage_gcompall_we <- coverage_aipw_point_oc_we <-
  coverage_aipw_point_all_we <- coverage_ipw_point_oc_we <- NULL


M_coverage_gcompOR_we <- M_coverage_gcompall_we <- M_coverage_aipw_oc_we <-
  M_coverage_aipw_all <- M_coverage_ipw_oc_we <- NULL

seq_th <- seq(0.1,0.9,by=0.1)
#seq_th <- seq(0.1,0.2,by=0.1)



coverage_fun <- function(point,se,true){
  
  k <- qnorm(0.975) 
  
  L <- point - k*se
  U <- point + k*se
  
  return(ifelse((L < true) & (true < U) ,1,0))
  
}


################################################################################

n <- 1000
itera <- 10000
k <- 0

################################################################################
# E[Y(1) - Y(0) | V_A = 1]
# E[Y(1) - Y(0) | V_A = 1, V_B = 0]
# E[Y(1) - Y(0) | V_A = 1, V_B = 1]
################################################################################

trt_contrast <- 1
seq_1 <- c("A1","A1B0","A1B1")

for(pop_index in seq_1 ){
  print(paste("########## Population considered,", pop_index ))
  k <- k+1 
  for(i in 1:itera){
    if( (i == 1) | (i %% 500)==0){print(paste("Iteration,", i, "out of", itera)) }
    data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
                              beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
                              heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
    
    
    #----get the target population
    target_pop <- target_pop_fun(pop_index,data)
    data_Vk <- data[target_pop==1,]
    
    gcomp_oc_we <-  gcomp_OC(formula_m0,formula_m1,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
    gcomp_all_we <- gcomp_all_WE(formula_m0,formula_m1,trt_contrast,data,pop_index)
    
    ipw_oc_we <-  IPW_OC(formula_me,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
    
    aipw_oc_we <-  aipw_OC(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index) #IMPORTANT data_VK
    aipw_all_we  <- aipw_all(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index,nonpara=T)
    

    gcomp_point_oc_we[i]  <- gcomp_oc_we$mu
    gcomp_point_all_we[i] <- gcomp_all_we$mu
    
    IPW_point_oc_we[i]  <- ipw_oc_we$mu
    
    aipw_point_oc_we[i]  <- aipw_oc_we$mu
    aipw_point_all[i]  <- aipw_all_we$mu
    

    se_oc_we[i]  <- gcomp_oc_we$se_if
    se_all_we[i]  <- gcomp_all_we$se_if
    
    ipw_se_oc_we[i] <- ipw_oc_we$se_if
    
    aipw_se_oc_we[i]  <- aipw_oc_we$se_if
    aipw_se_all[i]  <- aipw_all_we$se_if
    true[i] <- true_par_fun(pop_index,trt_contrast,data)
    
    mse_gcomp_point_oc_we[i]  <- (gcomp_point_oc_we[i] - true[i])^2
    mse_gcomp_point_all_we[i]  <- (gcomp_point_all_we[i] - true[i])^2
    

    mse_aipw_point_oc_we[i]  <- (aipw_point_oc_we[i] - true[i])^2
    mse_aipw_point_all_we[i]  <- (aipw_point_all[i] - true[i])^2
    
    mse_ipw_point_oc_we[i]  <- (IPW_point_oc_we[i] - true[i])^2

    
    coverage_gcompOR_we[i] <- coverage_fun(gcomp_point_oc_we[i],
                                       se_oc_we[i],true[i])
    
    coverage_gcompall_we[i] <- coverage_fun(gcomp_point_all_we[i],
                                        se_all_we[i],true[i])
    
    coverage_ipw_point_oc_we[i] <- coverage_fun(IPW_point_oc_we[i],
                                            ipw_se_oc_we[i],true[i])
    
    
    coverage_aipw_point_oc_we[i] <- coverage_fun(aipw_point_oc_we[i],
                                             aipw_se_oc_we[i],true[i])
    
    
    coverage_aipw_point_all_we[i] <- coverage_fun(aipw_point_all[i],
                                             aipw_se_all[i],true[i])
    
    
    
  }#end for (i)
  
  
  #bias
  bias_gcomp_point_oc_we[k] <- (mean(gcomp_point_oc_we) - mean(true))
  bias_gcomp_point_all_we[k] <- (mean(gcomp_point_all_we) - mean(true))
  
  bias_ipw_point_oc_we[k] <- (mean(IPW_point_oc_we) - mean(true))
  
  bias_aipw_point_oc_we[k] <- (mean(aipw_point_oc_we) - mean(true))
  bias_aipw_point_all_we[k] <- (mean(aipw_point_all) - mean(true))
  
  #Standard error
  se_gcomp_point_oc_we_sv[k] <- sd(gcomp_point_oc_we)
  se_gcomp_point_all_we_sv[k] <- sd(gcomp_point_all_we)
  
  se_gcomp_point_oc_we_if[k] <- mean(se_oc_we)
  se_gcomp_point_all_we_if[k] <- mean(se_all_we)
  
  
  se_ipw_point_oc_we_sv[k] <- sd(IPW_point_oc_we)
  se_ipw_point_oc_we_if[k] <- mean(ipw_se_oc_we)
  
  se_aipw_point_oc_we_sv[k] <- sd(aipw_point_oc_we)
  se_aipw_point_all_we_sv[k] <- sd(aipw_point_all)
  
  se_aipw_point_oc_we_if[k] <- mean(aipw_se_oc_we)
  se_aipw_point_all_we_if[k] <- mean(aipw_se_all)
  
  
  
  M_mse_gcomp_point_oc_we[k] <- mean(mse_gcomp_point_oc_we)
  M_mse_gcomp_point_all_we[k] <- mean(mse_gcomp_point_all_we)
  
  M_mse_aipw_point_oc_we[k] <- mean(mse_aipw_point_oc_we)
  M_mse_aipw_point_all_we[k] <- mean(mse_aipw_point_all_we)
  
  M_mse_ipw_point_oc_we[k] <- mean(mse_ipw_point_oc_we)
  

  M_coverage_gcompOR_we[k] <- mean(coverage_gcompOR_we)
  M_coverage_gcompall_we[k] <- mean(coverage_gcompall_we)
  
  M_coverage_ipw_oc_we[k] <- mean(coverage_ipw_point_oc_we)
  
  M_coverage_aipw_oc_we[k] <- mean(coverage_aipw_point_oc_we)
  M_coverage_aipw_all[k] <- mean(coverage_aipw_point_all_we)
  
  
  
}#end threshold



bias <- c(
          bias_gcomp_point_oc_we,
          bias_gcomp_point_all_we,

          bias_ipw_point_oc_we,

          bias_aipw_point_oc_we,
          bias_aipw_point_all_we)


se_sv <- c(
           se_gcomp_point_oc_we_sv,
           se_gcomp_point_all_we_sv,

           se_ipw_point_oc_we_sv,

           se_aipw_point_oc_we_sv,
           se_aipw_point_all_we_sv)

se_if <- c(
           se_gcomp_point_oc_we_if,
           se_gcomp_point_all_we_if,

           se_ipw_point_oc_we_if,

           se_aipw_point_oc_we_if,
           se_aipw_point_all_we_if)

mse <- c(
         M_mse_gcomp_point_oc_we,
         M_mse_gcomp_point_all_we,

         M_mse_ipw_point_oc_we,

         M_mse_aipw_point_oc_we,
         M_mse_aipw_point_all_we)

coverage <- c(
              M_coverage_gcompOR_we,
              M_coverage_gcompall_we,

              M_coverage_ipw_oc_we, 
              
              M_coverage_aipw_oc_we, 
              M_coverage_aipw_all)


method <- c(
            rep("Gcomp-OC-WE",length(seq_1)),
            rep("Gcomp-All",length(seq_1)),
            
            rep("IPW-OC-WE",length(seq_1)),
           
            rep("AIPW-OC-WE",length(seq_1)),
            rep("AIPW-All",length(seq_1)))

pop <- rep(seq_1,5)

data_plot_1 <- data.frame(bias,se_sv,se_if,mse,coverage,method,pop)

data_plot_1




################################################################################
# E[Y(2) - Y(0) | V_B = 1]
# E[Y(2) - Y(0) | V_A = 0, V_B = 1]
# E[Y(2) - Y(0) | V_A = 1, V_B = 1]
################################################################################

gcomp_point_oc_we <- gcomp_point_all_we <- true <- NULL
aipw_point_oc_we <- aipw_point_all <- NULL
gcomp_se_oc_we <- gcomp_se_all_we <- NULL
aipw_se_oc_we <- aipw_se_all <- NULL

IPW_point_oc_we <- ipw_se_oc_w <- ipw_se_oc_we <- NULL

se_oc_we <- se_all_w <- se_all_we <- NULL


mse_gcomp_point_oc_we <- mse_gcomp_point_all_we <- 
  mse_aipw_point_oc_we <- mse_aipw_point_all_we <- 
  mse_ipw_point_oc_we <- NULL

bias_gcomp_point_oc_we <- bias_gcomp_point_all_we <- 
  bias_aipw_point_oc_we <- bias_aipw_point_all_we <-     
  bias_ipw_point_oc_we <- NULL

se_ipw_point_oc_we_sv <- se_ipw_point_oc_we_if <- NULL

se_gcomp_point_oc_we_sv <- se_gcomp_point_all_we_sv <- 
  se_gcomp_point_oc_we_if <- se_gcomp_point_all_we_if <- NULL

se_aipw_point_oc_we_sv <- se_aipw_point_all_we_sv <- 
  se_aipw_point_oc_we_if <- se_aipw_point_all_we_if <- NULL


M_mse_gcomp_point_oc_we <- M_mse_gcomp_point_all_we <-
  M_mse_aipw_point_oc_we <- M_mse_aipw_point_all_we <- 
  M_mse_ipw_point_oc_we <- NULL



coverage_gcompOR_we <- coverage_gcompall_we <- coverage_aipw_point_oc_we <-
  coverage_aipw_point_all_we <- coverage_ipw_point_oc_we <- NULL


M_coverage_gcompOR_we <- M_coverage_gcompall_we <- M_coverage_aipw_oc_we <-
  M_coverage_aipw_all <- M_coverage_ipw_oc_we <- NULL


bias <- se_sv <- se_if <- mse <- coverage <- NULL
k <- 0

trt_contrast <- 2
seq_2 <- c("B1","A0B1","A1B1")


#!!!! seq_2 VS seq_1 above!!!
for(pop_index in seq_2 ){
  print(paste("########## Population considered,", pop_index ))
  k <- k+1 
  for(i in 1:itera){
    if( (i == 1) | (i %% 500)==0){print(paste("Iteration,", i, "out of", itera)) }
    data <- data_genera(n,tr_effectA,tr_effectB,beta_W_Y,beta_W_E,beta_E_Y,
                        beta_W_A,beta_W_B,beta_W_AB,beta_Vk_E,beta_E_Vk,
                        heteroTE,maxE,thresholdA,thresholdB,isDeterministic)
    
    
    #----get the target population
    target_pop <- target_pop_fun(pop_index,data)
    data_Vk <- data[target_pop==1,]
    
    gcomp_oc_we <-  gcomp_OC(formula_m0,formula_m1,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
    gcomp_all_we <- gcomp_all_WE(formula_m0,formula_m1,trt_contrast,data,pop_index)
    
    ipw_oc_we <-  IPW_OC(formula_me,trt_contrast,data_Vk,pop_index) #IMPORTANT data_VK
    
    aipw_oc_we <-  aipw_OC(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index) #IMPORTANT data_VK
    aipw_all_we  <- aipw_all(formula_m0,formula_m1,formula_me,trt_contrast,data,pop_index,nonpara=T)
    
    
    gcomp_point_oc_we[i]  <- gcomp_oc_we$mu
    gcomp_point_all_we[i] <- gcomp_all_we$mu
    
    IPW_point_oc_we[i]  <- ipw_oc_we$mu
    
    aipw_point_oc_we[i]  <- aipw_oc_we$mu
    aipw_point_all[i]  <- aipw_all_we$mu
    
    
    se_oc_we[i]  <- gcomp_oc_we$se_if
    se_all_we[i]  <- gcomp_all_we$se_if
    
    ipw_se_oc_we[i] <- ipw_oc_we$se_if
    
    aipw_se_oc_we[i]  <- aipw_oc_we$se_if
    aipw_se_all[i]  <- aipw_all_we$se_if
    true[i] <- true_par_fun(pop_index,trt_contrast,data)
    
    mse_gcomp_point_oc_we[i]  <- (gcomp_point_oc_we[i] - true[i])^2
    mse_gcomp_point_all_we[i]  <- (gcomp_point_all_we[i] - true[i])^2
    
    
    mse_aipw_point_oc_we[i]  <- (aipw_point_oc_we[i] - true[i])^2
    mse_aipw_point_all_we[i]  <- (aipw_point_all[i] - true[i])^2
    
    mse_ipw_point_oc_we[i]  <- (IPW_point_oc_we[i] - true[i])^2
    
    
    coverage_gcompOR_we[i] <- coverage_fun(gcomp_point_oc_we[i],
                                       se_oc_we[i],true[i])
    
    coverage_gcompall_we[i] <- coverage_fun(gcomp_point_all_we[i],
                                        se_all_we[i],true[i])
    
    coverage_ipw_point_oc_we[i] <- coverage_fun(IPW_point_oc_we[i],
                                            ipw_se_oc_we[i],true[i])
    
    
    coverage_aipw_point_oc_we[i] <- coverage_fun(aipw_point_oc_we[i],
                                             aipw_se_oc_we[i],true[i])
    
    
    coverage_aipw_point_all_we[i] <- coverage_fun(aipw_point_all[i],
                                              aipw_se_all[i],true[i])
    
    
    
  }#end for (i)
  
  
  #bias
  bias_gcomp_point_oc_we[k] <- (mean(gcomp_point_oc_we) - mean(true))
  bias_gcomp_point_all_we[k] <- (mean(gcomp_point_all_we) - mean(true))
  
  bias_ipw_point_oc_we[k] <- (mean(IPW_point_oc_we) - mean(true))
  
  bias_aipw_point_oc_we[k] <- (mean(aipw_point_oc_we) - mean(true))
  bias_aipw_point_all_we[k] <- (mean(aipw_point_all) - mean(true))
  
  #Standard error
  se_gcomp_point_oc_we_sv[k] <- sd(gcomp_point_oc_we)
  se_gcomp_point_all_we_sv[k] <- sd(gcomp_point_all_we)
  
  se_gcomp_point_oc_we_if[k] <- mean(se_oc_we)
  se_gcomp_point_all_we_if[k] <- mean(se_all_we)
  
  
  se_ipw_point_oc_we_sv[k] <- sd(IPW_point_oc_we)
  se_ipw_point_oc_we_if[k] <- mean(ipw_se_oc_we)
  
  se_aipw_point_oc_we_sv[k] <- sd(aipw_point_oc_we)
  se_aipw_point_all_we_sv[k] <- sd(aipw_point_all)
  
  se_aipw_point_oc_we_if[k] <- mean(aipw_se_oc_we)
  se_aipw_point_all_we_if[k] <- mean(aipw_se_all)
  
  
  
  M_mse_gcomp_point_oc_we[k] <- mean(mse_gcomp_point_oc_we)
  M_mse_gcomp_point_all_we[k] <- mean(mse_gcomp_point_all_we)
  
  M_mse_aipw_point_oc_we[k] <- mean(mse_aipw_point_oc_we)
  M_mse_aipw_point_all_we[k] <- mean(mse_aipw_point_all_we)
  
  M_mse_ipw_point_oc_we[k] <- mean(mse_ipw_point_oc_we)
  
  
  M_coverage_gcompOR_we[k] <- mean(coverage_gcompOR_we)
  M_coverage_gcompall_we[k] <- mean(coverage_gcompall_we)
  
  M_coverage_ipw_oc_we[k] <- mean(coverage_ipw_point_oc_we)
  
  M_coverage_aipw_oc_we[k] <- mean(coverage_aipw_point_oc_we)
  M_coverage_aipw_all[k] <- mean(coverage_aipw_point_all_we)
  
  
  
}#end threshold



bias <- c(
  bias_gcomp_point_oc_we,
  bias_gcomp_point_all_we,
  
  bias_ipw_point_oc_we,
  
  bias_aipw_point_oc_we,
  bias_aipw_point_all_we)


se_sv <- c(
  se_gcomp_point_oc_we_sv,
  se_gcomp_point_all_we_sv,
  
  se_ipw_point_oc_we_sv,
  
  se_aipw_point_oc_we_sv,
  se_aipw_point_all_we_sv)

se_if <- c(
  se_gcomp_point_oc_we_if,
  se_gcomp_point_all_we_if,
  
  se_ipw_point_oc_we_if,
  
  se_aipw_point_oc_we_if,
  se_aipw_point_all_we_if)

mse <- c(
  M_mse_gcomp_point_oc_we,
  M_mse_gcomp_point_all_we,
  
  M_mse_ipw_point_oc_we,
  
  M_mse_aipw_point_oc_we,
  M_mse_aipw_point_all_we)

coverage <- c(
  M_coverage_gcompOR_we,
  M_coverage_gcompall_we,
  
  M_coverage_ipw_oc_we, 
  
  M_coverage_aipw_oc_we, 
  M_coverage_aipw_all)


method <- c(
  rep("Gcomp-OC-WE",length(seq_2)),
  rep("Gcomp-All",length(seq_2)),
  
  rep("IPW-OC-WE",length(seq_2)),
  
  rep("AIPW-OC-WE",length(seq_2)),
  rep("AIPW-All",length(seq_2)))

pop <- rep(seq_2,5)


data_plot_1 <- data_plot_1 %>% arrange(pop, method)

data_plot_1_b <- t(round(data_plot_1[,1:5],2))
data_plot_1_b %>% kbl()

data_plot_2 <- data.frame(bias,se_sv,se_if,mse,coverage,method,pop)

data_plot_2 <- data_plot_2 %>% 
  # mutate(across(where(is.numeric), round, 2)) %>%
  arrange(pop, method)


data_plot_2_b <- t(round(data_plot_2[,1:5],2))
data_plot_2_b %>% kbl()



save.image("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24 - k2.RData")
