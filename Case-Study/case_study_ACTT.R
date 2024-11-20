rm(list=ls())

library(readr)
library(tidyverse)
library(gtsummary)


################################################################################

naive_OC <- function(data){
  
  data_Vk <- data[data$V_trA==1,]
  nk <- nrow(data_Vk)
  n1 <- sum(data_Vk$A==1)
  n0 <- sum(data_Vk$A==0)
  
  phi1 <- data_Vk$Y[data_Vk$A==1]
  
  phi0 <- data_Vk$Y[data_Vk$A==0]
  
  concATE <-  mean(phi1) - mean(phi0) 
  concATE_se <- sqrt(var(phi1)/n1 + var(phi0)/n0)
  
    return(list(concATE = concATE,
                concATE_se = concATE_se)
    )
  
}


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

#--------------------All controls
gcomp_ATE <- function(formula_m0,formula_m1,data){
  
  #----------controls
  
  #regress among controls all controls 
  m0_all_cont <- glm(formula_m0,data=data[data$A==0,])
  
  Y <- data$Y
  X <- model.matrix(formula_m0, data = data,drop = FALSE)
  
  xbeta <- (X) %*% coef(m0_all_cont) # same as predict(m0_only_cont, newdata = data)
  #e0 <- sum(data$V_trA*xbeta)/sum(data$V_trA)
  e0 <- sum(xbeta)/nrow(data)
  
  h0 <- X * (1-data$A) * as.numeric(t(Y -  xbeta) ) 
  #h1 <- data$V_trA*(xbeta - e0)
  h1 <- (xbeta - e0)
  
  #g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nrow(data)
  g11 <- ( t(X*(1-data$A)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  #g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
  g21 <- apply(-X,2,sum)/nrow(data)
  
  g22 <- nrow(data)
  g22inv <- 1/g22
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m0_all_cont))
  # same as 
  # phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2
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
  #e1 <- sum(data$V_trA*xbeta)/sum(data$V_trA) 
  e1 <- sum(xbeta)/nrow(data)
  
  # h0 <- X * data$V_trA * (data$A) * as.numeric(t(Y -  xbeta) )
  # h1 <- data$V_trA*(xbeta - e1)
  
  h0 <- X * (data$A) * as.numeric(t(Y -  xbeta) )
  h1 <- (xbeta - e1)
  
  #g11 <- ( t(X*data$V_trA*(data$A)) %*% X) / nrow(data)
  g11 <- ( t(X*(data$A)) %*% X) / nrow(data)
  g11inv <- solve(g11)
  
  g21 <- apply(-X*data$V_trA,2,sum)/nrow(data)
  
  g22 <- nrow(data)
  g22inv <- 1/g22
  
  phi_beta_hat <- g11inv %*% t(h0)
  # B <- t(h0) %*% h0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nrow(data)^2; diag(sandwich(m1_only_conc))
  # same as 
  # diag(phi_beta_hat %*% t(phi_beta_hat)/nrow(data)^2)
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





#-----------AIPW all
aipw_all_stoch <- function(formula_m0,formula_m1,formula_me,formula_vk,data,nonpara){
  
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
  
  
  
  #regress Vk
  fitVk <- glm(formula_vk, data = data, family = binomial(link = "logit")) #V_trA is NA
  
  #obtain probability of Vk
  pvk <- predict(fitVk, type="response", newdata=data) 
  
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
    #temp0 <- ( ((1-data$A)/(1-eps0) )*(data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
    temp0 <- ( ((1-data$A)/(1-eps0) ) * (pvk/pv) * (data$Y-p0_only_concC_m0_allC)  + (IVk/pv) *p0_only_concC_m0_allC )
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





################################################################################



################################################################################

f_boot <- function(data,formula_m0,formula_m1,formula_me,estimator,indices){
  
  d <- data[indices,]
  d_Vk <- d[d$V_trA==1,]
  
  #################### Naive
  naive_temp <- naive_OC(d)
  
  #################### OR-oc
  if(estimator =="or-oc"){ 
    
    or <- gcomp_OC(formula_m0,formula_m1,d_Vk)
    point <- or$mu
    se <- or$se_if
    
  }#end or-oc
  
  #################### OR-all
  if(estimator =="or-all"){ 
    
    or <- gcomp_all_WE(formula_m0,formula_m1,d)
    point <- or$mu
    se <- or$se_if
    
  }#end or-oc
  
  #################### IPW
  if(estimator =="ipw"){ 
    
    or <- IPW_OC(formula_me,d)
    point <- or$mu
    se <- or$se_if
    
  }#end ipw
  
  #################### DR-oc
  if(estimator =="dr-oc"){ 
    
    or <- aipw_OC(formula_m0,formula_m1,formula_me,d_Vk)
    point <- or$mu
    se <- or$se_if
    
  }#end ipw
  
  
  #################### DR-all
  if(estimator =="dr-all"){ 
    
    or <- aipw_all(formula_m0,formula_m1,formula_me,d,nonpara=T)
    point <- or$mu
    se <- or$se_if
    
  }#end ipw
  
  
  diff <- point - naive_temp$concATE
  
  rr <- naive_temp$concATE_se/se
  
  return(c(point,naive_temp$concATE,diff,se,naive_temp$concATE_se,rr))
  
}



################################################################################


diff_rr_boot <- function(data,R,statistic,formula_m0,formula_m1,formula_me,estimator){
  
  boot_res <- boot(data=data,
                   R=R,
                   statistic=f_boot,
                   formula_m0=formula_m0,
                   formula_m1=formula_m1,
                   formula_me=formula_me,
                   estimator=estimator)
  
  
  ####============= mu_one - mu_two
  
  mu_one <- boot_res$t0[1]
  mu_two <- boot_res$t0[2]
  mu_diff <- boot_res$t0[3]
  
  var_mu_one <- var(boot_res$t[,1])
  var_mu_two <- var(boot_res$t[,2])
  cov_mu <- cov(boot_res$t[,1],boot_res$t[,2])
  
  # var_diff <- var_mu_one + var_mu_two - 2*cov_mu
  var_diff <- var(boot_res$t[,3])
  diff <- mu_one-mu_two
  
  ci_diff <- c(diff-qnorm(0.975)*sqrt(var_diff), 
               diff+qnorm(0.975)*sqrt(var_diff))
  
  test_W <- (diff)/sqrt(var_diff)
  p_value_diff <- 2*pnorm(-abs(test_W))
  
  
  ####============= var_one / var_two
  
  sigma_one <- boot_res$t0[5]
  sigma_two <- boot_res$t0[4]
  
  var_sigma_one <- var(boot_res$t[,5])
  var_sigma_two <- var(boot_res$t[,4])
  cov_sigma <- cov(boot_res$t[,5],boot_res$t[,4])
  
  # var_rr <- rr_dm(sigma_one,sigma_two,
  #                 var_sigma_one,var_sigma_two,
  #                 cov_sigma,n)
  
  rr <- boot_res$t0[6]
  var_rr <- var(boot_res$t[,6])
  
  ci_rr <- c(rr-qnorm(0.975)*sqrt(var_rr), 
             rr+qnorm(0.975)*sqrt(var_rr))
  
  
  test_W_rr <- (boot_res$t0[6]-1)/sd(boot_res$t[,6])
  p_value_rr <- 2*pnorm(-abs(test_W_rr))
  
  
  return(list(point_est = boot_res$t0[1],
              point_naive = boot_res$t0[2],
              diff = boot_res$t0[3],
              se_diff = sqrt(var_diff),
              cov_mu = cov_mu,
              ci_diff = ci_diff,
              test_diff = test_W,
              p_value_diff = p_value_diff, 
              se_est = boot_res$t0[4],
              se_naive = boot_res$t0[5],
              rr = boot_res$t0[6],
              se_rr = sqrt(var_rr),
              cov_sigma = cov_sigma,
              ci_rr = ci_rr,
              test_rr = test_W_rr,
              p_value_rr = p_value_rr 
  )
  )
  
}


ACTT1 <- read_csv("~/Documents/NYU/Non concurrent/Ident-Estima-Parametric/data/ACTT1.csv")
ACTT1 <- ACTT1 %>% mutate(study="ACTT1")
ACTT2 <- read_csv("~/Documents/NYU/Non concurrent/Ident-Estima-Parametric/data/ACTT2.csv")
ACTT2 <- ACTT2 %>% mutate(study="ACTT2")

ACTT_full <- rbind(ACTT1,ACTT2)


#Placebo + RDV    Remdesivir is our CONTROL group

#Baricitinib + RDV is our Treatment group


ACTT_full <- ACTT_full %>% filter(TRTP != "Placebo") %>%
  mutate(trt_new = case_when(TRTP=="Placebo + RDV" | TRTP=="Remdesivir"~0,
                             T~1),
         V_trA = case_when(study=="ACTT1"~0,
                           T~1),
         E = 1:length(V_trA)/length(V_trA))
         

################################################################################
#Descriptive statistics

ACTT_table <- ACTT_full %>% ungroup %>% select(-c(USUBJID,TRTP,study)) %>%
  filter(complete.cases(ACTT_full))   %>% 
  # filter(RACE!="MULTIPLE") %>% # Positivity issue!!
  filter(STRATUM!="Mild-Moderate Disease") %>% # Positivity issue!!
  filter(COMORB1!="Unknown") %>% # Positivity issue!!
  filter(ETHNIC!="NOT REPORTED") %>%
  filter(ETHNIC!="UNKNOWN")


ACTT_table$RACE_2 <- ACTT_table$RACE
ACTT_table$RACE_2[which(ACTT_table$RACE_2 == "AMERICAN INDIAN OR ALASKA NATIVE")] <- "OTHER"
ACTT_table$RACE_2[which(ACTT_table$RACE_2 == "MULTIPLE")] <- "OTHER"
ACTT_table$RACE_2[which(ACTT_table$RACE_2 == "UNKNOWN")] <- "OTHER"
ACTT_table$RACE_2[which(ACTT_table$RACE_2 == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER")] <- "OTHER"


# ACTT_table$ETHNIC_2 <- ACTT_table$ETHNIC
# ACTT_table$ETHNIC_2[which(ACTT_table$ETHNIC_2 == "AMERICAN INDIAN OR ALASKA NATIVE")] <- "OTHER"
# ACTT_table$ETHNIC_2[which(ACTT_table$ETHNIC_2 == "MULTIPLE")] <- "OTHER"


# ACTT_table_temp <- ACTT_table %>% select(BDURSYMP, 
#                                          HYPFL, CADFL, CHFFL, CRDFL, CORFL, CLDFL, CKDFL, DIAB1FL,
#                                          DIAB2FL, OBESIFL, CANCERFL, IMMDFL, ASTHMAFL, COMORB1)
# 
# ACTT_table$sum_flags <- rowSums(ACTT_table_temp)
  
# ACTT_table %>%
#   tbl_summary(
#     by = trt_new,
#     statistic = list(
#       all_continuous() ~ "{mean} ({sd})",
#       all_categorical() ~ "{n} / {N} ({p}%)"
#     ),
#     digits = all_continuous() ~ 2,
#     missing_text = "(Missing)"
#   ) %>% add_p() %>% add_overall()





print_stuff <- function(point,se,D=2){
  
  point <- round(point,D)
  se <- round(se,D)
  
  k <- qnorm(0.975)
  ci <- round(point + c(-k*se,k*se),D)
  w <- (point-0)/se
  pvalue <- round(2*pnorm(-abs(w)),D)
  
  print <- paste("Point=",
                 point,
                 "SE=", 
                 se, 
                 "CI (", 
                 ci[1],
                 ";",
                 ci[2],
                 "); p-value",
                 pvalue, sep=" ") 
  
  return(print)
  
}


################################################################################
#STRATUM !

# formula_m0 <- Y ~ AGE + SEX + RACE + ETHNIC + BMI + REGION + STRATUM + 
#   BDURSYMP + HYPFL + CADFL + CHFFL + CRDFL + CORFL + CLDFL + CKDFL + DIAB1FL + 
#   DIAB2FL + OBESIFL + CANCERFL + IMMDFL + ASTHMAFL + COMORB1 + E

formula_m0 <- Y ~ AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E


# CORFL +  + DIAB1FL

# formula_m1 <- Y ~ AGE + SEX + RACE + ETHNIC + BMI + REGION + STRATUM + 
#   BDURSYMP + HYPFL + CADFL + CHFFL + CRDFL + CORFL + CLDFL + CKDFL + DIAB1FL + 
#   DIAB2FL + OBESIFL + CANCERFL + IMMDFL + ASTHMAFL + COMORB1  + E

formula_m1 <- formula_m0


# formula_me <- A ~ AGE + SEX + RACE + ETHNIC + BMI + REGION + STRATUM + 
#   BDURSYMP + HYPFL + CADFL + CHFFL + CRDFL + CORFL + CLDFL + CKDFL + DIAB1FL + 
#   DIAB2FL + OBESIFL + CANCERFL + IMMDFL + ASTHMAFL + COMORB1  + E

formula_me <- A ~ AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E


formula_vk <- V_trA ~ E
formula_vk <- V_trA ~ AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E

data <- ACTT_table
data$A <- data$trt_new
data$Y <- data$TTRECOV
data_Vk <- data[data$V_trA==1,]

#### - CONCATE - ####

#G-comp
gOC <- gcomp_OC(formula_m0,formula_m1,data_Vk)
gAll <- gcomp_all_WE(formula_m0,formula_m1,data)

#IPW
iOC <- IPW_OC(formula_me,data)

#AIPW
aOC <- aipw_OC(formula_m0,formula_m1,formula_me,data_Vk)
aAll <- aipw_all(formula_m0,formula_m1,formula_me,data,nonpara=T)
aAll_stoch <- aipw_all_stoch(formula_m0,formula_m1,formula_me,formula_vk,data,nonpara=T)

#Naive
nOC <- naive_OC(data)


#### - ATE - ####
nATE <- gcomp_ATE(formula_m0,formula_m1,data)

D <- 2
point <- round(c(gOC$mu,gAll$mu,iOC$mu,aOC$mu,aAll$mu,nOC$concATE),D)
sed <- c(gOC$se_if,gAll$se_if,iOC$se_if,aOC$se_if,aAll$se_if,nOC$concATE_se)
se <- round(c(gOC$se_if,gAll$se_if,iOC$se_if,aOC$se_if,aAll$se_if,nOC$concATE_se),D)
k <- qnorm(0.975)
ci_l <- round(point - c(k*se),D)
ci_u <- round(point + c(k*se),D)
w <- (point-0)/se
pvalue <- round(2*pnorm(-abs(w)),D)
ratio_se <- round(nOC$concATE_se/sed,D)
method <- c("gOC","gAll","iOC","aOC","aAll","nOC")


table_print <- data.frame(method,
                           point,
                           se,
                           paste("(",ci_l,";",ci_u,")",sep=""),
                           pvalue,
                           ratio_se
                           )

colnames(table_print) <- c("method","point","se","CI","pvalue","ratio")
table_print 


seq_est <- c("or-oc","or-all","ipw","dr-oc","dr-all")
i <- 0

res_rr <- res_rr_pv <- res_rr_lb <- res_rr_ub <- res_rr_se <- NULL
res_diff <- res_diff_pv <- res_diff_lb <- res_diff_ub <- res_diff_se <- NULL


for(est in seq_est){
i <- i + 1
  res <- diff_rr_boot(data=data,
                    R=100,
                    statistic=f_boot,
                    formula_m0=formula_m0,
                    formula_m1=formula_m1,
                    formula_me=formula_me,
                    estimator=est)
  
  res_rr[i] <- res$rr
  res_rr_se[i] <- res$se_rr
  res_rr_pv[i] <- res$p_value_rr
  res_rr_lb[i] <- res$ci_rr[1]
  res_rr_ub[i] <- res$ci_rr[2]
  
  res_diff[i] <- res$diff
  res_diff_se[i] <- res$se_diff
  res_diff_pv[i] <- res$p_value_diff
  res_diff_lb[i] <- res$ci_diff[1]
  res_diff_ub[i] <- res$ci_diff[2]

}

D <- 2
diff <- round(res_diff,D)
diff_pv <- res_diff_pv
diff_se <- round(res_diff_se,D)
ci_diff_l <- round(res_diff_lb,D)
ci_diff_u <- round(res_diff_ub,D)

rr <- round(res_rr,D)
rr_pv <- res_rr_pv
rr_se <- round(res_rr_se,D)
ci_rr_l <- round(res_rr_lb,D)
ci_rr_u <- round(res_rr_ub,D)

method <- seq_est



table_print_rr_diff <- data.frame(method,
                          diff,
                          diff_se,
                          paste("(",ci_diff_l,";",ci_diff_u,")",sep=""),
                          diff_pv,
                          rr,
                          rr_se,
                          paste("(",ci_rr_l,";",ci_rr_u,")",sep=""),
                          rr_pv
)

colnames(table_print_rr_diff) <- c("method","difference","se","CI","pvalue",
                                    "ratio","se","CI","pvalue")
table_print_rr_diff


fit_lm <- glm(Y ~ A + AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E,
              data=data)

summary(fit_lm)

fit_lm_int <- glm(Y ~ A + V_trA + AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E,
              data=data)

summary(fit_lm_int)


fit_lm_3 <- glm(Y ~ A + V_trA + E,
                data=data)

summary(fit_lm_3)


fit_lm_2 <- glm(Y ~ A + E,
              data=data)

summary(fit_lm_2)

# pate <- round(nATE$mu,D)
# seate <- round(nATE$se_if,D)
# k <- qnorm(0.975)
# ci_late <- round(pate - c(-k*seate),D)
# ci_uate <- round(pate + c(-k*seate),D)
# wate <- (pate-0)/seate
# pvalueate <- round(2*pnorm(-abs(wate)),D)
# 
# table_ate <- data.frame(pate,
#            seate,
#            paste("(",ci_late,";",ci_uate,")",sep=""),
#            pvalueate
# )
# 
# colnames(table_ate) <- c("point","se","CI","pvalue")
# table_ate
# 
# 
# 
# 
# ACTT2 %>% select(TTRECOV,TRTP) %>%
#   tbl_summary(
#     by = TRTP,
#     statistic = list(
#       all_continuous() ~ "{mean} ({sd})",
#       all_categorical() ~ "{n} / {N} ({p}%)"
#     ),
#     digits = all_continuous() ~ 2,
#     missing_text = "(Missing)"
#   ) %>% add_p() %>% add_overall()
# 
# 
# 
# print("ADD COVERAGE")









