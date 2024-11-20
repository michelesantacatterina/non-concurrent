
# library(parallel)
# library(SuperLearner)
# library(rpart)
# library(MASS)
# library(polspline)
# library(nnet)
# library(gam)
# library(splines)
# library(plotrix)
# library(earth)
# library(BayesTree)

rm(list=ls())

library(readr)
library(tidyverse)
library(gtsummary)
library(rARPACK) # if computing eigenvalues with this package (function eigs)

mmd.test = function(R,S,D.R,D.S,sig.meth='eig',num.reps=1e4,return.cutoff=FALSE){
  n = length(R)
  
  D.R.mat1 = matrix(rep(D.R,n),nrow=n)
  D.R.mat2 = matrix(rep(D.R,each=n),nrow=n)
  
  D.S.mat1 = matrix(rep(D.S,n),nrow=n)
  D.S.mat2 = matrix(rep(D.S,each=n),nrow=n)
  
  R.mat1 = matrix(rep(R,n),nrow=n)
  R.mat2 = matrix(rep(R,each=n),nrow=n)
  
  S.mat1 = matrix(rep(S,n),nrow=n)
  S.mat2 = matrix(rep(S,each=n),nrow=n)
  
  EE = ((2*(R.mat1-R.mat2)*(D.R.mat2-D.R.mat1) + 1 - (4*(R.mat1-R.mat2)^2-2)*D.R.mat1*D.R.mat2)*exp(-(R.mat1-R.mat2)^2)
        - ((2*(S.mat1-R.mat2)*(D.R.mat2-D.S.mat1) + 1 - (4*(S.mat1-R.mat2)^2-2)*D.S.mat1*D.R.mat2)*exp(-(S.mat1-R.mat2)^2))
        - ((2*(R.mat1-S.mat2)*(D.S.mat2-D.R.mat1) + 1 - (4*(R.mat1-S.mat2)^2-2)*D.R.mat1*D.S.mat2)*exp(-(R.mat1-S.mat2)^2))
        + (2*(S.mat1-S.mat2)*(D.S.mat2-D.S.mat1) + 1 - (4*(S.mat1-S.mat2)^2-2)*D.S.mat1*D.S.mat2)*exp(-(S.mat1-S.mat2)^2))
  
  
  if(sig.meth=='eig'){
    line.means = rowMeans(EE)
    EE.ctrd = EE - matrix(rep(line.means,n),nrow=n) - matrix(rep(line.means,each=n),nrow=n) + matrix(rep(mean(line.means),n^2),nrow=n)
    num.eigs = min(200,n)
    require('rARPACK')
    
    tmp = eigen(EE.ctrd)$values/n
    tmp = tmp[Im(tmp)==0]
    tmp = Re(tmp)
    tmp = sort(tmp,decreasing=TRUE)
    if(length(tmp)>num.eigs){
      tmp = tmp[1:num.eigs]
    }
    num.pos.eigs = min(num.eigs,length(tmp))
    
    draws=c(matrix(rnorm(num.reps*num.pos.eigs)^2-1,nrow=num.reps,ncol=num.pos.eigs)%*%cbind(tmp[1:num.pos.eigs]))
  }
  
  # U-statistic
  diag(EE) = 0
  est = (rbind(rep(1/(n-1),n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
  
  if(sig.meth=='eig'){
    pval = mean(draws>n*est)
  } else if(sig.meth=='var'){
    pval = pchisq(est/(2*var(D.R)/n)+1,df=1,lower.tail=FALSE)
  }
  
  return(if(!return.cutoff){
    c(est,pval)
  }else{
    c(est,pval,
      if(sig.meth=='eig'){
        quantile(draws,0.95)
      }else{
        2*var(D.R)*(qchisq(0.95,df=1)-1)})})
}

# Tests if A is used in the regression function, i.e. if
# E[Y|A,W] = E[Y|W] almost surely, or equivalently (under positivity) if
# E[Y|A=1,W] = E[Y|A=0,W] almost surely
# Returns an estimate of Psi, and a p-value
# est.psi.prob = function(formula_m0,data,sig.meth='var',est.g=FALSE,dynamic.bandwidth=FALSE){
#   
#   
#   #--- E[Y|A=0, W=w, E=e, Vk=1] = E[Y|A=0, W=w, E=e]
#   
#   data_Vk <- data[data$V_trA==1,]
#   n <- nrow(data)
#   
#   #--------------------------regress among all controls
#   m0_allC <- glm(formula_m0,data=data)
#   m0_onlyC <- glm(formula_m0,data=data_Vk)
#   # 
#   #predict among only concurrent controls
#   Qbar.est_all_Vk1 <- predict(m0_allC, newdata = data)
#   Qbar.est_only_Vk1 <- rep(0,n)
#   Qbar.est_only_Vk1[data$V_trA==1] <- predict(m0_onlyC, newdata = data_Vk)
#   
#   gg <- mean(data$V_trA)
# 
#   
#   
#   # 
#   # R = Qbar.est_only_Vk1
#   # D.R = data$V_trA/gg * (data$Y-R)
#   # D.R <- D.R[data$V_trA==1]
#   # 
#   # R <- R[data$V_trA==1]
#   # 
#   # S = Qbar.est_all_Vk1
#   # D.S = (data$Y-S)
#   # D.S <- D.S[data$V_trA==1]
#   # 
#   # S <- S[data$V_trA==1]
#   # 
#   # 
#   # if(dynamic.bandwidth){
#   #   inv.bw = 1/median(dist(c(R + D.R,S+D.S)))
#   #   print(1/inv.bw)
#   # }
#   # R = inv.bw*R
#   # S = inv.bw*S
#   # D.R = inv.bw*D.R
#   # D.S = inv.bw*D.S
#   # 
#   # return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
# 
#   
#   
#   # Plug-in estimate of blip
#   R = Qbar.est_only_Vk1 - Qbar.est_all_Vk1
#   R <- R[data$V_trA==1]
#   S = rep(0,n)
#   S = rep(0,nrow(data_Vk))
#   D.R_only <- data$V_trA/gg * (data$Y-Qbar.est_only_Vk1)
#   #D.R_all  <- (1-data$V_trA)/(1-gg) * (data$Y-Qbar.est_all_Vk1)
#   D.R_all  <- (data$Y-Qbar.est_all_Vk1)
#   D.R <- D.R_only - D.R_all
#   D.R <- D.R[data$V_trA==1]
#   D.S = rep(0,n)
#   D.S = rep(0,nrow(data_Vk))
#   
#   if(dynamic.bandwidth){
#     med = median(dist(c(R + D.R,S+D.S)))
#     print(med)
#     R = R/med
#     D.R = D.R/med
#   }
# 
#   return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
# }







# Tests if A is used in the regression function, i.e. if
# E[Y|A,W] = E[Y|W] almost surely, or equivalently (under positivity) if
# E[Y|A=1,W] = E[Y|A=0,W] almost surely
# Returns an estimate of Psi, and a p-value
est.psi.prob = function(formula_m0,data,sig.meth='var',est.g=FALSE,dynamic.bandwidth=FALSE){
  
  
  #--- E[Y|A=0, W=w, E=e, Vk=1] = E[Y|A=0, W=w, E=e] for all e s.t. Vk=1
  
  data_Vk <- data[data$V_trA==1,]
  n <- nrow(data)
  
  
  #-----Outcome part
  
      #--------------------------regress among all controls
      m0_allC <- glm(formula_m0,data=data[data$trt_new==0,])
      
      #--------------------------regress among only concurrent
      m0_onlyC <- glm(formula_m0,data=data_Vk[data_Vk$trt_new==0,])
      
      #prediction 
      #predict among only concurrent controls
      Qbar.est_only_Vk1                 <- rep(0,n)
      Qbar.est_only_Vk1[data$V_trA==1]  <- predict(m0_onlyC, newdata = data_Vk)
      
      #predict for all e s.t. Vk=1 
      #Recall that "for all e s.t. Vk=1" is the same as "e in the Vk=1 population" 
      #summary(data$E[which(data$V_trA==1)]) is the same as
      #summary(data_Vk$E)
      # Qbar.est_all_Vk1                  <- rep(0,n)
      # Qbar.est_all_Vk1[data$V_trA==1]   <- predict(m0_allC, newdata = data_Vk)
      Qbar.est_all_Vk1                  <- predict(m0_allC, newdata = data) # -----------------------------------> check
      
  #-----Treatment/concurrent part    
  #E[Y|A=0, W=w, E=e, Vk=1] = ( 1[A=0,Vk=1] / (P(A=0,Vk=1 | W, E)) ) Y
  #by chain rule: P(A=0,Vk=1 | W, E) = P(A=0 | Vk=1, W, E) P(Vk=1 | W, E)    
      
  #compute P(A=0 | Vk=1, W, E)
      fitprA <- glm( I(1-trt_new) ~ AGE + SEX + RACE_2 + BMI + REGION + STRATUM + COMORB1 + E
                     , data = data_Vk, family = "binomial" )
      temp_prA0 <- predict(fitprA,type="response")
      prA0 <- rep(1,n)
      prA0[data$V_trA==1] <- temp_prA0
      # prA0[data$V_trA==1] <- 0.5
      
  #compute P(Vk=1 | W, E)    
      # fitprV <- glm( V_trA ~ AGE + SEX + RACE_2 + BMI + REGION + STRATUM + COMORB1
      #                , data = data, family = "binomial" )
      # prV <- predict(fitprV,type="response")
      prV <- mean(data$V_trA) # I do not model this one -----------------------------------> check because P(Vk=1 | W, E) = I[E>t] = Vk
      
      
  #E[Y|A=0, W=w, E=e] for all e s.t. Vk=1 = ( 1[A=0] / (P(A=0 | W, E)) ) Y
  #this is the same as prA0
  
  # Plug-in estimate of blip
  Qbar.est_all_Vk1[data$V_trA == 0]  <- 0 # -----------------------------------> check (i believe this is what we mean with "for all e s.t. Vk=1")
  
  R = Qbar.est_only_Vk1 - Qbar.est_all_Vk1
  S = rep(0,n)
  D.R_only <- ( (1-data$trt_new) * data$V_trA)/(prA0*prV) * (data$Y-Qbar.est_only_Vk1)
  #D.R_only <- ( (1-data$trt_new))/(prA0) * (data$Y-Qbar.est_only_Vk1)
  D.R_all  <- ( (1-data$trt_new)/prA0 ) * (data$Y-Qbar.est_all_Vk1)
  D.R_all[data$V_trA == 0]  <- 0
  D.R <- D.R_only - D.R_all
  D.S = rep(0,n)
  
  if(dynamic.bandwidth){
    med = median(dist(c(R + D.R,S+D.S)))
    print(med)
    R = R/med
    D.R = D.R/med
  }
  
  return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
}










# Tests if A is used in the regression function, i.e. if
# E[Y|A,W] = E[Y|W] almost surely, or equivalently (under positivity) if
# E[Y|A=1,W] = E[Y|A=0,W] almost surely
# Returns an estimate of Psi, and a p-value
est.psi.prob_with_Ivan = function(formula_m0,data,sig.meth='var',est.g=FALSE,dynamic.bandwidth=FALSE){
  
  
  #--- E[Y|A=0, W=w, E=e, Vk=1] = E[Y|A=0, W=w, E=e] for all e s.t. Vk=1
  
  data_Vk <- data[data$V_trA==1,]
  n <- nrow(data)
  
  
  #-----Outcome part
  
  #--------------------------regress among all controls
  m0_allC <- glm(formula_m0,data=data[data$trt_new==0,])
  
  #--------------------------regress among only concurrent
  m0_onlyC <- glm(formula_m0,data=data_Vk[data_Vk$trt_new==0,])
  
  #prediction 
  #predict among only concurrent controls
  # Qbar.est_only_Vk1                 <- rep(0,n)
  Qbar.est_only_Vk1  <- predict(m0_onlyC, newdata = data_Vk)
  
  #predict for all e s.t. Vk=1 
  #Recall that "for all e s.t. Vk=1" is the same as "e in the Vk=1 population" 
  #summary(data$E[which(data$V_trA==1)]) is the same as
  #summary(data_Vk$E)
  Qbar.est_all_Vk1   <- predict(m0_allC, newdata = data_Vk)

  
  #-----Treatment/concurrent part    
  #E[Y|A=0, W=w, E=e, Vk=1] = ( 1[A=0,Vk=1] / (P(A=0,Vk=1 | W, E)) ) Y
  #by chain rule: P(A=0,Vk=1 | W, E) = P(A=0 | Vk=1, W, E) P(Vk=1 | W, E)    
  
  #compute P(A=0 | Vk=1, W, E)
  fitprA <- glm( I(1-trt_new) ~ AGE + SEX + RACE_2 + BMI + REGION + STRATUM + COMORB1 + E
                 , data = data_Vk, family = "binomial" )
  temp_prA0 <- predict(fitprA,type="response")
  # prA0 <- rep(1,n)
  prA0 <- temp_prA0
  # prA0[data$V_trA==1] <- 0.5
  
  #compute P(Vk=1 | W, E)    
  # fitprV <- glm( V_trA ~ AGE + SEX + RACE_2 + BMI + REGION + STRATUM + COMORB1
  #                , data = data, family = "binomial" )
  # prV <- predict(fitprV,type="response")
  prV <- mean(data$V_trA) # I do not model this one -----------------------------------> check because P(Vk=1 | W, E) = I[E>t] = Vk
  prV <- 1
  
  #E[Y|A=0, W=w, E=e] for all e s.t. Vk=1 = ( 1[A=0] / (P(A=0 | W, E)) ) Y
  fitprA_noVk <- glm( I(1-trt_new) ~ AGE + SEX + RACE_2 + BMI + REGION + STRATUM + COMORB1 + E
                 , data = data, family = "binomial" )
  temp_prA0_noVk <- predict(fitprA_noVk,newdata= data_Vk, type="response")
  prA0_noVk <- temp_prA0_noVk
  
  # Plug-in estimate of blip
  #Qbar.est_all_Vk1[data$V_trA == 0]  <- 0 # -----------------------------------> check (i believe this is what we mean with "for all e s.t. Vk=1")
  
  R = Qbar.est_only_Vk1 - Qbar.est_all_Vk1
  S = rep(0,length(R))
  D.R_only <- ( (1-data_Vk$trt_new) )/(prA0) * (data_Vk$Y-Qbar.est_only_Vk1)
  #D.R_only <- ( (1-data$trt_new))/(prA0) * (data$Y-Qbar.est_only_Vk1)
  D.R_all  <- ( (1-data_Vk$trt_new)/prA0_noVk ) * (data_Vk$Y-Qbar.est_all_Vk1)
  #D.R_all[data$V_trA == 0]  <- 0
  D.R <- D.R_only - D.R_all
  D.S = rep(0,length(D.R))
  
  if(dynamic.bandwidth){
    med = median(dist(c(R + D.R,S+D.S)))
    print(med)
    R = R/med
    D.R = D.R/med
  }
  
  return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
}



#####---------------------------------------------------------------------------


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

ACTT_full$Erand <- runif(nrow(ACTT_full),0,1)


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



# ACTT_table_A0 <- ACTT_table %>% filter(trt_new == 0)


ACTT_table  %>% select(c(AGE, SEX,  RACE_2,   BMI,  REGION,  STRATUM,  COMORB1,  E, V_trA, trt_new)) %>%
  tbl_summary(
    by = trt_new,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing_text = "(Missing)"
  ) %>% add_p() %>% add_overall()


#I only consider among controls!
#data <- ACTT_table %>% filter(trt_new==0)

data <- ACTT_table

################################################################################


data$Y <- data$TTRECOV


formula_m0 <- Y ~ AGE + SEX + RACE_2 +  BMI + REGION + STRATUM + COMORB1 + E
#formula_m0 <- Y ~ E
#formula_m0 <- Y ~ Erand

mmd.var = est.psi.prob_with_Ivan(formula_m0,data,sig.meth='var',est.g=FALSE)
mmd.var

mmd.eig = est.psi.prob_with_Ivan(formula_m0,data,sig.meth='eig',est.g=FALSE)
mmd.eig

