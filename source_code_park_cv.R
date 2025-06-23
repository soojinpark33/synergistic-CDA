## Clear workspace
rm(list = ls())

## Load packages
library(dplyr)
library(gt) #to draw table to show results at the end
library(parallel)
library(MASS)
library(rbw)
library(caret)
library(xgboost)
library(purrr)

#-----------------------------------------------------------------------------------------------#
# Simulation Study
#-----------------------------------------------------------------------------------------------#
# Data Generation
## Set important parameters
beta_m = 0.7
beta_rm = -0.2
beta_rma = 0.5

## Generate population data
invLogit = function(x){
  exp(x)/(1 + exp(x))
}

N = 1e+06
Xvec <- mvrnorm(N, mu = rep(0, 3), Sigma = diag(3))
X.mis <- cbind(exp(Xvec[,1]/2), Xvec[,2]*(1+exp(Xvec[,1]))^(-1)+10,
               (Xvec[,1]*Xvec[,3]/25+.6)^3)
C = rbinom(N, 1, 0.4)
pR1_C = invLogit(.5 - 0.5*C) # P(R=1|C)
R = rbinom(N, 1, pR1_C)
pA1_RC = invLogit(- .8 + 1*R + 1.5*C + Xvec[,1] + 0.2*Xvec[,2]-0.5*Xvec[,3] + R*Xvec[,3]) 
A = rbinom(N, 1, pA1_RC)
Z = -.5+.2*R + 1.2*A + .5*C -.5* Xvec[,1] +  0.7*Xvec[,2]+0.5*Xvec[,3] + rnorm(N, 0, 1)
pM1_RXC = invLogit(-1+2*R + 0.2*A + 1*C - Xvec[,1] - 0.2*Xvec[,2]+1.5*Xvec[,3] + R*Xvec[,2] +.5*Z) # P(M=1|R,X,C,U)
M = rbinom(N, 1, pM1_RXC)
Y = 1 - 0.5*R + beta_m*M + beta_rm*R*M + beta_rma*R*M*A + A + C -Xvec[,1] +0.5*Xvec[,2]-.5*Xvec[,3] -.5*Z + rnorm(N, 0, 1) 
popdata = data.frame( Y, R, M, A, C, Xvec, X.mis, Z)

## True value of IIE
YG0 <- (1 - 0.5*1 + beta_m*mean(M[R==0&C==0]) + beta_rm*1*mean(M[R==0&C==0]) + 
          beta_rma*1*mean(M[R==0&C==0])*mean(A[R==0&C==0]) + 1*mean(A[R==0&C==0]) + 
          mean(-Xvec[R==1&C==0,1] + 0.5*Xvec[R==1&C==0,2]-0.5*Xvec[R==1&C==0,3])-
          .5*(-.5+.2*1 + 1.2*mean(A[R==0&C==0]) +mean(c(-.5, 0.7, 0.5) %*% t(Xvec[R == 1 & C == 0, ]))))
delta_IIE = mean(Y[R==1&C==0]) - YG0
zeta_IIE = YG0 - mean(Y[R==0&C==0])


## simulation starts:
#::::::::::::::::::::::::::
# Function 1: GLM analysis
#::::::::::::::::::::::::::

  # Simulation 0: correct specification
  # Simulation 1: correct imputation only
  # Simulation 2: correct weighting only
  # Simulation 3: correct outcome and A specification
  # Simulation 4: Incorrect for all 

glm.func <- function(sim, data){ 
#:::::::::::: Pure Imputation
  #Fit models for Y, A, and M
      if (sim == 0 | sim == 1 | sim == 3){
        fit.y <- lm(Y~R+ A+ M + X1+ X2+ X3 + Z + C+ R:M:A + R:M, data=data)
      }
      else if (sim == 2 | sim == 4){
        fit.y <- lm(Y~R+ A+ M + X1.1+ X2.1+ X3.1+ Z+ C, data=data)
      }
      else{
        print("ERROR: choose between 0 to 4 for s0 to s4")
      }
      
    fit.m0 <- glm(M ~ C, data = subset(data, R == 0), family = binomial(logit))
    fit.a0 <- glm(A ~ C, data = subset(data, R == 0), family = binomial(logit))
    
  # Compute Observed outcome of Y given C
    wy0 <- lm(Y ~ C, data = subset(data, R == 0)) 
    wy1 <- lm(Y ~ C, data = subset(data, R == 1))
    
  # Predict M and A
    pm0 <- predict(fit.m0, data, type = "response")
    pa0 <- predict(fit.a0, data, type = "response")
    
  # Compute Mu after incorporating predicted values of M and A
    dat_1<- data
    dat_1[,"M"] <- rbinom(length(pm0), size = 1, prob = pm0)
    dat_1[,"A"] <- rbinom(length(pa0), size = 1, prob = pa0)
    data$muldm <- predict(fit.y, newdata = dat_1)
    
  # Compute nu
      if (sim == 0 | sim == 1){
        a <- lm(muldm ~ A + C + X1+ X2+ X3, data = subset(data, R == 1))
      }
      else if ( sim == 2| sim == 3 | sim == 4){
        a <- lm(muldm ~ A + C + X1.1+ X2.1+ X3.1, data = subset(data, R == 1))
      }
      else{
        print("ERROR: choose between 0 to 4 for s0 to s4")
      }

    nu <- predict(a, newdata= subset(dat_1,R==1))
    b <- lm(nu ~ C, data=subset(data,R==1))
    
  # Calculate the estimates of the disparity reduction and remaining
    delta_imp <- wy1$coef[1] - b$coef[1]
    zeta_imp <- b$coef[1] - wy0$coef[1]  
    
#::::::::::::: Pure Weighting
  # Fit regression models for M & A
    if (sim == 0 | sim == 2 ){
      fit.m1 <- glm(M ~ A + X1 + X2 + X3 + Z + C, data = subset(data, R == 1), family = binomial(logit))
    }
    else if (sim == 1 | sim == 3| sim == 4){
      fit.m1 <- glm(M ~ A + X1.1 + X2.1 + X3.1 + Z+ C, data = subset(data, R == 1), family = binomial(logit))
    }
    else{
      print("ERROR: choose between 0 to 4 for s0 to s4")
    }    
    
    if (sim == 0 | sim == 2| sim == 3){
      fit.a1 <- glm(A ~ X1 + X2 + X3 + C, data = subset(data, R == 1), family = binomial(logit))
    }
    else if (sim == 1 | sim == 4){
      fit.a1 <- glm(A ~ X1.1 + X2.1 + X3.1 + C, data = subset(data, R == 1), family = binomial(logit))
    }
    else{
      print("ERROR: choose between 0 to 4 for s0 to s4")
    }        
    
  # Calculate the conditional probabilities of M and A
    pm1 <- predict(fit.m1, data, type = "response")
    pa1 <- predict(fit.a1, data, type = "response")
    
    
  # Calculate weights
    WM <- rep(NA, nrow(data))
    ind11 <- data$R == 1 & data$M == 1
    ind10 <- data$R == 1 & data$M == 0
    WM[ind11] <- pm0[ind11]/pm1[ind11]
    WM[ind10] <- (1 - pm0[ind10])/(1-pm1[ind10])
    WM[data$R == 0] <- 0
    data$WM <- WM
    
    WA <- rep(NA, nrow(data))
    ind11a <- data$R == 1 & data$A == 1
    ind10a <- data$R == 1 & data$A == 0
    WA[ind11a] <- pa0[ind11a]/pa1[ind11a]
    WA[ind10a] <- (1-pa0[ind10a])/(1-pa1[ind10a])
    WA[data$R == 0] <- 0
    data$WA <- WA
    
  # Calculate weighted y
    wmu1xdm <- lm(Y ~ C, weight = WA*WM, data = subset(data, R == 1))
    
  # Calculate the estimates of the disparity reduction and remaining
    delta_wgt = wy1$coef[1] - wmu1xdm$coef[1]
    zeta_wgt = wmu1xdm$coef[1] - wy0$coef[1]
    
#::::::::::::: Imputation and Weighting
    wa_mu <- lm(muldm ~ C, weight = WA, data = subset(data, R == 1))
    delta_iw = wy1$coef[1] - coef(wa_mu)[1]
    zeta_iw = coef(wa_mu)[1] -wy0$coef[1]   
    
#::::::::::::: Triply Robust 
    wa_nu <- lm(nu ~  C , weight = WA, data = subset(data, R == 1))
    data$mu <- predict(fit.y, newdata=data)
    wawm_mu <- lm(mu ~ C, weight = WM*WA, data = subset(data, R == 1)) 
    
    delta_tr = wy1$coef[1] -(b$coef[1]+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))
    zeta_tr = (mean(nu)+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))-wy0$coef[1] 
    
    return(c(delta_imp, zeta_imp, delta_wgt, zeta_wgt, delta_iw, zeta_iw, delta_tr, zeta_iw))
}

#::::::::::::::::::::::::::
# Function 2: xg analysis
#::::::::::::::::::::::::::

xg.func <- function(sim, data){ 
  
  #:::::::::::: Pure Imputation
  #Fit model for Y
  selected_columns <- if (sim == 0 | sim == 1 | sim == 3) {
    c("R", "A", "M", "X1", "X2", "X3", "Z", "C")
  } else {
    c("R", "A", "M", "X1.1", "X2.1", "X3.1", "Z","C")
  }
  X <- as.matrix(data %>% dplyr::select(all_of(selected_columns)))
  Y <- data$Y
  
  # Set up training data for XGBoost (outcome model)
  dtrain_y <- xgb.DMatrix(data = X, label = Y)
  
  # Train XGBoost model for outcome prediction
  fit.y <- xgb.cv(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6),
                  data = dtrain_y, 
                  nfold = 5,
                  nrounds = 100, 
                  early_stopping_rounds = 10,
                  verbose = 0)
  best_nrounds <- fit.y$best_iteration
  
  fit.y <- xgboost(data = dtrain_y, 
                   max_depth = 6, 
                   nrounds = best_nrounds, 
                   objective = "reg:squarederror", 
                   lambda = 1, 
                   alpha = 0.1, 
                   verbose = 0)
  #Fit models for M & A
  fit.m0 <- glm(M ~ C, data = subset(data, R == 0), family = binomial(logit))
  fit.a0 <- glm(A ~ C, data = subset(data, R == 0), family = binomial(logit))

  # Compute observed outcomes for each group
  wy0 <- lm(Y ~ C, data = subset(data, R == 0)) 
  wy1 <- lm(Y ~ C, data = subset(data, R == 1)) 
  
  # Compute predicted values of M and A
  pm0 <- predict(fit.m0, data, type = "response")
  pa0 <- predict(fit.a0, data, type = "response")

  # Compute mu
  dat_1<- data
  PredictM <- rbinom(length(pm0), size = 1, prob = pm0)
  dat_1[, "M"] <- (PredictM)
  PredictA <- rbinom(length(pa0), size = 1, prob = pa0)
  dat_1[, "A"] <- (PredictA)
  data$muldm <- predict(fit.y, newdata = as.matrix(dat_1 %>% 
                                                     dplyr::select(all_of(selected_columns))))
  
  # Compute nu
  selected_columns1 <- if (sim == 0 | sim == 1) {
    c( "A", "X1", "X2", "X3", "C")
  } else {
    c( "A", "X1.1", "X2.1", "X3.1", "C")
  }
  
  A_final <- as.matrix(data %>% filter(R == 1) %>% dplyr::select(all_of(selected_columns1)))
  Y_final <- data %>% filter(R == 1) %>% pull(muldm)
  
  # Train the final outcome model on data with mu
  dtrain_final <- xgb.DMatrix(data = A_final, label = Y_final)
  fit_final <- xgb.cv(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6),
                      data = dtrain_final, 
                      nfold = 5,
                      nrounds = 100, 
                      early_stopping_rounds = 10,
                      verbose = 0)
  best_nrounds <- fit_final$best_iteration
  
  fit_final <- xgboost(data = dtrain_final, 
                       max_depth = 6, 
                       nrounds = best_nrounds, 
                       objective = "reg:squarederror", 
                       lambda = 1, 
                       alpha = 0.1, 
                       verbose = 0)

  X_A1C0 <- as.matrix(dat_1 %>% dplyr::select(all_of(selected_columns1)))
  data$nu <- predict(fit_final, newdata= X_A1C0)
  b <- lm(nu ~ C, data=subset(data,R==1))
  
  # Calculate the estimates of the disparity reduction and remaining
  delta_imp <- wy1$coef[1] - b$coef[1]
  zeta_imp <- b$coef[1] - wy0$coef[1]  
  
  
  #::::::::::::: Pure Weighting
  # Fit regression models for M
  # Prepare data matrices for XGBoost (for R == 1 subset)
  data_r1 <- data %>% filter(R == 1)
  
  # covariate matrix for predicting M 
  selected_columns2 <- if (sim == 0 | sim == 2) {
    c( "A", "X1", "X2", "X3", "Z", "C")
  } else {
    c( "A", "X1.1", "X2.1", "X3.1", "Z", "C")
  }
  M_covariates <- as.matrix(data_r1 %>% dplyr::select(all_of(selected_columns2)))
  M_outcome <- data_r1$M
  
  # Train XGBoost model to predict M
  dtrain_m <- xgb.DMatrix(data = M_covariates, label = M_outcome)
  pos_weight1 <- sum(data_r1$M == 0) / sum(data_r1$M == 1)
  fit.m1 <- xgb.cv(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 6),
                   data = dtrain_m, 
                   nfold = 5,
                   nrounds = 100, 
                   early_stopping_rounds = 10,
                   verbose = 0)
  best_nrounds <- fit.m1$best_iteration
  
  fit.m1 <- xgboost(data = dtrain_m, 
                    max_depth = 6, 
                    nrounds = best_nrounds, 
                    scale_pos_weight = pos_weight1,
                    objective = "binary:logistic", 
                    #lambda = 1, 
                    #alpha = 0.1, 
                    verbose = 0)
  
  # covariate matrix for predicting X
  selected_columns3 <- if (sim == 0 | sim == 2 | sim == 3) {
    c(  "X1", "X2", "X3", "C")
  } else {
    c(  "X1.1", "X2.1", "X3.1", "C")
  }
  
  A_covariates <- as.matrix(data_r1 %>% dplyr::select(all_of(selected_columns3)))
  A_outcome <- data_r1$A
  
  # Train XGBoost model to predict X
  pos_weight <- sum(data_r1$A == 0) / sum(data_r1$A == 1)
  dtrain_a <- xgb.DMatrix(data = A_covariates, label = A_outcome)
  fit.a1 <- xgb.cv(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 6),
                   data = dtrain_a, 
                   nfold = 5,
                   nrounds = 100, 
                   early_stopping_rounds = 10,
                   verbose = 0)
  best_nrounds <- fit.a1$best_iteration
  fit.a1 <- xgboost(data = dtrain_a, 
                    max_depth = 6, 
                    nrounds = best_nrounds, 
                    scale_pos_weight = pos_weight,
                    objective = "binary:logistic", 
                    #lambda = 1, 
                    #alpha = 0.1, 
                    verbose = 0)
  
  # Predict probabilities for M and A for the entire dataset (even outside R == 1)
  X_full_M <- as.matrix(data %>% dplyr::select(all_of(selected_columns2)))
  
  ## sim switch for x full covariates
  X_full_A <- as.matrix(data %>% dplyr::select(all_of(selected_columns3)))
  
  # Predict probabilities for M using XGBoost
  pm1 <- predict(fit.m1, newdata = X_full_M, type = "response")
  
  # Predict probabilities for X using XGBoost
  pa1 <- predict(fit.a1, newdata = X_full_A, type = "response")
  
  # Initialize weights
  WM <- rep(NA, nrow(data))
  WA <- rep(NA, nrow(data))
  
  # Subset indices for calculating weights
  ind11 <- data$R == 1 & data$M == 1
  ind10 <- data$R == 1 & data$M == 0
  
  # Calculate weights for WM based on predicted probabilities
  WM[ind11] <- pm0[ind11] / pm1[ind11]
  WM[ind10] <- (1 - pm0[ind10]) / (1 - pm1[ind10])
  WM[data$R == 0] <- 0
  data$WM <- WM
  
  # Subset indices for calculating weights
  ind11a <- data$R == 1 & data$A == 1
  ind10a <- data$R == 1 & data$A == 0
  
  # Calculate weights for WA based on predicted probabilities
  WA[ind11a] <- pa0[ind11a] / pa1[ind11a]
  WA[ind10a] <- (1-pa0[ind10a]) / (1 - pa1[ind10a])
  WA[data$R == 0] <- 0
  data$WA <- WA
  
  # Calculate weighted Y using the weighted least squares (lm)
  wmu1xdm <- lm(Y ~ C, weights = WA * WM, data = subset(data, R == 1))
  delta_wgt <- wy1$coef[1] - wmu1xdm$coef[1]
  zeta_wgt <- wmu1xdm$coef[1] - wy0$coef[1]
  
  #::::::::::::: Imputation and Weighting
  wa_mu <- lm(muldm ~ C, weight = WA, data = data)
  delta_iw = wy1$coef[1] - coef(wa_mu)[1]
  zeta_iw = coef(wa_mu)[1] -wy0$coef[1]   
  
  #::::::::::::: Triply Robust 
  wa_nu <- lm(nu ~  C , weights=WA, data = subset(data, R == 1))
  data$mu <- predict(fit.y, newdata=dtrain_y)
  wawm_mu <- lm(mu ~ C, weight = WM*WA, data = subset(data, R == 1)) 
  
  
  delta_tr = wy1$coef[1] -(b$coef[1]+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))
  zeta_tr = (b$coef[1]+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))-wy0$coef[1] 
  
  return(c(delta_imp, zeta_imp, delta_wgt, zeta_wgt, delta_iw, zeta_iw, delta_tr, zeta_iw))
}
#::::::::::::::::::::::::::
# Function 3: xg anlysis with cf
#::::::::::::::::::::::::::

xg.func.cf <- function(sim, data, K = 2) {
  # number of cross-fitting folds
  folds <- createFolds(data$Y, K)
  main_list <- vector(mode = "list", K)

  for (k in 1:K) {
    # Split data
    train_data <- data[-folds[[k]], ]
    
    # pure Imputation ::::::::
    # Define X matrix for Y
    selected_columns <- if (sim == 0 | sim == 1 | sim == 3) {
      c("R", "A", "M", "X1", "X2", "X3", "Z", "C")
    } else {
      c("R", "A", "M", "X1.1", "X2.1", "X3.1", "Z","C")
    }
    X <- as.matrix(train_data %>% dplyr::select(all_of(selected_columns)))
    Y <- as.numeric(train_data$Y)
    
    # Train outcome model for imputation
    dtrain_y <- xgb.DMatrix(data = X, label = Y)
    fit.y <- xgb.cv(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6),
                    data = dtrain_y, 
                    nfold = 5,
                    nrounds = 100, 
                    early_stopping_rounds = 10,
                    verbose = 0)
    best_nrounds <- fit.y$best_iteration
    
    fit.y <- xgboost(data = dtrain_y, 
                     max_depth = 6, 
                     nrounds = best_nrounds, 
                     objective = "reg:squarederror", 
                     lambda = 1, 
                     alpha = 0.1, 
                     verbose = 0)
    
    # Fit models for M and A
    fit_m0 <- glm(M ~ C, data = subset(train_data, R == 0), family = binomial(logit))
    fit_a0 <- glm(A ~ C, data = subset(train_data, R == 0), family = binomial(logit))
    pm0 <- predict(fit_m0, train_data, type = "response")
    pa0 <- predict(fit_a0, train_data, type = "response")
    
    # Predict outcome
    train_data_1<- train_data
    PredictM <- rbinom(length(pm0), size = 1, prob = pm0)
    train_data_1[, "M"] <- (PredictM)
    PredictA <- rbinom(length(pa0), size = 1, prob = pa0)
    train_data_1[, "A"] <- (PredictA)

    # For latter use (weighting-then-imputation estimator)
    fit_m0_d <- glm(M ~ C, data = subset(data, R == 0), family = binomial(logit))
    fit_a0_d <- glm(A ~ C, data = subset(data, R == 0), family = binomial(logit))
    
    pm0_d <- predict(fit_m0_d, data, type = "response")
    pa0_d <- predict(fit_a0_d, data, type = "response")
    data_1<- data
    PredictM_d <- rbinom(length(pm0_d), size = 1, prob = pm0_d)
    data_1[, "M"] <- (PredictM_d)
    PredictA_d <- rbinom(length(pa0_d), size = 1, prob = pa0_d)
    data_1[, "A"] <- (PredictA_d)

    # Compute Mu
    train_data$muldm <- predict(fit.y, newdata = as.matrix(train_data_1 %>% 
                                                           dplyr::select(all_of(selected_columns))))
    data$muldm <- predict(fit.y, newdata = as.matrix(data_1 %>% 
                                                             dplyr::select(all_of(selected_columns))))
    data$mu <- predict(fit.y, newdata=as.matrix(data %>% 
                                                  dplyr::select(all_of(selected_columns))))
    
    # Compute Nu 
    selected_columns1 <- if (sim == 0 | sim == 1) {
      c( "A", "X1", "X2", "X3", "C")
    } else {
      c( "A", "X1.1", "X2.1", "X3.1", "C")
    }
    
    X_final <- as.matrix(train_data %>% filter(R == 1) %>% dplyr::select(all_of(selected_columns1)))
    Y_final <- train_data %>% filter(R == 1) %>% pull(muldm)
    
    # Train the final outcome model on data with predicted M
    dtrain_final <- xgb.DMatrix(data = X_final, label = Y_final)
    fit_final <- xgb.cv(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6),
                        data = dtrain_final, 
                        nfold = 5,
                        nrounds = 100, 
                        early_stopping_rounds = 10,
                        verbose = 0)
    best_nrounds <- fit_final$best_iteration
    
    fit_final <- xgboost(data = dtrain_final, 
                         max_depth = 6, 
                         nrounds = best_nrounds, 
                         objective = "reg:squarederror", 
                         lambda = 1, 
                         alpha = 0.1, 
                         verbose = 0)
    
    X_A1C0 <- as.matrix(data_1 %>% dplyr::select(all_of(selected_columns1)))
    data$nu <- predict(fit_final, newdata= X_A1C0)
    
    #Weighting ::::::::
    # Weighting step for M
    train_data_r1 <- train_data %>% filter(R == 1)
    
    selected_columns2 <- if (sim == 0 | sim == 2) {
      c( "A", "X1", "X2", "X3", "Z", "C")
    } else {
      c( "A", "X1.1", "X2.1", "X3.1", "Z","C")
    }
    
    X_mediator <- as.matrix(train_data_r1 %>% dplyr::select(all_of(selected_columns2)))
    M_mediator <- train_data_r1$M
    
    # Train XGBoost model to predict M
    pos_weight1 <- sum(train_data_r1$M == 0) / sum(train_data_r1$M == 1)
    dtrain_m <- xgb.DMatrix(data = X_mediator, label = M_mediator)
    fit.m1 <- xgb.cv(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 6),
                     data = dtrain_m,
                     nfold = 5,
                     nrounds = 100,
                     early_stopping_rounds = 10,
                     verbose = 0)
    best_nrounds <- fit.m1$best_iteration

    fit.m1 <- xgboost(data = dtrain_m,
                      max_depth = 6,
                      nrounds = best_nrounds,
                      objective = "binary:logistic",
                      scale_pos_weight = pos_weight1,
                      verbose = 0)
    
    # Weighting step for A
    selected_columns3 <- if (sim == 0 | sim == 2 | sim == 3) {
      c(  "X1", "X2", "X3", "C")
    } else {
      c(  "X1.1", "X2.1", "X3.1", "C")
    }
    
    X_covariates <- as.matrix(train_data_r1 %>% dplyr::select(all_of(selected_columns3)))
    A_target <- train_data_r1$A
    
    dtrain_a <- xgb.DMatrix(data = X_covariates, label = A_target)
    pos_weight <- sum(train_data_r1$A == 0) / sum(train_data_r1$A == 1)
    fit.a1 <- xgb.cv(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 6),
                     data = dtrain_a,
                     nfold = 5,
                     nrounds = 100,
                     early_stopping_rounds = 10,
                     verbose = 0)
    best_nrounds <- fit.a1$best_iteration
    fit.a1 <- xgboost(data = dtrain_a,
                      max_depth = 6,
                      nrounds = best_nrounds,
                      scale_pos_weight = pos_weight,
                      objective = "binary:logistic",
                      verbose = 0)
    
    # Predict probabilities for M and X for the entire dataset (even outside R == 1)
     X_full_M <- as.matrix(data %>% dplyr::select(all_of(selected_columns2)))

     X_full_A <- as.matrix(data %>% dplyr::select(all_of(selected_columns3)))
    
    # Predict probabilities for M using XGBoost
    pm1 <- predict(fit.m1, newdata = X_full_M)
    
    # Predict probabilities for X using XGBoost
    pa1 <- predict(fit.a1, newdata = X_full_A)
    
    # Calculate weights for WM 
    WM <- rep(NA, nrow(data))
    ind11 <- data$R == 1 & data$M == 1
    ind10 <- data$R == 1 & data$M == 0
    WM[ind11] <- pm0[ind11] / pm1[ind11]
    WM[ind10] <- (1 - pm0[ind10]) / (1 - pm1[ind10])
    WM[data$R == 0] <- 0
    data$WM <- WM
    
    # Calculate weights for WX
    WA <- rep(NA, nrow(data))
    ind11a <- data$R == 1 & data$A == 1
    ind10a <- data$R == 1 & data$A == 0
    WA[ind11a] <- pa0[ind11a] / pa1[ind11a]
    WA[ind10a] <- (1-pa0[ind10a]) / (1 - pa1[ind10a])
    WA[data$R == 0] <- 0
    data$WA <- WA
    
    main_list[[k]] <- data[folds[[k]],]
  }
  
  main <- reduce(main_list, bind_rows) 
  # Calculate disparity estimates 
  wy0 <- lm(Y ~ C, data = subset(main, R == 0))
  wy1 <- lm(Y ~ C, data = subset(main, R == 1))
  
  #:::::::::: Imputation 
  b <- lm(nu ~ C, data = subset(main, R==1))
  delta_imp <- wy1$coef[1] - b$coef[1]
  zeta_imp <- b$coef[1] - wy0$coef[1] 
  
  #:::::::::: Weighting
  wmu1xdm <- lm(Y ~ C, weights = WA * WM, data = subset(main, R == 1))
  delta_wgt <- wy1$coef[1] - wmu1xdm$coef[1]
  zeta_wgt <- wmu1xdm$coef[1] - wy0$coef[1]
  
  #::::::::::::: Imputation and Weighting
  wa_mu <- lm(muldm ~ C, weight = WA, data = subset(main, R == 1))
  delta_iw = wy1$coef[1] - coef(wa_mu)[1]
  zeta_iw = coef(wa_mu)[1] -wy0$coef[1] 
  
  #::::::::::::: Triply Robust 
  wa_nu <- lm(nu ~  C , weights=WA, data = subset(main, R == 1))
  wawm_mu <- lm(mu ~ C, weight = WM*WA, data = subset(main, R == 1)) 
  
  delta_tr = wy1$coef[1] -(b$coef[1]+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))
  zeta_tr = (b$coef[1]+ (coef(wa_mu)[1]- coef(wa_nu)[1]) +(wmu1xdm$coef[1]- coef(wawm_mu)[1]))-wy0$coef[1] 
  
  return(c(delta_imp, zeta_imp, delta_wgt, zeta_wgt, delta_iw, zeta_iw, delta_tr, zeta_tr))
}


