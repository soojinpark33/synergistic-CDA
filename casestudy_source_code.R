# clear
    rm(list = ls())
# set system language
    Sys.setenv(LANGUAGE = "en")

# Call in packages
    library(readxl)
    library(dplyr)
    library(boot)
    library(parallel)
    library(dplyr)
    library(tidyverse)
    library(rlang)
    library(xgboost)
    library(causal.decomp)
    library(purrr)
    library(rbw)
    library(caret)
    

#::::::::::::::::::::::::::
# Function 1: GLM analysis
#::::::::::::::::::::::::::    

    calculate_disparities_glm <- function(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data, comparison, reference, cluster_id = "SCH_ID", n_bootstrap = 1000, seed=123) {
      set.seed = seed
      disparity_calculation_glm <- function(data_boot) {
        
        # Subset data for the reference and comparison groups
        data_r0 <- data_boot %>% filter(!!sym(rvar) == reference)
        data_r1 <- data_boot %>% filter(!!sym(rvar) == comparison)
        
        # Handle empty covariates
        covariates <- if (is.null(covariates) || length(covariates) == 0) "1" else paste(covariates, collapse = " + ")
        
        # Fit initial models for Y
        wy0_form <- as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + ")))
        wy0 <- lm(wy0_form, data = data_r0)
        wy1 <- lm(wy0_form, data = data_r1)
        
        # Fit models for M and A
        fit.m0 <- glm(as.formula(paste0(mvar, " ~ ", paste(covariates, collapse = " + "))), data = data_r0, family = binomial(logit))
        fit.a0 <- glm(as.formula(paste0(avar, " ~ ", paste(covariates, collapse = " + "))), data = data_r0, family = binomial(logit))
        
        pm0 <- predict(fit.m0, data_boot, type = "response")
        pa0 <- predict(fit.a0, data_boot, type = "response")
        
        # Fit model for Y
        fit.y <- lm(y_form, data = data_boot)
        
        # Predict outcome after incorporating predicted values of M and A (mu)
        dat_1 <- data_boot
        dat_1[[mvar]] <- as.factor(rbinom(length(pm0), size = 1, prob = pm0)) 
        dat_1[[avar]] <- as.factor(rbinom(length(pa0), size = 1, prob = pa0)) 
        data_boot$muldm <- predict(fit.y, newdata = dat_1)
        
        # Fit final model
        fit_final <- lm(as.formula(paste0("muldm", " ~ ", paste(c(avar,cvars,xvars), collapse = " + "))), data = data_boot %>% filter(!!sym(rvar) == comparison))
        data_boot$nu <- predict(fit_final, newdata =dat_1)
        b <- lm(as.formula(paste0("nu", " ~ ", paste(cvars, collapse = " + "))), data=data_boot%>% filter(!!sym(rvar) == comparison))
        
        # Disparity Reduction Estimates
        delta_imp <- wy1$coef[1] - b$coef[1]
        zeta_imp <- b$coef[1] - wy0$coef[1]
        
        # Pure Weighting
        fit.m1 <- glm(m_form, data = data_r1, family = binomial(logit))
        fit.a1 <- glm(a_form, data = data_r1, family = binomial(logit))
        
        pm1 <- predict(fit.m1, data_boot, type = "response")
        pa1 <- predict(fit.a1, data_boot, type = "response")
        
        # Calculate weights
        WM <- ifelse(data_boot[[rvar]] == comparison & data_boot[[mvar]] == 1, pm0 / pm1,
                     ifelse(data_boot[[rvar]] == comparison & data_boot[[mvar]] == 0, (1 - pm0) / (1 - pm1), 0))

        lower_bound <- quantile(WM[WM > 0], 0.01, na.rm = TRUE)  # Bottom 1%
        upper_bound <- quantile(WM[WM > 0], 0.99, na.rm = TRUE)  # Top 99%
        
        WM[WM > upper_bound] <- upper_bound
        WM[WM < lower_bound & WM > 0] <- lower_bound
        data_boot$WM <- WM
 
        WA <- ifelse(data_boot[[rvar]] == comparison & data_boot[[avar]] == 1, pa0 / pa1,
                     ifelse(data_boot[[rvar]] == comparison & data_boot[[avar]] == 0, (1-pa0) / (1 - pa1), 0))
        lower_bounda <- quantile(WA[WA > 0], 0.01, na.rm = TRUE)  # Bottom 1%
        upper_bounda <- quantile(WA[WA > 0], 0.99, na.rm = TRUE)  # Top 99%
        
        WA[WA > upper_bound] <- upper_bounda
        WA[WA < lower_bound & WA > 0] <- lower_bounda
        data_boot$WA <- WA
        
        # Weighted Y Calculation
        wmu1xdm <- lm(as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + "))), weights = WA * WM, 
                      data = data_boot %>% filter(!!sym(rvar) == comparison))
        delta_wgt <- wy1$coef[1] - wmu1xdm$coef[1]
        zeta_wgt <- wmu1xdm$coef[1] - wy0$coef[1]
        
        # Imputation and Weighting
        wx_mu <- lm(as.formula(paste0("muldm", " ~ ", paste(cvars, collapse = " + "))), weights = WA, data = data_boot %>% filter(!!sym(rvar) == comparison))
        delta_iw <- wy1$coef[1] - coef(wx_mu)[1]
        zeta_iw <- coef(wx_mu)[1] - wy0$coef[1]
        
        # Triply Robust Estimation
        wx_nu <- lm(as.formula(paste0("nu", " ~ ", paste(cvars, collapse = " + "))), weights = WA, data = data_boot %>% filter(!!sym(rvar) == comparison))
        data_boot$mu <- predict(fit.y, newdata = data_boot, type = "response")
        wxwm_mu <- lm(as.formula(paste0("mu", " ~ ", paste(cvars, collapse = " + "))), weights = WM* WA, data = data_boot %>% filter(!!sym(rvar) == comparison))
        
        delta_tr <- wy1$coef[1] - (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1]))
        zeta_tr <- (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1])) - wy0$coef[1]
        
        # Return all deltas and zetas
        return(c(delta_tr, zeta_tr))
      }

      # Formulas for the models
          # A ~ X + C
          a_rhs <- c(cvars, xvars) %>%
            reduce(~ expr(!!.x + !!.y))
          a_form <- as.formula(expr(!!sym(avar) ~ !!a_rhs))
          
          # M ~ X + A + Z + C
          m_rhs <- c(cvars, xvars, sym(avar), zvars) %>%
            reduce(~ expr(!!.x + !!.y))
          m_form <- as.formula(expr(!!sym(mvar) ~ !!m_rhs))
          
          # Y ~ R + X + A + Z + M + C
          interaction_terms <- exprs(!!sym(rvar) * !!sym(mvar), !!sym(rvar) * !!sym(avar))
          y_rhs <- c(cvars, sym(rvar), xvars, sym(avar), zvars, sym(mvar), interaction_terms) %>%
            reduce(~ expr(!!.x + !!.y))
          y_form <- as.formula(expr(!!sym(yvar) ~ !!y_rhs))
      
      # Perform clustered bootstrap
      clusters <- unique(data[[cluster_id]])
      glm_results <- replicate(n_bootstrap, {
        sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
        boot_data <- data %>% filter(!!sym(cluster_id) %in% sampled_clusters)
        disparity_calculation_glm(boot_data)
      }, simplify = TRUE)
    
      # Calculate means of all delta and zeta estimates
      delta_glm <- mean(glm_results[1, ])
      zeta_glm <- mean(glm_results[2, ])

      # Calculate standard errors from the bootstrap results
      delta_glm_se <- sd(glm_results[1, ])
      zeta_glm_se <- sd(glm_results[2, ])

      # Return the results with standard errors for all deltas and zetas
      return(list(
        delta_glm = delta_glm, zeta_glm = zeta_glm, delta_glm_se = delta_glm_se, zeta_glm_se = zeta_glm_se
      ))
    }
    
#::::::::::::::::::::::::::
# Function 2: XGBoost
#::::::::::::::::::::::::::  
   
    calculate_disparities_xg <- function(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data, comparison, reference, cluster_id = "SCH_ID", n_bootstrap = 1000, seed=123) {
      set.seed = seed

      disparity_calculation_xg <- function(data_boot) {
        
        # Subset data for the reference and comparison groups
        data_r0 <- data_boot %>% filter(!!sym(rvar) == reference)
        data_r1 <- data_boot %>% filter(!!sym(rvar) == comparison)
        
        # Handle empty covariates
        covariates <- if (is.null(covariates) || length(covariates) == 0) "1" else paste(covariates, collapse = " + ")
        
        # Fit initial models for Y
        wy0_form <- as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + ")))
        wy0 <- lm(wy0_form, data = data_r0) 
        wy1 <- lm(wy0_form, data = data_r1) 
        
        # Fit models for M and A
        fit.m0 <- glm(as.formula(paste0(mvar, " ~ ", paste(covariates, collapse = " + "))), data = data_r0, family = binomial(logit))
        fit.a0 <- glm(as.formula(paste0(avar, " ~ ", paste(covariates, collapse = " + "))), data = data_r0, family = binomial(logit))
        pm0 <- predict(fit.m0, data_boot, type = "response")
        pa0 <- predict(fit.a0, data_boot, type = "response")
        
        # set up training data for XGBoost
        X <- model.matrix(y_form, data = data_boot)[,-1] 
        Y <- data_boot[[yvar]]
        dtrain_y <- xgb.DMatrix(data = X, label = Y)
        
        # Train XGBoost model for outcome prediction
        fit.y <- xgboost(data = dtrain_y, max_depth = 6, nrounds = 100, 
                         objective = "reg:squarederror", early_stopping_rounds = 10, verbose = 0)
        
        # Predict outcome after incorporating predicted values of M and A (mu)
        dat_1 <- data_boot
        dat_1[[mvar]] <- rbinom(length(pm0), size = 1, prob = pm0)
        dat_1[[avar]] <- rbinom(length(pa0), size = 1, prob = pa0)
        X_new <- model.matrix(y_form, data = dat_1)[,-1] 
        colnames(X_new)[colnames(X_new) == mvar] <- paste0(mvar, "1")  # Adjust mediator column name if needed
        colnames(X_new)[colnames(X_new) == avar] <- paste0(avar, "1")  # Adjust mediator column name if needed
        data_boot$muldm <- predict(fit.y, newdata = X_new)
        
        # Compute outcomes for the comparison group and train final outcome model 
        X_final <- model.matrix(y_form2, data = data_r1)[,-1] 
        Y_final <- data_boot %>% filter(!!sym(rvar) == comparison) %>% pull(muldm)
        dtrain_final <- xgb.DMatrix(data = X_final, label = Y_final)
        fit_final <- xgboost(data = dtrain_final, max_depth = 6, nrounds = 100, 
                             objective = "reg:squarederror", early_stopping_rounds = 10, verbose = 0)
        
        X_X1C0 <- model.matrix(y_form2, data = dat_1)[,-1] 
        colnames(X_X1C0)[colnames(X_X1C0) == avar] <- paste0(avar, "1")  
        data_boot$nu <- predict(fit_final, newdata = X_X1C0)
        b <- lm(as.formula(paste0("nu", " ~ ", paste(cvars, collapse = " + "))), data=data_boot%>% filter(!!sym(rvar) == comparison))
        
        # Disparity Reduction Estimates
        delta_imp <- wy1$coef[1] - b$coef[1]
        zeta_imp <- b$coef[1] - wy0$coef[1]
        
        # Pure Weighting
        X_mediator <- model.matrix(m_form, data = data_r1)[,-1] 
        M_mediator <- as.numeric(as.character(data_r1[[mvar]]))
        dtrain_m <- xgb.DMatrix(data = X_mediator, label = M_mediator)
        fit.m1 <- xgboost(data = dtrain_m, max_depth = 6, nrounds = 100, 
                          objective = "binary:logistic", early_stopping_rounds = 10, verbose = 0)
        
        X_covariates <- model.matrix(a_form, data = data_r1)[,-1]
        A_target <- as.numeric(as.character(data_r1[[avar]]))
        dtrain_a <- xgb.DMatrix(data = X_covariates, label = A_target)
        fit.a1 <- xgboost(data = dtrain_a, max_depth = 6, nrounds = 100, 
                          objective = "binary:logistic", early_stopping_rounds = 10, verbose = 0)
        
        # Predict M and A for the full dataset
        pm1 <- predict(fit.m1, newdata = model.matrix(m_form, data = data_boot)[,-1])
        pa1 <- predict(fit.a1, newdata = model.matrix(a_form, data = data_boot)[,-1])
        
        # Calculate Weights (WM and WX)
        WM <- rep(NA, nrow(data_boot))
        ind11 <- data_boot[[rvar]] == comparison & data_boot[[mvar]] == 1
        ind10 <- data_boot[[rvar]] == comparison & data_boot[[mvar]] == 0
        WM[ind11] <- pm0[ind11] / pm1[ind11]
        WM[ind10] <- (1 - pm0[ind10]) / (1 - pm1[ind10])
        WM[data_boot[[rvar]] != comparison] <- 0
        
        lower_bound <- quantile(WM[WM > 0], 0.01, na.rm = TRUE)  # Bottom 1%
        upper_bound <- quantile(WM[WM > 0], 0.99, na.rm = TRUE)  # Top 99%
        
        WM[WM > upper_bound] <- upper_bound
        WM[WM < lower_bound & WM > 0] <- lower_bound
        data_boot$WM <- WM
        
        WA <- rep(NA, nrow(data_boot))
        ind11a <- data_boot[[rvar]] == comparison & data_boot[[avar]] == 1
        ind10a <- data_boot[[rvar]] == comparison & data_boot[[avar]] == 0
        WA[ind11a] <- pa0[ind11a] / pa1[ind11a]
        WA[ind10a] <-  (1 - pa0[ind10a]) / (1 - pa1[ind10a])
        WA[data_boot[[rvar]] != comparison] <- 0
        
        lower_bounda <- quantile(WA[WA > 0], 0.01, na.rm = TRUE)  # Bottom 1%
        upper_bounda <- quantile(WA[WA > 0], 0.99, na.rm = TRUE)  # Top 99%
        
        WA[WA > upper_bound] <- upper_bounda
        WA[WA < lower_bound & WA > 0] <- lower_bounda
        data_boot$WA <- WA
        
        
        # Weighted Y Calculation
        wmu1xdm <- lm(as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + "))), weights = WA * WM, 
                      data = data_boot%>% filter(!!sym(rvar) == comparison))
        delta_wgt <- wy1$coef[1] - wmu1xdm$coef[1]
        zeta_wgt <- wmu1xdm$coef[1] - wy0$coef[1]
        
        # Imputation and Weighting
        wx_mu <- lm(as.formula(paste0("muldm" , " ~ ", paste(cvars, collapse = " + "))), weight = WA, data = data_boot%>% filter(!!sym(rvar) == comparison) )
        delta_iw <- wy1$coef[1] - coef(wx_mu)[1]
        zeta_iw <- coef(wx_mu)[1] - wy0$coef[1]
        
        # Triply Robust Estimation
        wx_nu <- lm(as.formula(paste0("nu" , " ~ ", paste(cvars, collapse = " + "))), weight = WA, data = data_boot%>% filter(!!sym(rvar) == comparison))
        data_boot$mu <- predict(fit.y, newdata = dtrain_y)
        wxwm_mu <- lm(as.formula(paste0("mu" , " ~ ", paste(cvars, collapse = " + "))), weight = WA, data = data_boot%>% filter(!!sym(rvar) == comparison) )
        
        delta_tr <- wy1$coef[1] - (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1]))
        zeta_tr <- (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1])) - wy0$coef[1]
        
        # Return all deltas and zetas
        return(c(delta_tr, zeta_tr))
      }

      # Formulas for the models
      # A ~ X + C
      a_rhs <- c(cvars, xvars) %>%
        reduce(~ expr(!!.x + !!.y))
      a_form <- as.formula(expr(!!sym(avar) ~ !!a_rhs))
      
      # M ~ X + A + Z + C
      m_rhs <- c(cvars, xvars, sym(avar), zvars) %>%
        reduce(~ expr(!!.x + !!.y))
      m_form <- as.formula(expr(!!sym(mvar) ~ !!m_rhs))
      
      # Y ~ R + X + A + Z + M + C
      y_rhs <- c(cvars, sym(rvar), xvars, sym(avar), zvars, sym(mvar)) %>%
        reduce(~ expr(!!.x + !!.y))
      y_form <- as.formula(expr(!!sym(yvar) ~ !!y_rhs))
      
      # Y ~ (R) + X + A + C
      y_rhs2 <- c(cvars, xvars, sym(avar)) %>%
        reduce(~ expr(!!.x + !!.y))
      y_form2 <- as.formula(expr(!!sym(yvar) ~ !!y_rhs2))
      
      
      # Perform clustered bootstrap
      clusters <- unique(data[[cluster_id]])
      
      xg_results <- replicate(n_bootstrap, {
        sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
        boot_data <- data %>% filter(!!sym(cluster_id) %in% sampled_clusters)
        disparity_calculation_xg(boot_data)
      }, simplify = TRUE)
      
      # Calculate means of all delta and zeta estimates
      delta_xg <- mean(xg_results[1, ])
      zeta_xg <- mean(xg_results[2, ])
      
      # Calculate standard errors from the bootstrap results
      delta_xg_se <- sd(xg_results[1, ])
      zeta_xg_se <- sd(xg_results[2, ])
      
      
      # Return the results with standard errors for all deltas and zetas
      return(list(
        delta_xg = delta_xg, zeta_xg = zeta_xg, delta_xg_se = delta_xg_se, zeta_xg_se = zeta_xg_se
      ))
    }
    
#:::::::::::::::::::::::::::::::::::::::
# Function 3: XGBoost with crossfitting
#:::::::::::::::::::::::::::::::::::::::   
     
    calculate_disparities_xgcf <- function(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data, comparison, reference, cluster_id = "SCH_ID", n_bootstrap = 1000, seed=123) {
      set.seed = seed
      disparity_calculation_xgcf <- function(boot_data, K = 5) {
        
        # # Subset data for the reference and comparison groups
        data_r0 <- boot_data %>% filter(!!sym(rvar) == reference)
        data_r1 <- boot_data %>% filter(!!sym(rvar) == comparison)
        
        # Handle empty covariates
        covariates <- if (is.null(covariates) || length(covariates) == 0) "1" else paste(covariates, collapse = " + ")
        
        # number of cross-fitting folds
        folds <- createFolds(boot_data[[yvar]], K) 
        main_list <- vector(mode = "list", K)
        
        for (k in 1:K) {
          
          # Split data
          train_data <- boot_data[-folds[[k]], ]
          train_r0 <- train_data %>% filter(!!sym(rvar) == reference)
          train_r1 <- train_data %>% filter(!!sym(rvar) == comparison)
          
          # Fit models for M and A
          fit.m0 <- glm(as.formula(paste0(mvar, " ~ ", paste(covariates, collapse = " + "))), data = train_r0, family = binomial(logit))
          fit.a0 <- glm(as.formula(paste0(avar, " ~ ", paste(covariates, collapse = " + "))), data = train_r0, family = binomial(logit))
          pm0 <- predict(fit.m0, train_data, type = "response")
          pa0 <- predict(fit.a0, train_data, type = "response")
          
          # Set up training data for Y
          X <- model.matrix(y_form, data = train_data)[,-1] 
          Y <- train_data[[yvar]]
          dtrain_y <- xgb.DMatrix(data = X, label = Y)
          
          # Train XGBoost model for outcome prediction
          fit.y <- xgboost(data = dtrain_y, max_depth = 6, nrounds = 100, 
                           objective = "reg:squarederror", early_stopping_rounds = 10, verbose = 0)
          
          # Compute mu after incorporating predicted values of M & A using training data
          dat_1 <- train_data
          dat_1[[mvar]] <- rbinom(length(pm0), size = 1, prob = pm0)
          dat_1[[avar]] <- rbinom(length(pa0), size = 1, prob = pa0)
          X_new <- model.matrix(y_form, data = dat_1)[,-1] 
          colnames(X_new)[colnames(X_new) == mvar] <- paste0(mvar, "1")  
          colnames(X_new)[colnames(X_new) == avar] <- paste0(avar, "1")  
          train_data$muldm <- predict(fit.y, newdata = X_new)
          
          # Prepare data for imputation-then-weighting estimator
          fit.m0_d <- glm(m_form, data = data_r0, family = binomial(logit))
          fit.a0_d <- glm(m_form, data = data_r0, family = binomial(logit))
          pm0_d <- predict(fit.m0_d, boot_data, type = "response")
          pa0_d <- predict(fit.a0_d, boot_data, type = "response")
          
          # Compute mu after incorporating predicted values of M & A using validation data
          dat_1d <- boot_data
          dat_1d[[mvar]] <- rbinom(length(pm0_d), size = 1, prob = pm0_d)
          dat_1d[[avar]] <- rbinom(length(pa0_d), size = 1, prob = pa0_d)
          X_new_d <- model.matrix(y_form, data = dat_1d)[,-1] 
          colnames(X_new_d)[colnames(X_new_d) == mvar] <- paste0(mvar, "1")  
          colnames(X_new_d)[colnames(X_new_d) == avar] <- paste0(avar, "1")  
          boot_data$muldm <- predict(fit.y, newdata = X_new_d)
          boot_data$mu <- predict(fit.y, newdata = model.matrix(y_form, data = boot_data)[,-1])
          
          # Fit model for nu 
          X_final <- model.matrix(y_form2, data = train_r1)[,-1] 
          Y_final <- train_data %>% filter(!!sym(rvar) == comparison) %>% pull(muldm)
          dtrain_final <- xgb.DMatrix(data = X_final, label = Y_final)
          fit_final <- xgboost(data = dtrain_final, max_depth = 6, nrounds = 100, 
                               objective = "reg:squarederror", early_stopping_rounds = 10, verbose = 0)
          
          X_X1C0 <- model.matrix(y_form2, data = dat_1d)[,-1] 
          colnames(X_X1C0)[colnames(X_X1C0) == avar] <- paste0(avar, "1")
          boot_data$nu <- predict(fit_final, newdata = X_X1C0)
          
          # Fit models for M and A to compute Weights
          M_covariates <- model.matrix(m_form, data = train_r1)[,-1] 
          M_target <- as.numeric(as.character(train_r1[[mvar]]))
          dtrain_m <- xgb.DMatrix(data = M_covariates, label = M_target)
          pos_weight1 <- sum(train_r1[[mvar]] == 0) / sum(train_r1[[mvar]] == 1)
          fit.m1 <- xgboost(data = dtrain_m, max_depth = 6, nrounds = 100, scale_pos_weight = pos_weight1,
                            objective = "binary:logistic", early_stopping_rounds = 10, verbose = 0)
          
          A_covariates <- model.matrix(a_form, data = train_r1)[,-1]
          A_target <- as.numeric(as.character(train_r1[[avar]]))
          dtrain_a <- xgb.DMatrix(data = A_covariates, label = A_target)
          pos_weight <- sum(train_r1[[avar]] == 0) / sum(train_r1[[avar]] == 1)
          fit.a1 <- xgboost(data = dtrain_a, max_depth = 6, nrounds = 100, scale_pos_weight = pos_weight,
                            objective = "binary:logistic", early_stopping_rounds = 10, verbose = 0)
          
          # Predict M and A for the full dataset
          pm1 <- predict(fit.m1, newdata = model.matrix(m_form, data = boot_data)[,-1])
          pa1 <- predict(fit.a1, newdata = model.matrix(a_form, data = boot_data)[,-1])
          
          # Calculate Weights (WM and WA)
          WM <- rep(NA, nrow(boot_data))
          ind11 <- boot_data[[rvar]] == comparison & boot_data[[mvar]] == 1
          ind10 <- boot_data[[rvar]] == comparison & boot_data[[mvar]] == 0
          WM[ind11] <- pm0[ind11] / pm1[ind11]
          WM[ind10] <- (1 - pm0[ind10]) / (1 - pm1[ind10])
          WM[boot_data[[rvar]] != comparison] <- 0

          lower_bound <- quantile(WM[WM > 0], 0.01, na.rm = TRUE)  # Bottom 1%
          upper_bound <- quantile(WM[WM > 0], 0.99, na.rm = TRUE)  # Top 99%
          
          WM[WM > upper_bound] <- upper_bound
          WM[WM < lower_bound & WM > 0] <- lower_bound
          boot_data$WM <- WM
          
          WA <- rep(NA, nrow(boot_data))
          ind11a <- boot_data[[rvar]] == comparison & boot_data[[avar]] == 1
          ind10a <- boot_data[[rvar]] == comparison & boot_data[[avar]] == 0
          WA[ind11a] <- pa0[ind11a] / pa1[ind11a]
          WA[ind10a] <- (1 - pa0[ind10a]) / (1 - pa1[ind10a])
          WA[boot_data[[rvar]] != comparison] <- 0
          
          lower_bounda <- quantile(WA[WA > 0], 0.01, na.rm = TRUE)  # Bottom 1%
          upper_bounda <- quantile(WA[WA > 0], 0.99, na.rm = TRUE)  # Top 99%
          
          WA[WA > upper_bound] <- upper_bounda
          WA[WA < lower_bound & WA > 0] <- lower_bounda
          boot_data$WA <- WA
          
          main_list[[k]] <- boot_data[folds[[k]],]
        }
        
        main <- reduce(main_list, bind_rows) 
        
        # Calculate disparity estimates 
        wy0_form <- as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + ")))
        main_r0 <- main %>% filter(!!sym(rvar) == reference)
        main_r1 <- main %>% filter(!!sym(rvar) == comparison)
        wy0 <- lm(wy0_form, data = main_r0)
        wy1 <- lm(wy0_form, data = main_r1)
        
        #:::::::::: Imputation 
        b <- lm(as.formula(paste0("nu", " ~ ", paste(cvars, collapse = " + "))), main_r1)
        delta_imp <- wy1$coef[1] - b$coef[1]
        zeta_imp <- b$coef[1] - wy0$coef[1] 
        
        #:::::::::: Weighting 
        wmu1xdm <- lm(as.formula(paste0(yvar, " ~ ", paste(cvars, collapse = " + "))), weights = WA * WM, 
                      data = main%>% filter(!!sym(rvar) == comparison))
        delta_wgt <- wy1$coef[1] - wmu1xdm$coef[1]
        zeta_wgt <- wmu1xdm$coef[1] - wy0$coef[1]
        
        # Imputation and Weighting
        wx_mu <- lm(as.formula(paste0("muldm", " ~ ", paste(cvars, collapse = " + "))), weight = WA, 
                    data = main%>% filter(!!sym(rvar) == comparison))
        delta_iw <- wy1$coef[1] - coef(wx_mu)[1]
        zeta_iw <- coef(wx_mu)[1] - wy0$coef[1]
        
        # Triply Robust Estimation
        wx_nu <- lm(paste0("nu", " ~ ", paste(cvars, collapse = " + ")), weights = WA, data = main%>% filter(!!sym(rvar) == comparison))
        wxwm_mu <- lm(paste0("mu", " ~ ", paste(cvars, collapse = " + ")), weight = WM * WA, data = main%>% filter(!!sym(rvar) == comparison))
        
        delta_tr <- wy1$coef[1] - (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1]))
        zeta_tr <- (b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1])) - wy0$coef[1]
        return(c(delta_tr , zeta_tr ))
        
      } 
      # Formulas for the models
      # A ~ X + C
      a_rhs <- c(cvars, xvars) %>%
        reduce(~ expr(!!.x + !!.y))
      a_form <- as.formula(expr(!!sym(avar) ~ !!a_rhs))
      
      # M ~ X + A + Z + C
      m_rhs <- c(cvars, xvars, sym(avar), zvars) %>%
        reduce(~ expr(!!.x + !!.y))
      m_form <- as.formula(expr(!!sym(mvar) ~ !!m_rhs))
      
      # Y ~ R + X + A + Z + M + C
      y_rhs <- c(cvars, sym(rvar), xvars, sym(avar), zvars, sym(mvar)) %>%
        reduce(~ expr(!!.x + !!.y))
      y_form <- as.formula(expr(!!sym(yvar) ~ !!y_rhs))
      
      # Y ~ (R) + X + A + C
      y_rhs2 <- c(cvars, xvars, sym(avar)) %>%
        reduce(~ expr(!!.x + !!.y))
      y_form2 <- as.formula(expr(!!sym(yvar) ~ !!y_rhs2))
      
      
      # Perform clustered bootstrap
      clusters <- unique(data[[cluster_id]])
      
      xgcf_results <- replicate(n_bootstrap, {
        sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
        boot_data1 <- data %>% filter(!!sym(cluster_id) %in% sampled_clusters)
        disparity_calculation_xgcf(boot_data1, K=5)
      }, simplify = TRUE)
      
      
      # Calculate means of all delta and zeta estimates
      delta_xgcf <- mean(xgcf_results[1, ])
      zeta_xgcf <- mean(xgcf_results[2, ])
      
      # Calculate standard errors from the bootstrap results
      delta_xgcf_se <- sd(xgcf_results[1, ])
      zeta_xgcf_se <- sd(xgcf_results[2, ])
      
      
      # Return the results with standard errors for all deltas and zetas
      return(list(
        delta_xgcf = delta_xgcf, zeta_xgcf = zeta_xgcf, delta_xgcf_se = delta_xgcf_se, zeta_xgcf_se = zeta_xgcf_se
      ))
    }
    