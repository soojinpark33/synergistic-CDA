# Call in packages
    library(xtable)
    library(gt) 
    
    source("casestudy_source_code.R")

#::::::::::::::::::::::::::
# Data preprocessing
#::::::::::::::::::::::::::

data2 <- read_excel("Z:/data2.xlsx")
#View(data2)

    # School ID
        names(data2)[names(data2) == "Cam"] <- "SCH_ID"
    
    # Baseline covariates
        data2$races <- as.factor(data2$races)
        data2$races <- relevel(data2$races, ref = "White")
        data2$native <- as.factor(data2$native)
        data2$female <- scale(as.numeric(data2$female), center=TRUE)#center baseline cov
        data2$native <- scale(as.numeric(data2$native), center=TRUE)#center baseline cov
        
    # Groups
        quantiles <- quantile(data2$SES_9, probs = c(0.2, 0.8), na.rm = TRUE)
        data2$SES_Q5 <- with(data2, ifelse(SES_9 <= quantiles[1], "low SES",
                                           ifelse(SES_9 >= quantiles[2], "high SES", "middle SES")))
        
    # treatment
        data2$TREAT1 <- (1- data2$TREAT1) 
        data2$TREAT1 <- as.factor(data2$TREAT1)
    
    #Aggregate SEC 
        data2 <- data2 %>%
          group_by(SCH_ID) %>%
          mutate(SEC = mean(SES_9, na.rm = TRUE),
          Atest = mean(test_9, na.rm = TRUE))%>%
          ungroup()
    
    #Generate dummy SEC for Q5
        # Recode SEC into a dummy variable where values in the top 20% are 1
            data2$SEC_H20 <- as.factor(ifelse(data2$SEC > quantile(data2$SEC, 0.8), 1, 0))
        # Recode SEC into a dummy variable where values in the top 40% are 1
            data2$SEC_H40 <- as.factor(ifelse(data2$SEC > quantile(data2$SEC, 0.6), 1, 0))
        # Recode Atest into a dummy variable where values in the top 40% are 1
            data2$Atest_H40 <- as.factor(ifelse(data2$Atest > quantile(data2$Atest, 0.6), 1, 0))
    
    # Make formulas
        cvars <- exprs(native, female)
        rvar <- "races"
        xvars <- exprs(SES_9,gpa_8)
        avar <- "Atest_H40"
        zvars <- exprs(SEC_H40,mteff_t, identity_9, utility_9, interest_9, engage_9,
                       colexp_9, fricollege_9, control_9, locale_9, climate_9)
        mvar <- "TREAT1"
        yvar <- "test_11"
        
#::::::::::::::::::::::::::        
# Replicate Table 1
#::::::::::::::::::::::::::
#:::::::: Initial disparity
fit.lm0 <- lm(formula=test_11 ~ native + female + races, data = data2)
#:::::::: Disparity Reduction and Remaining
fit.m <- glm(TREAT1 ~ races+ native + female, data = data2, family = binomial(logit))
fit.y <- lm(test_11 ~ native + female + races + SES_9 + Atest_H40 +gpa_8  +SEC_H40+
              mteff_t + identity_9 + utility_9 + interest_9 + engage_9 +
              colexp_9 + fricollege_9 + control_9 + locale_9 + climate_9 +
              TREAT1 + TREAT1*races, data = data2)
smiRes <- smi(fit.m = fit.m, fit.y = fit.y, sims = 50, conf.level = .95,
              covariates = c("female","native"), treat = "races", seed = 227) #Null covariates are not allowed!!

fit.m1 <- glm(Atest_H40 ~ races+ native + female, data = data2, family = binomial(logit))
fit.y1 <- lm(test_11 ~native + female + races + SES_9 + gpa_8 + Atest_H40 +
               Atest_H40*races, data = data2)
smiRes1 <- smi(fit.m = fit.m1, fit.y = fit.y1, sims = 500, conf.level = .95,
              covariates = c("female","native"), treat = "races", seed = 227) #Null covariates are not allowed!!

#::::::::::::::::::::::::::        
# Replicate Table 2
#::::::::::::::::::::::::::
#::::::::::::::::::::::::::
# Disparity Reduction/Remaining
#::::::::::::::::::::::::::

comparison_group <- "Hispanic"  # You can change this dynamically
reference_group <- "White"  # You can change this dynamically

    result.glm.h <- calculate_disparities_glm(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                              comparison=comparison_group, reference=reference_group,
                                    cluster_id = "SCH_ID", n_bootstrap = 500)

    result.xg.h <- calculate_disparities_xg(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                            comparison=comparison_group, reference=reference_group,
                                        cluster_id = "SCH_ID", n_bootstrap = 500)

    result.xgcf.h <- calculate_disparities_xgcf(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                                comparison=comparison_group, reference=reference_group,
                                                  cluster_id = "SCH_ID", n_bootstrap = 500)

comparison_group <- "Black"  # You can change this dynamically
reference_group <- "White"  # You can change this dynamically


result.glm.b <- calculate_disparities_glm(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                          comparison=comparison_group, reference=reference_group,
                                          cluster_id = "SCH_ID", n_bootstrap = 500)

result.xg.b <- calculate_disparities_xg(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                        comparison=comparison_group, reference=reference_group,
                                        cluster_id = "SCH_ID", n_bootstrap = 500)

result.xgcf.b <- calculate_disparities_xgcf(covariates=NULL, cvars, rvar, xvars, zvars, avar, mvar, yvar, data2,
                                            comparison=comparison_group, reference=reference_group,
                                            cluster_id = "SCH_ID", n_bootstrap = 500)

save.image(file = "all.RData")


