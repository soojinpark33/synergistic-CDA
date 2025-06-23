
# R system language in English 
    Sys.setenv(LANG = "en")

# Set directory
    #getwd()
    #setwd("Z:/")

# call in packages and source file
    library(parallel)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
    library(grid)
    
    source("source_code_park_cv.R")


# simulation settings 
    n.sim <- 1 # replication - change to 1000    
    n.cores <- 1 # # of cores - change to 8
    # n.sample <- 1000 # sample size
      
      ### simulation 
    
    run_simulation <- function(n.sample, n.sim = n.sim, n.cores = n.cores) {
      
      combined.all <- function(index, n.sample) {
       
        # sample from popdata - sample #: 1000
        idx <- sample(1:N, n.sample, replace = FALSE)
        data <-popdata[idx, ]
        
  
        sim.results.glm <- mclapply(0:4, function(sim) glm.func(sim, data = data), mc.cores = n.cores)
        sim.results.xg <- mclapply(0:4, function(sim) xg.func(sim, data = data), mc.cores = n.cores) 
        sim.results.xg.cf <- mclapply(0:4, function(sim) xg.func.cf(sim, data = data), mc.cores = n.cores)
        
        data.glm <- do.call("rbind", sim.results.glm)  
        data.xg <- do.call("rbind", sim.results.xg)  
        data.xg.cf <- do.call("rbind", sim.results.xg.cf)  
        
        colnames(data.glm) <- c("delta_imp", "zeta_imp", "delta_wgt", "zeta_wgt", "delta_iw", "zeta_iw", "delta_tr", "zeta_tr")
        colnames(data.xg) <- c("delta_imp", "zeta_imp", "delta_wgt", "zeta_wgt", "delta_iw", "zeta_iw", "delta_tr", "zeta_tr")
        colnames(data.xg.cf) <- c("delta_imp", "zeta_imp", "delta_wgt", "zeta_wgt", "delta_iw", "zeta_iw", "delta_tr", "zeta_tr")
        
        return(list(data_glm = data.glm, data_xg = data.xg, data_xg.cf = data.xg.cf))
        #return(list(data_glm = data.glm))
      }
      
      # Run the simulation
      system.time(results.all <- mclapply(1:n.sim, combined.all, n.sample, mc.cores = n.cores))
      
      # Unpack the results into separate lists for each dataset type
      data.glm <- lapply(results.all, `[[`, "data_glm")
      data.xg <- lapply(results.all, `[[`, "data_xg")
      data.xg.cf <- lapply(results.all, `[[`, "data_xg.cf")
      
      # reshape the result to be grouped by simulation
      
      reshape.by.method <- function(data.list, n_sim) {
        reshaped.data <- vector("list", length = 5)  # Create a list for 5 methods (0 to 4)
        
        for (i in 0:4) {
          reshaped.data[[i + 1]] <- do.call(rbind, lapply(data.list, function(x) x[i + 1, ]))
        }
        
        names(reshaped.data) <- paste0("data.sim.", 0:4)
        
        return(reshaped.data)
      }
      
      reshaped.data.glm <- reshape.by.method(data.glm, n.sim)
      reshaped.data.xg <- reshape.by.method(data.xg, n.sim) 
      reshaped.data.xg.cf <- reshape.by.method(data.xg.cf, n.sim)  
      
      # save the result from simulations
      filename <- commandArgs(trailingOnly = T)[1]
      write.csv(reshaped.data.glm, file = paste0("glm/glm_data_", n.sample,"_",filename,".csv"))
      write.csv(reshaped.data.xg, file = paste0("xg/xg_data_", n.sample,"_",filename,".csv"))
      write.csv(reshaped.data.xg.cf, file = paste0("xg_cf/xg_data_cf_", n.sample,"_",filename,".csv"))
    }
    
      run_simulation(n.sample = 500, n.sim = n.sim, n.cores = n.cores)
      run_simulation(n.sample = 1000, n.sim = n.sim, n.cores = n.cores)
      run_simulation(n.sample = 2000, n.sim = n.sim, n.cores = n.cores)
      
      save.image("simulation_study.RData")
      