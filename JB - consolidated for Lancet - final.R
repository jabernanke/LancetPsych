# This R code accompanies the manuscript: 
#
# Structural brain measures among children with and without ADHD in the 
# Adolescent Brain and Cognitive Development Study cohort: 
# a cross-sectional US population-based study
#  
# It was developed by Joel Bernanke, Le Chang, and Jordan Dworkin.
#
# This code requires that you have access to the ABCD Study 2.0.1 release
# RDS file, which is publicly available.

# Libraries ---------------------------------------------------------------

library(mice)
library(lme4)
library(lmerTest)
library(emmeans)
library(EMAtools)
library(psych)
library(foreach)
library(doParallel)
library(abind)
library(Amelia)


# Function to find the mode of a factor -----------------------------------

mode_factor <- function(fac) {
  
  counts <- table(fac)
  modus <- names(counts)[max(counts) == counts][1]
  return(modus)
  
}


# Function to calculate effect sizes like ENIGMA --------------------------

calc_enig_e_s <- function(model, df) {
  
  t <- summary(model)$coefficients[2, 4]

  resd <- unlist(VarCorr(model))
  ser <- which(names(resd) == "serial")
  R <- as.numeric(resd[ser] / (sum(resd) + sigma(model)^2))
  
  no <- nrow(df)
  ni <- length(unique(df$serial))
  
  n1 <- table(df$adhd)[1]
  n2 <- table(df$adhd)[2]
  
  k <- nrow(summary(model)$coefficients)

  d = (t * (1 + (ni / no) * R) * sqrt(1 - R) * (n1 + n2)) / (sqrt(n1 * n2) * sqrt(no - k))
  return(as.numeric(d))
}


# Function to perform mixed effects modeling of brain measures ------------

# version = "full" for all covariates
# version = "reduced for just the Enigma covariates
model <- function(df, version = "reduced") {
  
  model_out <- data.frame(brain_var = character(),
                          beta = numeric(),
                          beta_se = numeric(),
                          beta_df = numeric(),
                          eff_size = numeric(),
                          eff_size_alt = numeric(),
                          lower_ci = numeric(),
                          upper_ci = numeric(),
                          p_value = numeric(),
                          adhd_mean = numeric(),
                          adhd_sd = numeric(),
                          no_adhd_mean = numeric(),
                          no_adhd_sd = numeric(),
                          diff_mean = numeric())
  
  # count the number of participants with and without adhd
  adhd_y <- sum(df$adhd == "yes")
  adhd_n <- sum(df$adhd == "no")

  for (i in 1:ncol(df)) {

    if (i == 1) j <- 0

    if(grepl(pattern = ".avg$|etiv", x = names(df)[i])) {

      j <- j + 1

      if (version == "full") {

        if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {

          m <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "race_eth", "educ", "income", 
                                  "married", "pds_p", "tbx_comp", "motion", "com_lumped",
                                  "(1|serial/family)"), 
                                collapse = " + ")))

        } else {

          m <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "race_eth", "educ", "income", 
                                  "married", "pds_p", "tbx_comp", "motion", "com_lumped", "scale(etiv)",
                                  "(1|serial/family)"), 
                                collapse = " + ")))
        }

      } else if (version == "reduced") {
        
        if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {

          m <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "(1|serial/family)"), 
                                collapse = " + ")))

        } else {

          m <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "scale(etiv)", 
                                  "(1|serial/family)"), 
                                collapse = " + ")))
        }
      }

      model <- lmer(formula = m, data = df, 
                    control = lmerControl(optimizer = "bobyqa"))  

      df_emmeans <- as.data.frame(emmeans(object = model, specs = "adhd"))
                  
      model_out[j, "brain_var"] <- names(df)[i]
      model_out[j, "beta"] <- coef(summary(model))[2, 1]
      model_out[j, "beta_se"] <- coef(summary(model))[2, 2]
      model_out[j, "beta_df"] <- coef(summary(model))[2, 3]
      model_out[j, "eff_size"] <- calc_enig_e_s(model =  model, df = df)
      model_out[j, "eff_size_alt"] <- lme.dscore(mod = model, data = df, type = "lme4")[1, 3]
      model_out[j, "lower_ci"] <- cohen.d.ci(d = model_out[j, "eff_size"], n1 = adhd_n, n2 = adhd_y)[1, 1]
      model_out[j, "upper_ci"] <- cohen.d.ci(d = model_out[j, "eff_size"], n1 = adhd_n, n2 = adhd_y)[1, 3]
      model_out[j, "p_value"] <- coef(summary(model))[2, 5]
      model_out[j, "adhd_mean"] <- df_emmeans[2, 2]
      model_out[j, "adhd_sd"] <- sd(resid(model)[df$adhd == "yes"])
      model_out[j, "no_adhd_mean"] <- df_emmeans[1, 2]
      model_out[j, "no_adhd_sd"] <- sd(resid(model)[df$adhd == "no"])
      model_out[j, "diff_mean"] <- model_out[j, "adhd_mean"] - model_out[j, "no_adhd_mean"]
      
    }
  }

  return(model_out)
   
}


# Function to return p-value on the adhd:sex interaction term -------------

model_int_sex <- function(df) {
  
  model_out <- data.frame(brain_var = character(),
                          int_p_value = numeric())

  for (i in 1:ncol(df)) {

    if (i == 1) j <- 0

    if(grepl(pattern = ".avg$|etiv", x = names(df)[i])) {

      j <- j + 1

      if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {

        m <- as.formula(paste(paste(names(df)[i], "~"), 
                        paste(c("adhd*sex", "age", "race_eth", "educ", "income", 
                                "married", "pds_p", "tbx_comp", "motion", "com_lumped",
                                "(1|serial/family)"), 
                              collapse = " + ")))

      } else {

        m <- as.formula(paste(paste(names(df)[i], "~"), 
                        paste(c("adhd*sex", "age", "race_eth", "educ", "income", 
                                "married", "pds_p", "tbx_comp", "motion", "com_lumped", "scale(etiv)",
                                "(1|serial/family)"), 
                              collapse = " + ")))
      }

      model <- lmer(formula = m, data = df, 
                    control = lmerControl(optimizer = "bobyqa"))  

      model_out[j, "brain_var"] <- names(df)[i]
      model_out[j, "int_p_value"] <- coef(summary(model))[nrow(coef(summary(model))), 5]
      
    }
  }

  return(model_out)
  
}


# Function to return p-value on the adhd:com interaction term -------------

model_int_com <- function(df) {
  
  model_out <- data.frame(brain_var = character(),
                          int_p_value = numeric())

  for (i in 1:ncol(df)) {

    if (i == 1) j <- 0

    if(grepl(pattern = ".avg$|etiv", x = names(df)[i])) {

      j <- j + 1

      if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {

        m <- as.formula(paste(paste(names(df)[i], "~"), 
                        paste(c("adhd*com_lumped", "age", "race_eth", "educ", "income", 
                                "married", "pds_p", "tbx_comp", "motion", "sex",
                                "(1|serial/family)"), 
                              collapse = " + ")))

      } else {

        m <- as.formula(paste(paste(names(df)[i], "~"), 
                        paste(c("adhd*com_lumped", "age", "race_eth", "educ", "income", 
                                "married", "pds_p", "tbx_comp", "motion", "sex", "scale(etiv)",
                                "(1|serial/family)"), 
                              collapse = " + ")))
      }

      model <- lmer(formula = m, data = df, 
                    control = lmerControl(optimizer = "bobyqa"))  

      model_out[j, "brain_var"] <- names(df)[i]
      model_out[j, "int_p_value"] <- coef(summary(model))[nrow(coef(summary(model))), 5]
      
    }
  }

  return(model_out)
  
}


# Function to return LR test p-value (adhd:race) --------------------------

model_int_race <- function(df) {
  
  model_out <- data.frame(brain_var = character(),
                          int_p_value = numeric())

  for (i in 1:ncol(df)) {

    if (i == 1) j <- 0

    if(grepl(pattern = ".avg$|etiv", x = names(df)[i])) {

      j <- j + 1

      if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {

        m_1 <- as.formula(paste(paste(names(df)[i], "~"), 
                        paste(c("adhd*race_eth", "age", "sex", "educ", "income", 
                                "married", "pds_p", "tbx_comp", "motion", "com_lumped",
                                "(1|serial/family)"), 
                              collapse = " + ")))
        
        m_2 <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "race_eth", "educ", "income", 
                                  "married", "pds_p", "tbx_comp", "motion", "com_lumped",
                                  "(1|serial/family)"), 
                                collapse = " + ")))

      } else {

        m_1 <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd*race_eth", "age", "sex", "educ", "income", 
                                  "married", "pds_p", "tbx_comp", "motion", "com_lumped", "scale(etiv)",
                                  "(1|serial/family)"), 
                                collapse = " + ")))
          
        m_2 <- as.formula(paste(paste(names(df)[i], "~"), 
                          paste(c("adhd", "age", "sex", "race_eth", "educ", "income", 
                                  "married", "pds_p", "tbx_comp", "motion", "com_lumped", "scale(etiv)",
                                  "(1|serial/family)"), 
                                collapse = " + ")))
      }

      model_1 <- lmer(formula = m_1, data = df, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa"))  
      
      model_2 <- lmer(formula = m_2, data = df, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa"))
      
      temp <- anova(model_1, model_2)

      model_out[j, "brain_var"] <- names(df)[i]
      model_out[j, "int_p_value"] <- temp$'Pr(>Chisq)'[2]
      
    }
  }

  return(model_out)
  
}


# Function to calculate FDR p-values by area, thickness, and volume -------

make_fdr_group <- function(df, raw_p) {
  
  out <- vector(length = nrow(df))
  
  out[grepl(pattern = "vol_|etiv", x = df$brain_var)] <-
    p.adjust(p = df[, raw_p][grepl(pattern = "vol_|etiv", x = df$brain_var)], 
             method = "BH")

  out[grepl(pattern = "thick_", x = df$brain_var)] <-
    p.adjust(p = df[, raw_p][grepl(pattern = "thick_", x = df$brain_var)], 
             method = "BH")

  out[grepl(pattern = "area_", x = df$brain_var)] <-
    p.adjust(p = df[, raw_p][grepl(pattern = "area_", x = df$brain_var)], 
             method = "BH")

  return(out)  
  
}

# Function to simulate "true" adhd diagnoses ------------------------------

sim_diag <- function(df, pct_fp, pct_fn) {

  # create subjectkeys
  df$subjectkey <- c(1:nrow(df))
  
  # create obs(erved) and "true" adhd variables
  names(df)[names(df) == "adhd"] <- "adhd_obs" 
  df$adhd_true <- df$adhd_obs
  
  # re-classify some ADHD diagnoses as false positives
  obs_index <- df$subjectkey[df$adhd_obs == "yes"]
  obs_index_fp <- sort(sample(x = obs_index, size = length(obs_index) * pct_fp, replace = FALSE))
  
  # re-classify some ADHD diagnoses as false negatives
  obs_index <- df$subjectkey[df$adhd_obs == "no"]
  obs_index_fn <- sort(sample(x = obs_index, size = length(obs_index) * pct_fn, replace = FALSE))
  
  # replace the observed with the "true" cases
  df$adhd_true[obs_index_fp] <- "no"
  df$adhd_true[obs_index_fn] <- "yes"
  
  return(df)

}


# Function to make a list of models, one for each brain variable ----------

make_model_list <- function(df, type = "sim") {
  
  # if called on simulated data, "adhd_true" should be the outcome variable;
  # if called on the ABCD data, "adhd_obs" should be the outcome variable
  adhd_type <- ifelse(type == "sim", "adhd_true", "adhd_obs")
  
  for (i in 1:ncol(df)) {
   
    if (i == 1) {
      j <- 0
      out <- list()
    }
    
    if(grepl(pattern = ".avg$|etiv", x = names(df)[i])) {
      
      j <- j + 1
      
      if(grepl(pattern = "thick_|etiv", x = names(df)[i])) {
        
        m <- as.formula(paste(paste(names(df)[i], "~"), 
                              paste(c(adhd_type, "age", "sex", "(1|serial/family)"), 
                                    collapse = " + ")))
        
      } else {
        
        m <- as.formula(paste(paste(names(df)[i], "~"), 
                              paste(c(adhd_type, "age", "sex", "scale(etiv)", 
                                      "(1|serial/family)"), 
                                    collapse = " + ")))
      }
      
      model <- lmer(formula = m, data = df, 
                    control = lmerControl(optimizer = "bobyqa")) 
      
      out[[j]] <- model
      
    }
  }
  
  return(out)
  
}


# Function to calculate new coefficients ----------------------------------

new_coef <- function(df, model, enigma_eff_size) {
  
  coefs <- summary(model)$coefficients
  
  t <- coefs[2, 4]
  se <- coefs[2, 2]
  
  resd <- unlist(VarCorr(model))
  ser <- which(names(resd) == "serial")
  R <- as.numeric(resd[ser] / (sum(resd) + sigma(model)^2))
  
  no <- nrow(df)
  ni <- length(unique(df$serial))
  
  n1 <- table(df$adhd_obs)[1]
  n2 <- table(df$adhd_obs)[2]
  
  k <- nrow(summary(model)$coefficients)
  
  new_t <- enigma_eff_size * (sqrt(n1 * n2) * sqrt(no - k)) / ((1 + (ni/no) * R) * sqrt(1 - R) * (n1 + n2))
  new_coef <- new_t * se
  
  return(as.numeric(new_coef))
  
}


# Function to extract, change model parameters ----------------------------

get_mod_params <- function(df, model, enigma_eff_size) {
  
  coefs <- summary(model)$coefficients
  coefs <- coefs[, 1]
  coefs[2] <- new_coef(df, model, enigma_eff_size)
  
  rand_effs <- summary(model)$varcor
  fam_ser <- sqrt(as.numeric(rand_effs$`family:serial`))
  ser <- sqrt(as.numeric(rand_effs$`serial`))
  resid <- summary(model)$sigma
  rand <- c(fam_ser, ser, resid)
  names(rand) <- c("fs_rsd", "serial_rsd", "residual_sd")
  coefs <- t(c(coefs, rand))
  colnames(coefs)[grepl("adhd",colnames(coefs))] <- "adhdyes"
  
  return(as.data.frame(coefs))
}


# Function to simulate brain data -----------------------------------------

# e_s_fixed allows for a single effect size to be assigned for all brain
# variables; if NULL, Enigma effect sizes are used
sim_brain_data <- function(df, m_list, e_s_fixed = NULL) {
  
  for (i in 1:length(m_list)) {
    
    if (is.null(e_s_fixed)) {
      
      var <- paste(formula(m_list[[i]])[2])
      
      if (any(var %in% enigma_eff_size$brain_var)) {
        
        e_s <- enigma_eff_size[enigma_eff_size$brain_var == var, "eff_size"] 
      
      } else e_s <- NA
      
      # some variables might be NA in the look up table, some might not be in 
      # lookup table at all; default missing effect sizes to the value below
      if (is.na(e_s)) e_s <- -0.1
    
    }
    
    mod_params <- get_mod_params(df = df, model = m_list[[i]], enigma_eff_size = e_s)
    ser_res <- rnorm(n = length(unique(df$serial)), mean = 0, sd = mod_params$serial_rsd)
    fam_res <- rnorm(n = length(unique(df$family)), mean = 0, sd = mod_params$fs_rsd)
    
    if (any(grepl(pattern = "etiv", x = formula(m_list[[i]])[3]))) {
      
      sim_brain <- mod_params$`(Intercept)` + 
        (as.numeric(df$adhd_true) - 1) * mod_params$adhdyes + 
        df$age * mod_params$age + 
        (as.numeric(df$sex) - 1) * mod_params$sexM + 
        scale(df$etiv) * mod_params$`scale(etiv)` +
        ser_res[as.numeric(df$serial)] + 
        fam_res[as.numeric(df$family)] +
        rnorm(n = nrow(df), mean = 0, sd = mod_params$residual_sd) 
      
    } else {
      
      sim_brain <- mod_params$`(Intercept)` + 
        (as.numeric(df$adhd_true) - 1) * mod_params$adhdyes + 
        df$age * mod_params$age + 
        (as.numeric(df$sex) - 1) * mod_params$sexM + 
        ser_res[as.numeric(df$serial)] + 
        fam_res[as.numeric(df$family)] +
        rnorm(n = nrow(df), mean = 0, sd = mod_params$residual_sd) 
      
    }
    
    if (i == 1) {
      
      out <- sim_brain
      names_out <- paste(formula(m_list[[i]])[2])
      
    } else {
      
      out <- cbind(out, sim_brain)
      names_out <- c(names_out, paste(formula(m_list[[i]])[2]))
      
    }
    
  }
  
  out <- as.data.frame(out)
  names(out) <- names_out
  
  names(out)[names(out) == "etiv"] <- "etiv_sim"
  out$etiv_obs <- df$etiv
  
  out <- cbind(df[, c("subjectkey", "adhd_obs", "adhd_true", "age", "sex", "family", "serial")], out)
  return(out)
  
}


# Function to perform mixed effects modeling in simulated data ------------

sim_model <- function(df) {
  
  model_out <- data.frame(brain_var = character(),
                          eff_size_obs = numeric(),
                          p_value_obs = numeric(),
                          eff_size_sim = numeric(),
                          p_value_sim = numeric())
  
  for (i in 1:ncol(df)) {
    
    if (i == 1) j <- 0
    
    if(grepl(pattern = ".avg$|etiv_sim", x = names(df)[i])) {
      
      j <- j + 1
      
      if(grepl(pattern = "thick|etiv_sim", x = names(df)[i])) {
        
        m_obs <- as.formula(paste(paste(names(df)[i], "~"), 
                                  paste(c("adhd_obs", "age", "sex", "(1|serial/family)"), 
                                        collapse = " + ")))
        
        m_sim <- as.formula(paste(paste(names(df)[i], "~"), 
                                  paste(c("adhd_true", "age", "sex", "(1|serial/family)"), 
                                        collapse = " + ")))
        
      } else {
        
        m_obs <- as.formula(paste(paste(names(df)[i], "~"), 
                                  paste(c("adhd_obs", "age", "sex", "scale(etiv_obs)", 
                                          "(1|serial/family)"), 
                                        collapse = " + ")))
        
        m_sim <- as.formula(paste(paste(names(df)[i], "~"), 
                                  paste(c("adhd_true", "age", "sex", "scale(etiv_obs)", 
                                          "(1|serial/family)"), 
                                        collapse = " + ")))
      }
      
      lme_obs <- lmer(formula = m_obs, data = df, 
                      control = lmerControl(optimizer = "bobyqa"))  
      
      df_emmeans_obs <- as.data.frame(emmeans(object = lme_obs, specs = "adhd_obs"))   
      
      lme_sim <- lmer(formula = m_sim, data = df, 
                      control = lmerControl(optimizer = "bobyqa"))  
      
      df_emmeans_sim <- as.data.frame(emmeans(object = lme_sim, specs = "adhd_true"))
      
      model_out[j, "brain_var"] <- names(df)[i]
      model_out[j, "eff_size_obs"] <- calc_enig_e_s(model = lme_obs, df = df)
      model_out[j, "p_value_obs"] <- coef(summary(lme_obs))[2, 5]
      model_out[j, "eff_size_sim"] <- calc_enig_e_s(model = lme_sim, df = df)
      model_out[j, "p_value_sim"] <- coef(summary(lme_sim))[2, 5]
      
    }
  }
  
  return(model_out)
  
}


# Function to collate effect sizes from the simulation --------------------

make_df_eff_sizes <- function(out_sim_results) {
  
  for (i in 1:length(out_sim_results)) {
  
    if (i == 1) {
      
      df_eff_sizes <- as.data.frame(out_sim_results[[1]]$brain_var)
      names(df_eff_sizes) <- "brain_var"
      
      for (j in df_eff_sizes$brain_var) {
      
        if (j %in% enigma_eff_size$brain_var) {
        
          df_eff_sizes$enigma_e_s[df_eff_sizes$brain_var == j] <- enigma_eff_size$eff_size[enigma_eff_size$brain_var == j]
        
        }
      }
    }
    
    df_eff_sizes[, paste("sim_e_s_", i, sep = "")] <- round(out_sim_results[[i]]$eff_size_sim, digits = 3)
    df_eff_sizes[, paste("obs_e_s_", i, sep = "")] <- round(out_sim_results[[i]]$eff_size_obs, digits = 3)
    
  }
  
  df_eff_sizes$sim_min <- apply(X = df_eff_sizes[, grepl(pattern = "sim", x = names(df_eff_sizes))],
                                  MARGIN = 1, FUN = range)[1, ]
  df_eff_sizes$sim_max <- apply(X = df_eff_sizes[, grepl(pattern = "sim", x = names(df_eff_sizes))],
                                  MARGIN = 1, FUN = range)[2, ]
  df_eff_sizes$sim_avg <- apply(X = df_eff_sizes[, grepl(pattern = "sim", x = names(df_eff_sizes))],
                                MARGIN = 1, FUN = mean)
  
  df_eff_sizes$obs_min <- apply(X = df_eff_sizes[, grepl(pattern = "obs", x = names(df_eff_sizes))],
                                  MARGIN = 1, FUN = range)[1, ]
  df_eff_sizes$obs_max <- apply(X = df_eff_sizes[, grepl(pattern = "obs", x = names(df_eff_sizes))],
                                  MARGIN = 1, FUN = range)[2, ]
  df_eff_sizes$obs_avg <- apply(X = df_eff_sizes[, grepl(pattern = "obs", x = names(df_eff_sizes))],
                                MARGIN = 1, FUN = mean)
  

  
  return(df_eff_sizes)
  
}


# Function to collate p-values from the simulation ------------------------

make_df_p_vals <- function(out_sim_results) {
  
  for (i in 1:length(out_sim_results)) {
    
    if (i == 1) {
      
      df_p_vals <- as.data.frame(out_sim_results[[1]]$brain_var)
      names(df_p_vals) <- "brain_var"
      
      for (j in df_p_vals$brain_var) {
      
        if (j %in% enigma_eff_size$brain_var) {
        
          df_p_vals$enigma_e_s[df_p_vals$brain_var == j] <- enigma_eff_size$eff_size[enigma_eff_size$brain_var == j]
        
        }
      }
    }
   
    df_p_vals[, paste("raw_", i, sep = "")] <- out_sim_results[[i]]$p_value_obs
    df_p_vals[, paste("fdr_", i, sep = "")] <- p.adjust(out_sim_results[[i]]$p_value_obs, method = "BH")
    df_p_vals[, paste("fdr_group_", i, sep = "")] <- NA
    
    df_p_vals[grepl(pattern = "vol_|etiv", x = df_p_vals$brain_var), paste("fdr_group_", i, sep = "")] <-
      p.adjust(out_sim_results[[i]]$p_value_obs[grepl(pattern = "vol_|etiv", x = df_p_vals$brain_var)], method = "BH")
    
    df_p_vals[grepl(pattern = "thick_", x = df_p_vals$brain_var), paste("fdr_group_", i, sep = "")] <-
      p.adjust(out_sim_results[[i]]$p_value_obs[grepl(pattern = "thick_", x = df_p_vals$brain_var)], method = "BH")
   
    df_p_vals[grepl(pattern = "area_", x = df_p_vals$brain_var), paste("fdr_group_", i, sep = "")] <-
      p.adjust(out_sim_results[[i]]$p_value_obs[grepl(pattern = "area_", x = df_p_vals$brain_var)], method = "BH")
    
  }

  df_p_vals$fdr_group_min <- apply(X = df_p_vals[, grepl(pattern = "fdr_group", x = names(df_p_vals))],
                                   MARGIN = 1, FUN = range)[1, ]
  df_p_vals$fdr_group_max <- apply(X = df_p_vals[, grepl(pattern = "fdr_group", x = names(df_p_vals))],
                                   MARGIN = 1, FUN = range)[2, ]
  
  df_p_vals$group_always_sig <- df_p_vals$fdr_group_max < 0.05
  
  return(df_p_vals)

}


# Function to generate a distribution of significant results from  --------

make_df_distro_sig <- function(df_p_vals) { 

  df_distro_sig <- data.frame(fdr = apply(X = df_p_vals[, grepl(pattern = "fdr_[1-9]", x = names(df_p_vals))], 
                                           MARGIN = 2, FUN = function(x) { sum(as.numeric(x < 0.05)) }),
                               
                              fdr_group = apply(X = df_p_vals[, grepl(pattern = "fdr_group_[1-9]", x = names(df_p_vals))], 
                                                 MARGIN = 2, FUN = function(x) { sum(as.numeric(x < 0.05)) }))

  return(df_distro_sig)
  
}
  

# Store the Enigma effect sizes (children) in a dataframe -----------------

# no ventraldc in Enigma results; set to median of the Enigma effect sizes
# for subcortical volume
enigma_eff_size <- data.frame(
  
  brain_var = c("etiv",                                  
  "smri_thick_cort.desikan_bankssts.avg",                
  "smri_thick_cort.desikan_caudalanteriorcingulate.avg", 
  "smri_thick_cort.desikan_caudalmiddlefrontal.avg",     
  "smri_thick_cort.desikan_cuneus.avg",                  
  "smri_thick_cort.desikan_entorhinal.avg",              
  "smri_thick_cort.desikan_fusiform.avg",                
  "smri_thick_cort.desikan_inferiorparietal.avg",        
  "smri_thick_cort.desikan_inferiortemporal.avg",        
  "smri_thick_cort.desikan_isthmuscingulate.avg",        
  "smri_thick_cort.desikan_lateraloccipital.avg",        
  "smri_thick_cort.desikan_lateralorbitofrontal.avg",    
  "smri_thick_cort.desikan_lingual.avg",                 
  "smri_thick_cort.desikan_medialorbitofrontal.avg",     
  "smri_thick_cort.desikan_middletemporal.avg",          
  "smri_thick_cort.desikan_parahippocampal.avg",         
  "smri_thick_cort.desikan_paracentral.avg",             
  "smri_thick_cort.desikan_parsopercularis.avg",         
  "smri_thick_cort.desikan_parsorbitalis.avg",           
  "smri_thick_cort.desikan_parstriangularis.avg",        
  "smri_thick_cort.desikan_pericalcarine.avg",           
  "smri_thick_cort.desikan_postcentral.avg",             
  "smri_thick_cort.desikan_posteriorcingulate.avg",      
  "smri_thick_cort.desikan_precentral.avg",              
  "smri_thick_cort.desikan_precuneus.avg",               
  "smri_thick_cort.desikan_rostralanteriorcingulate.avg",
  "smri_thick_cort.desikan_rostralmiddlefrontal.avg",    
  "smri_thick_cort.desikan_superiorfrontal.avg",         
  "smri_thick_cort.desikan_superiorparietal.avg",        
  "smri_thick_cort.desikan_superiortemporal.avg",        
  "smri_thick_cort.desikan_supramarginal.avg",           
  "smri_thick_cort.desikan_frontalpole.avg",             
  "smri_thick_cort.desikan_temporalpole.avg",            
  "smri_thick_cort.desikan_transversetemporal.avg",      
  "smri_thick_cort.desikan_insula.avg",                  
  "smri_thick_cort.desikan_mean.avg",                    
  "smri_area_cort.desikan_bankssts.avg",                 
  "smri_area_cort.desikan_caudalanteriorcingulate.avg",  
  "smri_area_cort.desikan_caudalmiddlefrontal.avg",      
  "smri_area_cort.desikan_cuneus.avg",                   
  "smri_area_cort.desikan_entorhinal.avg",               
  "smri_area_cort.desikan_fusiform.avg",                 
  "smri_area_cort.desikan_inferiorparietal.avg",         
  "smri_area_cort.desikan_inferiortemporal.avg",         
  "smri_area_cort.desikan_isthmuscingulate.avg",         
  "smri_area_cort.desikan_lateraloccipital.avg",         
  "smri_area_cort.desikan_lateralorbitofrontal.avg",     
  "smri_area_cort.desikan_lingual.avg",                  
  "smri_area_cort.desikan_medialorbitofrontal.avg",      
  "smri_area_cort.desikan_middletemporal.avg",           
  "smri_area_cort.desikan_parahippocampal.avg",          
  "smri_area_cort.desikan_paracentral.avg",              
  "smri_area_cort.desikan_parsopercularis.avg",          
  "smri_area_cort.desikan_parsorbitalis.avg",            
  "smri_area_cort.desikan_parstriangularis.avg",         
  "smri_area_cort.desikan_pericalcarine.avg",            
  "smri_area_cort.desikan_postcentral.avg",              
  "smri_area_cort.desikan_posteriorcingulate.avg",       
  "smri_area_cort.desikan_precentral.avg",               
  "smri_area_cort.desikan_precuneus.avg",                
  "smri_area_cort.desikan_rostralanteriorcingulate.avg", 
  "smri_area_cort.desikan_rostralmiddlefrontal.avg",     
  "smri_area_cort.desikan_superiorfrontal.avg",          
  "smri_area_cort.desikan_superiorparietal.avg",         
  "smri_area_cort.desikan_superiortemporal.avg",         
  "smri_area_cort.desikan_supramarginal.avg",            
  "smri_area_cort.desikan_frontalpole.avg",              
  "smri_area_cort.desikan_temporalpole.avg",             
  "smri_area_cort.desikan_transversetemporal.avg",       
  "smri_area_cort.desikan_insula.avg",                   
  "smri_area_cort.desikan_total.avg",                    
  "smri_vol_subcort.aseg_thalamus.proper.avg",           
  "smri_vol_subcort.aseg_caudate.avg",                   
  "smri_vol_subcort.aseg_putamen.avg",                   
  "smri_vol_subcort.aseg_pallidum.avg",                  
  "smri_vol_subcort.aseg_hippocampus.avg",               
  "smri_vol_subcort.aseg_amygdala.avg",                  
  "smri_vol_subcort.aseg_accumbens.area.avg",            
  "smri_vol_subcort.aseg_ventraldc.avg"),
  
  eff_size = c(-0.14,
  -0.06,
  -0.02,
  -0.07,
  -0.02,
  -0.09,
  -0.17,
  -0.08,
  -0.08,
  0.02,
  -0.1,
  -0.05,
  -0.08,
  0.02,
  -0.07,
  -0.15,
  -0.09,
  -0.07,
  -0.05,
  0,
  -0.04,
  -0.08,
  0,
  -0.16,
  -0.1,
  0.06,
  -0.03,
  -0.01,
  -0.07,
  -0.05,
  -0.07,
  -0.02,
  -0.18,
  0.06,
  -0.09,
  -0.05,
  -0.1,
  -0.08,
  -0.15,
  -0.06,
  -0.05,
  -0.13,
  -0.12,
  -0.12,
  -0.13,
  -0.12,
  -0.17,
  -0.09,
  -0.16,
  -0.13,
  -0.04,
  -0.07,
  -0.09,
  -0.07,
  -0.1,
  -0.04,
  -0.1,
  -0.16,
  -0.1,
  -0.12,
  -0.16,
  -0.13,
  -0.19,
  -0.12,
  -0.15,
  -0.13,
  -0.05,
  -0.1,
  -0.07,
  -0.12,
  -0.21,
  0.01,
  -0.13,
  -0.18,
  -0.01,
  -0.12,
  -0.18,
  -0.19,
  -0.13))


# Load the ABCD 2.0.1 pre-made Rds file -----------------------------------

df <- readRDS(file = "../abcd201/nda2.0.1.Rds")


# Subset to baseline observations -----------------------------------------

df <- subset(x = df, subset = (eventname == "baseline_year_1_arm_1"))


# Drop participants with missing or poor data -----------------------------

# missing all relevant imaging data
df <- df[apply(X = df[, grep("smri_area.*desikan|smri_thick.*desikan|smri.*vol_subcort", names(df))], 
               MARGIN = 1, FUN = function(x) { !all(is.na(x))}), ]

# unacceptable T1 or FreeSurfer QC score
df <- subset(x = df, subset = (df$iqc_t1_1_qc_score == "Pass" & df$fsqc_qc == "accept"))

# missing ADHD current diagnosis
df <- subset(x = df, subset = (!is.na(df$ksads_14_853_p)))

# missing site, scanner serial number, or family information
df <- subset(x = df, subset = (!is.na(df$site_id_l) & 
                               !is.na(df$mri_info_device.serial.number) &
                               !is.na(df$rel_family_id)))


# Fix the motion QC score -------------------------------------------------

# there should only be integer values but some decimal values are present;
# we round up to fix
df$fsqc_qu_motion <- ceiling(x = df$fsqc_qu_motion)


# Merge and fix the pubertal development scores ---------------------------

# pubertal development scores are stored in separate variables for boys
# and girls; merge them by adding them, or mark them NA if both are NA

pds_vars_p <- c("pubertdev_ss_male_category_p", "pubertdev_ss_female_category_p")
df$pds_p <- apply(X = df[, pds_vars_p], MARGIN = 1, FUN = function(x) {
  ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE)) })
rm(pds_vars_p)

# some pubertal development scores are listed as 5 (range 1-4); mark these
# as NA
df$pds_p[df$pds_p == 5] <- NA


# Generate a lumped anxiety/depression variable ---------------------------

# diagnoses must be current or in partial remission; not past
# marked as 1 if any disorder is present; otherwise marked as 0
# variable is considered missing only if all relevant diagnoses are missing

# anxiety disorders
# note: mutism variables are not actually present in the dataset
# PD, agoraphobia, separation, social, phobia, GAD (all present)
# "ksads_9_867_p" specific phobia
com_anx_vars <- c("ksads_5_857_p", "ksads_6_859_p", "ksads_7_861_p",
                  "ksads_8_863_p", "ksads_9_867_p", "ksads_10_869_p")
df$com_anx <- apply(X = df[, com_anx_vars], MARGIN = 1, FUN = function(x) {
  if(all(is.na(x))) return(NA)
  else if(any(x == 1, na.rm = TRUE)) return(1)
  else return(0)
})
df$com_anx <- factor(df$com_anx, levels = c(0, 1), labels = c("no", "yes"))
rm(com_anx_vars)

# OCD
# present
df$com_ocd <- factor(df$ksads_11_917_p, levels = c(0, 1), labels = c("no", "yes"))

# depressive disorders
# PDD present, PDD partial remission, MDD present, MDD partial remission
com_dep_vars <- c("ksads_1_843_p", "ksads_1_844_p", "ksads_1_840_p", "ksads_1_841_p")
df$com_dep <- apply(X = df[, com_dep_vars], MARGIN = 1, FUN = function(x) {
  if(all(is.na(x))) return(NA)
  else if(any(x == 1, na.rm = TRUE)) return(1)
  else return(0)
})
df$com_dep <- factor(df$com_dep, levels = c(0, 1), labels = c("no", "yes"))
rm(com_dep_vars)

# make a single co-morbidity variable with anxiety, OCD, and depression
com_vars <- c("com_anx", "com_ocd", "com_dep")
df$com_lumped <- apply(X = df[, com_vars], MARGIN = 1, FUN = function(x) {
  if(all(is.na(x))) return(NA)
  else if(any(x == "yes", na.rm = TRUE)) return(1)
  else return(0)
})
df$com_lumped <- factor(df$com_lumped, levels = c(0, 1), labels = c("no", "yes"))
rm(com_vars)


# Generate a lumped exclusion variable for healthy controls ---------------

hc_excl <- c("ksads_1_840_p", # Major Depressive Disorder Present
"ksads_1_841_p", # Major Depressive Disorder Current in Partial Remission
"ksads_2_835_p", # Bipolar II Disorder currently hypomanic
"ksads_2_836_p", # Bipolar II Disorder currently depressed
"ksads_2_831_p", # Bipolar I Disorder current episode depressed
"ksads_2_832_p", # Bipolar I Disorder currently hypomanic
"ksads_2_830_p", # Bipolar I Disorder current episode manic
"ksads_10_869_p", # Generalized Anxiety Disorder Present
"ksads_11_917_p", # Obsessive-Compulsive Disorder Present
"ksads_13_929_p", # Anorexia Nervosa Binge eating/purging subtype PRESENT
"ksads_13_930_p", # Anorexia Nervosa Binge eating/purging subtype PRESENT IN PARTIAL REMISSION
"ksads_13_932_p", # Anorexia Nervosa Restricting subtype PRESENT
"ksads_13_933_p", # Anorexia Nervosa Restricting subtype PRESENT IN PARTIAL REMISSION
"ksads_13_935_p", # Bulimia Nervosa PRESENT
"ksads_13_937_p", # Bulimia Nervosa Partial Remission Present
"ksads_15_901_p", # Oppositional Defiant Disorder Present
"ksads_16_897_p", # Conduct Disorder present childhood onset
"ksads_21_921_p") # Post-Traumatic Stress Disorder PRESENT

df$hc_excl <- apply(X = df[, hc_excl], MARGIN = 1, FUN = function(x) {
  if(all(is.na(x))) return(NA)
  else if(any(x == 1, na.rm = TRUE)) return(1)
  else return(0)
})
df$hc_excl <- factor(df$hc_excl, levels = c(0, 1), labels = c("no", "yes"))

# clean up
rm(hc_excl)


# Subset the data to the relevant predictors and outcomes -----------------

# predictors
pred <- c("ksads_14_853_p",   # ADHD current (parent report)
          "cbcl_scr_dsm5_adhd_t",  # CBCL DSM 5 ADHD subscale t score
          "bpmt_scr_attention_t", # BPM attention subscale t score
          "interview_age",    # age
          "sex",              # sex
          "race_ethnicity",   # race/ethnicity (White, Black, Hispanic, Asian, other)
          "high.educ",        # parental highest education
          "household.income", # parental combined income
          "married",          # parental marital status
          "pds_p",            # pubertal development scale (generated, parent report)
          "ksads_15_901_p",   # ODD
          "com_lumped",       # depression/anxiety diagnosis (generated, parent report)
          "hc_excl",          # exclusion criteria met for hcs
          "nihtbx_totalcomp_fc", # nih toolbox composite fully corrected t score
          "fsqc_qu_motion",   # motion (fixed)
          "site_id_l",        # site
          "mri_info_device.serial.number", # MRI serial number
          "rel_family_id")    # family

brain_vars <- names(df)[grep("smri_area.*desikan|smri_thick.*desikan|smri.*vol_subcort", names(df))]

# do not include ventricles, white matter, hypointensities, some brainstem 
# measures, and some summary measures
drop <- c("smri_vol_subcort.aseg_lateral.ventricle.lh",
          "smri_vol_subcort.aseg_lateral.ventricle.rh",
          "smri_vol_subcort.aseg_inf.lat.vent.lh",
          "smri_vol_subcort.aseg_inf.lat.vent.rh",
          "smri_vol_subcort.aseg_3rd.ventricle",
          "smri_vol_subcort.aseg_4th.ventricle",
          "smri_vol_subcort.aseg_cerebellum.white.matter.lh",
          "smri_vol_subcort.aseg_cerebellum.white.matter.rh",
          "smri_vol_subcort.aseg_cerebellum.cortex.lh",  
          "smri_vol_subcort.aseg_cerebellum.cortex.rh",
          "smri_vol_subcort.aseg_brain.stem",
          "smri_vol_subcort.aseg_csf",
          "smri_vol_subcort.aseg_cerebral.white.matter.lh",
          "smri_vol_subcort.aseg_cerebral.white.matter.rh",
          "smri_vol_subcort.aseg_wm.hypointensities",
          "smri_vol_subcort.aseg_wm.hypointensities.lh",
          "smri_vol_subcort.aseg_wm.hypointensities.rh",
          "smri_vol_subcort.aseg_cc.posterior",                 
          "smri_vol_subcort.aseg_cc.mid.posterior",             
          "smri_vol_subcort.aseg_cc.central",                   
          "smri_vol_subcort.aseg_cc.mid.anterior",              
          "smri_vol_subcort.aseg_cc.anterior",
          "smri_vol_subcort.aseg_wholebrain",
          "smri_vol_subcort.aseg_latventricles",
          "smri_vol_subcort.aseg_allventricles",
          "smri_vol_subcort.aseg_supratentorialvolume",      
          "smri_vol_subcort.aseg_subcorticalgrayvolume")

brain_vars <- brain_vars[!(brain_vars %in% drop)]

# merge the predictors (pred) and outcomes (brain_vars)
keep <- c(pred, brain_vars)

# subset the dataset to the relevant outcome and predictor variables
df <- subset(x = df, select = keep)

# clean up
rm(pred, brain_vars, drop, keep)


# Average left and right hemisphere measures ------------------------------

# loop over all variables
for (i in 1:ncol(df)) {
  
  # identify if a variable ends in lh ("left hemisphere")
  if (grepl(pattern = "lh$", x = names(df)[i])) {
    
    # generate a new variable name; delete .lh and append .avg
    new.var <- paste0(substr(x = names(df)[i], start = 1, stop = nchar(names(df)[i])-3), ".avg")
    
    # identify the corresponding right hemisphere variable; delete .lh and append .rh
    rh.var <- paste0(substr(x = names(df)[i], start = 1, stop = nchar(names(df)[i])-3), ".rh")
    
    # take the average
    df[, new.var] <- (df[, i] + df[, rh.var]) / 2
    
  }
}

rm (i, new.var, rh.var)


# Clean up variable names for simplicity ----------------------------------

# rename non-brain variables
names(df)[1:18] <- c("adhd", "cbcl", "bpm", "age", "sex", "race_eth", "educ", 
                     "income", "married", "pds_p", "odd", "com_lumped", 
                     "hc_excl", "tbx_comp", "motion", "site", "serial", "family")

# rename estimated total intracranial volume (etiv)
names(df)[names(df) == "smri_vol_subcort.aseg_intracranialvolume"] <- "etiv"


# Format the variables ----------------------------------------------------

df$adhd <- factor(x = df$adhd, levels = c(0, 1), labels = c("no", "yes"))
df$odd <- factor(x = df$odd, levels = c(0,1), labels = c("no", "yes"))
df$pds_p <- as.factor(df$pds_p)
df$motion <- as.factor(df$motion)
df$family <- as.factor(df$family)


# Impute missing demographic data -----------------------------------------

# save the original dataframe with missing data
df_orig <- df

# do not use serial, family, or imaging data for imputation; we've already 
# excluded participants with missing adhd and site data, so no imputation 
# will be done for those variables
imp_vars <- c("adhd", "age", "sex", "race_eth", "educ", "income", "married",
              "pds_p", "com_lumped", "tbx_comp", "motion", "site")

# create a dataframe with only complete rows (Enigma covariates)
df_orig_comp_reduced <- df_orig[complete.cases(df_orig[, c("age", "sex")]), ]

# creata a dataframe with only complete rows (all covariates)
df_orig_comp_full <- df_orig[complete.cases(df_orig[, imp_vars]), ]

# create a new dataset to hold the imputed data
df_imp_merge <- df

# impute missing values among and using the imp_vars above
out_mice <- mice(data = df_orig[, imp_vars], m = 15, seed = 123)

# generate a "merged" imputed dataset by taking the mean the imputed results for
# numeric variables and the mode for factors

# stack the imputed datasets into an array (all data are now characters)
out_mice_stacked <- abind(complete(data = out_mice, action = "all"), along = 3)

# numeric and ordinal variables - replace the NAs in the original data with the
# (rounded) mean
df_imp_merge$pds_p[is.na(df_imp_merge$pds_p)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$pds_p)), "pds_p", ], MARGIN = 1, FUN = function(x) round(mean(as.numeric(x))))
df_imp_merge$tbx_comp[is.na(df_imp_merge$tbx_comp)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$tbx_comp)), "tbx_comp", ], MARGIN = 1, FUN = function(x) round(mean(as.numeric(x))))
df_imp_merge$motion[is.na(df_imp_merge$v)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$motion)), "motion", ], MARGIN = 1, FUN = function(x) round(mean(as.numeric(x))))

# factors - replace the NAs in the original data with the mode
df_imp_merge$sex[is.na(df_imp_merge$sex)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$sex)), "sex", ], MARGIN = 1, FUN = mode_factor)
df_imp_merge$race_eth[is.na(df_imp_merge$race_eth)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$race_eth)), "race_eth", ], MARGIN = 1, FUN = mode_factor)
df_imp_merge$educ[is.na(df_imp_merge$educ)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$educ)), "educ", ], MARGIN = 1, FUN = mode_factor)
df_imp_merge$income[is.na(df_imp_merge$income)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$income)), "income", ], MARGIN = 1, FUN = mode_factor)
df_imp_merge$married[is.na(df_imp_merge$married)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$married)), "married", ], MARGIN = 1, FUN = mode_factor)
df_imp_merge$com_lumped[is.na(df_imp_merge$com_lumped)] <- apply(X = out_mice_stacked[which(is.na(df_imp_merge$com_lumped)), "com_lumped", ], MARGIN = 1, FUN = mode_factor)

# create a list of the m imputed datasets; add imputed data back to the
# original dataframe
for (i in 1:out_mice$m) {
  if (i == 1) list_df_imp <- list()
  temp <- df_orig
  temp[, imp_vars] <- complete(data = out_mice, action = i)
  list_df_imp[[i]] <- temp
}

# clean up
rm(df, imp_vars, out_mice_stacked, i, temp)


# Main analysis: model calls ----------------------------------------------

out_imp_merge_full <- model(df = df_imp_merge, version = "full")
out_imp_merge_full$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_full, raw_p = "p_value")

out_imp_merge_reduced <- model(df = df_imp_merge, version = "reduced")
out_imp_merge_reduced$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_reduced, raw_p = "p_value")

# write.csv(x = out_imp_merge_full, file = "../tables/out_imp_merge_full.csv")
# write.csv(x = out_imp_merge_reduced, file = "../tables/out_imp_merge_reduced.csv")


# Main analysis: interaction models (sex, comorbidity) --------------------

out_imp_merge_full_sex <- model_int_sex(df = df_imp_merge)
out_imp_merge_full_sex$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_full_sex, raw_p = "int_p_value")

out_imp_merge_full_com <- model_int_com(df = df_imp_merge)
out_imp_merge_full_com$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_full_com, raw_p = "int_p_value")


# Main analysis: LR test (race/ethnicity) ---------------------------------

out_imp_merge_full_race <- model_int_race(df = df_imp_merge)
out_imp_merge_full_race$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_full_race, raw_p = "int_p_value")


# Sensitivity analysis: imputed versus complete data ----------------------

out_orig_comp_full <- model(df = df_orig_comp_full, version = "full")
out_orig_comp_full$p_value_fdr_group <- make_fdr_group(df = out_orig_comp_full, raw_p = "p_value")

out_orig_comp_reduced <- model(df = df_orig_comp_reduced, version = "reduced")
out_orig_comp_reduced$p_value_fdr_group <- make_fdr_group(df = out_orig_comp_reduced, raw_p = "p_value")

# write.csv(x = out_orig_comp_full, file = "../tables/out_orig_comp_full.csv")
# write.csv(x = out_orig_comp_reduced, file = "../tables/out_orig_comp_reduced.csv")


# Sensitivity analysis: "healthy" controls ------------------------------

df_imp_merge_hc_excl <- df_imp_merge
df_imp_merge_hc_excl <- df_imp_merge_hc_excl[!(df_imp_merge_hc_excl$adhd == "no" & df_imp_merge_hc_excl$hc_excl == "yes"), ]

out_imp_merge_hc_excl_full <- model(df = df_imp_merge_hc_excl, version = "full")
out_imp_merge_hc_excl_full$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_hc_excl_full, raw_p = "p_value")

# out_imp_merge_hc_excl_reduced <- model(df = df_imp_merge_hc_excl, version = "reduced")
# out_imp_merge_hc_excl_reduced$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_hc_excl_reduced, raw_p = "p_value")

# write.csv(x = out_imp_merge_hc_excl_full, file = "../tables/out_imp_merge_hc_excl_full.csv")
# write.csv(x = out_imp_merge_hc_excl_reduced, file = "../tables/out_imp_merge_hc_excl_reduced.csv")


# Sensitivity analysis: CBCL instead of KSADS -----------------------------

df_imp_merge_cbcl <- df_imp_merge
df_imp_merge_cbcl$adhd <- factor(x = as.numeric(df_imp_merge_cbcl$cbcl >=65), levels = c(0, 1), labels = c("no", "yes"))
df_imp_merge_cbcl <- df_imp_merge_cbcl[!(is.na(df_imp_merge_cbcl$adhd)),]

out_imp_merge_cbcl_full <- model(df = df_imp_merge_cbcl, version = "full")
out_imp_merge_cbcl_full$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_cbcl_full, raw_p = "p_value")

# out_imp_merge_cbcl_reduced <- model(df = df_imp_merge_cbcl, version = "reduced")
# out_imp_merge_cbcl_reduced$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_cbcl_reduced, raw_p = "p_value")

# write.csv(x = out_imp_merge_cbcl_full, file = "../tables/out_imp_merge_cbcl_full.csv")
# write.csv(x = out_imp_merge_cbcl_reduced, file = "../tables/out_imp_merge_cbcl_reduced.csv")


# Sensitivity analysis: CBCL + BPM instead of KSADS -----------------------

df_imp_merge_bpm <- df_imp_merge
df_imp_merge_bpm <- df_imp_merge_bpm[!(is.na(df_imp_merge_bpm$bpm)),]

out_imp_merge_bpm_full_ksads <- model(df = df_imp_merge_bpm, version = "full")
out_imp_merge_bpm_full_ksads$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_bpm_full_ksads, raw_p = "p_value")

df_imp_merge_bpm$adhd <- factor(x = as.numeric(df_imp_merge_bpm$adhd =="yes" & df_imp_merge_bpm$bpm >=65), levels = c(0, 1), labels = c("no", "yes"))

out_imp_merge_bpm_full_both_k <- model(df = df_imp_merge_bpm, version = "full")
out_imp_merge_bpm_full_both_k$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_bpm_full_both_k, raw_p = "p_value")

df_imp_merge_bpm$adhd <- factor(x = as.numeric(df_imp_merge_bpm$cbcl >=65), levels = c(0, 1), labels = c("no", "yes"))

out_imp_merge_bpm_full_cbcl <- model(df = df_imp_merge_bpm, version = "full")
out_imp_merge_bpm_full_cbcl$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_bpm_full_ksads, raw_p = "p_value")

df_imp_merge_bpm$adhd <- factor(x = as.numeric(df_imp_merge_bpm$cbcl >=65 & df_imp_merge_bpm$bpm >=65), levels = c(0, 1), labels = c("no", "yes"))

out_imp_merge_bpm_full_both_c <- model(df = df_imp_merge_bpm, version = "full")
out_imp_merge_bpm_full_both_c$p_value_fdr_group <- make_fdr_group(df = out_imp_merge_bpm_full_both_c, raw_p = "p_value")


# Create model list for simulations ---------------------------------------

df_sim <- sim_diag(df = df_imp_merge, pct_fp = 0, pct_fn = 0)
model_list <- make_model_list(df_sim)  
rm(df_sim)

# Run a no misclassification simulation -----------------------------------

registerDoParallel(4)
set.seed(123)

out_none <- foreach(i = 1:1000) %dopar% {
  
  df_sim <- sim_diag(df = df_imp_merge, pct_fp = 0, pct_fn = 0)
  df_sim_brain <- sim_brain_data(df = df_sim, m_list = model_list)
  sim_out <- sim_model(df = df_sim_brain)
  return(sim_out)
  
  }

registerDoParallel()

save(out_none, file = "../rdata/sim_none.RData")


# Run the false-positive only simulation ----------------------------------

registerDoParallel(4)
set.seed(123)

out_fp <- foreach(i = 1:1000) %dopar% {

  df_sim <- sim_diag(df = df_imp_merge, pct_fp = .25, pct_fn = 0)
  df_sim_brain <- sim_brain_data(df = df_sim, m_list = model_list)
  sim_out <- sim_model(df = df_sim_brain)
  return(sim_out)
  
}

registerDoParallel()

save(out_fp, file = "../rdata/sim_fp.RData")


# Run the false-negative only simulation ----------------------------------
  
registerDoParallel(4)
set.seed(123)

out_fn <- foreach(i = 1:1000) %dopar% {

  df_sim <- sim_diag(df = df_imp_merge, pct_fp = 0, pct_fn = 0.1)
  df_sim_brain <- sim_brain_data(df = df_sim, m_list = model_list)
  sim_out <- sim_model(df = df_sim_brain)
  return(sim_out)
    
}

registerDoParallel()

save(out_fn, file = "../rdata/sim_fn.RData")

  
# Run the combined simulation ---------------------------------------------

registerDoParallel(4)
set.seed(123)
  
out_comb <- foreach(i = 1:100) %dopar% {   

  df_sim <- sim_diag(df = df_imp_merge, pct_fp = .25, pct_fn = 0.1)
  df_sim_brain <- sim_brain_data(df = df_sim, m_list = model_list)
  sim_out <- sim_model(df = df_sim_brain)
  return(sim_out)
  
}

registerDoParallel()

# save(out_comb, file = "../rdata/sim_comb.RData")


# Generate a dataframe with the effect sizes form the simulation ----------

df_none_eff_sizes <- make_df_eff_sizes(out_none)
df_fp_eff_sizes <- make_df_eff_sizes(out_fp)
df_fn_eff_sizes <- make_df_eff_sizes(out_fn)
df_comb_eff_sizes <- make_df_eff_sizes(out_comb)

# df_fp_eff_sizes[, c("brain_var", "enigma_e_s", "sim_min","sim_max", "sim_avg", "obs_min", "obs_max", "obs_avg")]


# Generate a dataframe with the p-values from the simulation --------------

df_none_p_vals <- make_df_p_vals(out_none)
df_fp_p_vals <- make_df_p_vals(out_fp)
df_fn_p_vals <- make_df_p_vals(out_fn)
df_comb_p_vals <- make_df_p_vals(out_comb)

# df_fp_p_vals[, c("brain_var", "enigma_e_s", "fdr_group_min","fdr_group_max", "group_always_sig")]


# Generate a dataframe with the distribution of significant results -------

df_none_distro_sig <- make_df_distro_sig(df_p_vals = df_none_p_vals)
df_fp_distro_sig <- make_df_distro_sig(df_p_vals = df_fp_p_vals)
df_fn_distro_sig <- make_df_distro_sig(df_p_vals = df_fn_p_vals)
df_comb_distro_sig <- make_df_distro_sig(df_p_vals = df_comb_p_vals)
