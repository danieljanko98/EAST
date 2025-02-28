 
variableRecoding <- function(data) {
  
  data$Congruency <- ifelse(data$Congruency == -1, -0.5, 0.5) # Congruent = -0.5
  data$PrevAcc <- ifelse(data$PrevAcc == -1, -0.5, 0.5) # Correct = -0.5
  data$PrevSubjAcc <- ifelse(data$PrevSubjAcc == -1, -0.5, 0.5) # Correct = -0.5
  data$Error <- ifelse(data$Error == -1, -0.5, 0.5) # Correct = -0.5
  data$TrialType <- ifelse(data$TrialType == -1, -0.5, 0.5) # Simple = -0.5
  data$SleepOrder <- ifelse(data$SleepOrder == 2, -0.5, 0.5) # WR first = -0.5
  data$TSD <- ifelse(data$TSD == 0, -0.5, 0.5) # WR = -0.5
  data <<- data
  
}


modelExploration <- function(model, inform_priors) {
  df_model_main <- as.data.frame(model)
  var_names_log <- colnames(df_model_main)[2:12]
  log_model_figure <- mcmc_intervals(model,
                                     regex_pars = var_names_log,
                                     prob_outer = .95,
                                     prob = .8,
                                     point_size = 2.5,
                                     inner_size = 1,
                                     point_est = "mean")
  print(log_model_figure)
  
  # Getting the probability of observing effect between -0.01 and 0.01 given the posterior distributions
  effects <- fixef(model)
  # Probability for Intercept, beta1, beta2, beta3
  probs <- vector()
  for (i in 1:12) {
    probs[i] <- round(pnorm(0.01, mean = effects[i,1], sd = effects[i,2]) - pnorm(-0.01, mean = effects[i,1], sd = effects[i,2]), 3)
  }
  
  # calculating the  Kullback-leibler divergence for all estimates but only those that had informative priors will be used later
  # calculaating for all due to easier implementation 
  plot_list_RT <- list()
  kl_results <- c()
  for (j in 1:12) {
    P <- as.numeric(str_extract_all(inform_priors[j,1], "[0-9\\.]+")[[1]])
    prior_dist <- rnorm(40000, P[1], P[2]) 
    post_dist <- df_model_main[,j]
    
    prior_df <- data.frame(value = prior_dist, dist = "Prior")
    posterior_df <- data.frame(value = post_dist, dist = "Posterior")
    df <- rbind(prior_df, posterior_df)
    plot_list_RT[[j]] <- ggplot(df, aes(x = value, fill = dist)) +
      geom_density(alpha = 0.5) +
      labs(title = paste("Prior vs Posterior Distributions", inform_priors[j,3])) +
      theme_minimal() +
      scale_fill_manual(values = c("blue", "red"))
    
    
    mu1 <- P[1]
    sigma1 <- P[2]
    mu2 <- effects[j,1]
    sigma2 <- effects[j,2]
    kl_results[j] <- round(KL_divergence(mu1, sigma1, mu2, sigma2),2)
  }
  
  ppCheck <<- pp_check(model, ndraws = 100, type = 'dens_overlay')
  
  model_figure <<- log_model_figure
  plot_list_RT <<- plot_list_RT
  
  
  results_overview <- rbind(probs, kl_results)
  colnames(results_overview) <- c('Intercept', c(inform_priors[2:12, 3]))
  rownames(results_overview) <- c('Probabilities', 'KL_divergence')
  Results_RT <<- results_overview
}

LogModelExploration <- function(model, inform_priors) {
  df_model_main <- as.data.frame(model)
  var_names_log <- colnames(df_model_main)[2:12]
  log_model_figure <- mcmc_intervals(model,
                                     regex_pars = var_names_log,
                                     prob_outer = .95,
                                     prob = .8,
                                     point_size = 2.5,
                                     inner_size = 1,
                                     point_est = "mean")
  print(log_model_figure)
  
  # Getting the probability of observing effect between -0.01 and 0.01 given the posterior distributions
  effects <- fixef(model)
  # Probability for Intercept, beta1, beta2, beta3
  probs <- vector()
  for (i in 1:11) {
    probs[i] <- round(pnorm(0.01, mean = effects[i,1], sd = effects[i,2]) - pnorm(-0.01, mean = effects[i,1], sd = effects[i,2]), 3)
  }
  
  # calculating the  Kullback-leibler divergence for all estimates but only those that had informative priors will be used later
  # calculating for all due to easier implementation 
  plot_list_ACC <- list()
  kl_results <- c()
  for (j in 1:11) {
    P <- as.numeric(str_extract_all(inform_priors[j,1], "[0-9\\.]+")[[1]])
    prior_dist <- rnorm(40000, P[1], P[2]) 
    post_dist <- df_model_main[,j]
    
    prior_df <- data.frame(value = prior_dist, dist = "Prior")
    posterior_df <- data.frame(value = post_dist, dist = "Posterior")
    df <- rbind(prior_df, posterior_df)
    plot_list_RT[[j]] <- ggplot(df, aes(x = value, fill = dist)) +
      geom_density(alpha = 0.5) +
      labs(title = paste("Prior vs Posterior Distributions", inform_priors[j,3])) +
      theme_minimal() +
      scale_fill_manual(values = c("blue", "red"))
    
    
    mu1 <- P[1]
    sigma1 <- P[2]
    mu2 <- effects[j,1]
    sigma2 <- effects[j,2]
    kl_results[j] <- round(KL_divergence(mu1, sigma1, mu2, sigma2),2)
  }
  
  ppCheck <<- pp_check(model, ndraws = 100, type = 'dens_overlay')
  
  log_model_figure <<- log_model_figure
  plot_list_RT <<- plot_list_RT
  
  
  results_overview <- rbind(probs, kl_results)
  results_overview <- rbind(results_overview, effects[,"Estimate"])
  colnames(results_overview) <- c('Intercept', c(inform_priors[2:11, 3]))
  rownames(results_overview) <- c('Probabilities', 'KL_divergence', 'Posterior Estimates')
  Results_ACC <<- results_overview
}

KL_divergence <- function(mu1, sigma1, mu2, sigma2) {
  term1 <- log(sigma2 / sigma1)
  term2 <- (sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2)
  term3 <- -0.5
  return(term1 + term2 + term3)
}

# the timeSeries function specifies the number of elements in each of the time series included in the data set. Since we have concatenated data, we need to 
# tell the model where one time series end and the other one begins. Otherwise the model will treat the whole data as one time series.
timeSeries <- function(data) {
  
  data$TSD <- ifelse(data$TSD == '0', 1, 2)
  data <- data %>% arrange(SubInfo, TSD)
  subj <- unique(data$SubInfo)
  count_vector <- c()
  vec <- seq(1, 60, by = 2)
  for (i in 1:length(subj)) {
    for (j in 1:2) {
      a <- c(0,1)
      count_vector[vec[i] + a[j]] <- data %>%
        filter(SubInfo == subj[i], TSD == j) %>%
        nrow()
    }
  }
  data <- as.data.frame(data)
  tsd_vector <<- rep(c(0, 1), length.out = 60)
  data$TSD <- ifelse(data$TSD == '1', 0, 1)
  
  data <<- data
  count_vector <<- count_vector
}

sessionSpec <- function(data) {
  for (i in 1:nrow(data)) {
    if (data$TSD[i] == 0 & data$SleepOrder[i] == 1) {
      data$session[i] <- 1
    }
    if (data$TSD[i] == 1 & data$SleepOrder[i] == 1) {
      data$session[i] <- 0
    }
    if (data$TSD[i] == 0 & data$SleepOrder[i] == 2) {
      data$session[i] <- 0
    }
    if (data$TSD[i] == 1 & data$SleepOrder[i] == 2) {
      data$session[i] <- 1
    }
    
  }
  
  session_vector <- rep(c(0, 1), length.out = 60)
  
  data <<- data
  session_vector <<- session_vector
  
  
  
  
  
}





