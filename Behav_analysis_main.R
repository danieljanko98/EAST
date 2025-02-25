#######################################################################
#                         EAST Main Model Fitting                     #
#######################################################################

# This script fits a linear Bayesian model on response time and Bayesian logistic regression on Accuracy. 
# Further processing of the results nad extraction of specific matrics is also included to allow for interpretation of the results. 

packages <- c("brms", "bayesplot", "hypr", "rstan", "ggplot2", "tidyr", "dplyr", 'lme4', 'stringr', 'ggridges')

# Function to check if packages are installed, install them if not, and load them
load_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      if (!require(package, character.only = TRUE)) {
        stop(paste("Package", package, "not found and installation failed. Please install it manually."))
      }
    }
    library(package, character.only = TRUE)
  }
}
load_packages(packages)


source('/Users/danieljanko/Desktop/Projects/EAST/functions/functions_EAST.R')
setwd('/Users/danieljanko/Desktop/Projects/EAST/')

#######################################################################
#                      Response Time Linear Model                     #
#######################################################################

data <- read.csv('EAST_Data_long.csv')

# Recoding some of the variables for a better interpretation
variableRecoding(data)

# Extracting posteriors from the pilot analysis to create informative priors
# Model used for the prior elicitation - RTlog ~ 1 + Congruency + Error + PrevAcc * PrevSubjAcc + (1 || SubInfo)
load("Bayes_model_full_RT.RData")
model_prior <- model8
rm(model8)
summary(model_prior)

# Setting up priors for the new model
inform_priors <- c(prior(normal(6.46,0.07), class = Intercept, lb = 4.6, ub = 8),
                   prior(normal(0.04, 0.02), class = b, coef = Congruency),
                   prior(normal(0.08, 0.02), class = b, coef = Error),
                   prior(normal(0, 2), class = b, coef = SleepOrder),
                   prior(normal(0.02, 0.02), class = b, coef = TrialType),
                   prior(normal(0, 2), class = b, coef = TSD),
                   prior(normal(0, 0.03), class = b, coef = PrevAcc),
                   prior(normal(0.14, 0.03), class = b, coef = PrevSubjAcc),
                   prior(normal(0, 2), class = b, coef = TSD:PrevAcc),
                   prior(normal(0, 2), class = b, coef = TSD:PrevSubjAcc),
                   prior(normal(-0.11, 0.04), class = b, coef = PrevAcc:PrevSubjAcc),
                   prior(normal(0, 2), class = b, coef = TSD:PrevAcc:PrevSubjAcc),
                   prior(normal(0.16, 0.05), class = sd, coef = Intercept, group = SubInfo),
                   prior(normal(0, 2), class = sigma)
)

main_model <- brm(RTlog ~ 1 + Congruency + Error + SleepOrder + TrialType + TSD * PrevAcc * PrevSubjAcc + (1 || SubInfo), data = data,
                  prior = inform_priors,
                  control = list(adapt_delta = 0.99, max_treedepth = 12),
                  save_pars = save_pars(all = TRUE),
                  iter = 20000)

summary(main_model)
save(main_model, file = "Bayes_model_main.RData")


load("Bayes_model_main.RData")
model <- main_model
rm(main_model)

# Exploring the results
modelExploration(model, inform_priors)

# Exploring the different outputs
print(ppCheck)
print(model_figure)
print(Results_RT)
# Specify which plot to show (indices based on the ordering of the priors)
print(plot_list_RT[1])

#######################################################################
#                        Accuracy Logistic Model                      #
#######################################################################

# Fitting the logistic regression
load("Bayes_model_log_Pilot.RData")
summary(log_model_pilot)
# Specifying the priors manually - important to specify them in the same order the model is later specified
inform_priors <- c(prior(normal(1.5,0.30), class = Intercept, coef = Intercept),
                   prior(normal(-0.15, 0.1), class = b, coef = Congruency),
                   prior(normal(-0.71, 0.1), class = b, coef = TrialType),
                   prior(normal(0, 5), class = b, coef = SleepOrder),
                   prior(normal(0, 5), class = b, coef = TSD),
                   prior(normal(0.02, 0.15), class = b, coef = PrevAcc),
                   prior(normal(-0.26, 0.15), class = b, coef = PrevSubjAcc),
                   prior(normal(0, 5), class = b, coef = TSD:PrevAcc),
                   prior(normal(0, 5), class = b, coef = TSD:PrevSubjAcc),
                   prior(normal(-0.47, 0.3), class = b, coef = PrevAcc:PrevSubjAcc),
                   prior(normal(0, 5), class = b, coef = TSD:PrevAcc:PrevSubjAcc),
                   prior(normal(0.65, 0.25), class = sd, coef = Intercept, group = SubInfo)
)

log_model_main <- brm(ACC ~ 1 + Congruency + TrialType + SleepOrder + TSD * PrevAcc * PrevSubjAcc + (1 || SubInfo),
                      prior = inform_priors,
                      data = data,
                      iter = 20000,
                      family = bernoulli(link = logit))
# Saving the model
save(log_model_main, file = "Bayes_model_log_main.RData")
# Load the model
load('Bayes_model_log_main.RData')
model <- log_model_main
LogModelExploration(model, inform_priors)

# Exploring the different outputs
print(Results_ACC)
# Specify which plot to show (indices based on the ordering of the priors)
print(plot_list_ACC[1])
print(ppCheck)
print(log_model_figure)






