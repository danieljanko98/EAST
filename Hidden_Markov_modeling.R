#######################################################################
#                         EAST Hidden Markov States                   #
#######################################################################

# Author: Daniel Janko 

# Goal
  # The main goal of this analysis is to see whether there are different states subjects are in that 
  # would predict their performance. Additionally, we want to see whether switching between states is affected
  # by sleep deprivation. We will also explore how different factors (PrevAcc, PrevSubjAcc etc.) influence performance 
  # in different states.

packages <- c("depmixS4", "dplyr")

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
setwd("/Users/danieljanko/Desktop/Projects/EAST")
data <- read.csv('EAST_Data_long.csv')


#######################################################################
#                   Preparing data for model fitting.                 #
#######################################################################

# The function below creates supporting data structures that are needed for a proper model fitting
timeSeries(data)

# Fitting a model that assumes 2 states, the starting state depends on the TSD and switching probability between
# states also depends on TSD
model <- depmix(list(ACC ~ 1, RTlog ~ 1),
                transition = ~ TSD,
                prior = ~ TSD,
                data = data, 
                nstates = 2,
                ntimes = count_vector,
                family = list(multinomial(), gaussian()), 
                initdata = data.frame(TSD = tsd_vector),
                instart = runif(4))

hmm_fit <- fit(model)
summary(hmm_fit)

# Fitting a different model that assumes the starting state depends on the session (first/second)

# Need to create a new variable session

sessionSpec(data)

model2 <- depmix(list(ACC ~ 1, RTlog ~ 1),
                transition = ~ TSD,
                prior = ~ session,
                data = data, 
                nstates = 2,
                ntimes = count_vector,
                family = list(multinomial(), gaussian()), 
                initdata = data.frame(session = session_vector),
                instart = runif(4))
hmm_fit2 <- fit(model2)
summary(hmm_fit2)

# Work in progress. More complex models will be added







