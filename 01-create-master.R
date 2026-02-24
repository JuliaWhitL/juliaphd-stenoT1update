## Author: Julia Whitman
## Date updated: Feb 18 2026
#----------------------------------------------------------------------

### README ---
# The goal is to create a list of 'top 5' learners to inform the ones we should
# include in our Superlearner library for the Steno risk engine update.

## INSTRUCTIONS
# 1) Run the code below to generate the 'play' data 

# 2) Construct & tune learners that predict the 5-year risk of CVD (event== 1) in 
# the presence of non-CVD death (event== 2). Event==0 indicates right censoring.
# Optimal handling of outcome and covariates (number of splines, interactions, etc.) is dealer's choice, 
# but with specific covariate requirements detailed below.

# 3) Repeat for 10-year CVD risk

# 4) Send top performers (max 5) to Julia or upload the code to our shared 
# Github repo: https://github.com/JuliaWhitL/juliaphd-stenoT1update.


### COVARIATE REQUIREMENTS ---
# All models should include the following covariates, which were pre-selected in the original Steno analysis:
# sex (1 for female, 0 for male)
# age (years) 
# diabetes_duration (years)
# value_SBP (mmHg)
# value_LDL (mmol/l) 
# value_HBA1C (mmol/mol) 
# value_Smoking (1 for yes, 0 for no)
# value_motion (aka "no regular exercise": 1 for yes, 0 for no) 
# macro_Albuminuria (1 for yes, 0 for no) 
# micro_Albuminuria (1 for yes, 0 for no) 
# Interaction term between age_cat and log2(eGFR):
  # log2(eGFR) for age <40 years (mL/min/1.73m2)
  # log2(eGFR) for age >=40 years (mL/min/1.73m2)

### NOTE: SUBMITTING YOUR LEARNER ---
# Your submission should be a fitted model object that can be passed directly to
# riskRegression::Score(list("My Model" = your_fit_object), ...).
# Score internally calls predictRisk() on the object, so your model class must 
# have a compatible predictRisk S3 method.

# Standard riskRegression models (e.g. CSC, FGR) handle this automatically.

# For custom learners you must write a predictRisk S3 method for your model class. 
# This method must:
#   - Accept arguments: object, newdata, times, cause, ...
#   - Return a matrix of predicted risks with nrow = nrow(newdata) and ncol = length(times) 
    # i.e., 1 row per person, 1 column per prediction horizon (5- and 10-year models can be created & assessed separately if preferred)
#   - Predict cause-specific cumulative incidence 

# Minimal example for a custom learner:
# predictRisk.myModel <- function(object, newdata, times, cause, ...) {
#   # generate predicted cumulative incidence at each time point for each row of newdata
#   # return matrix: rows = observations, columns = time points
# }

# Before submitting, verify with:
  # predictRisk(your_fit, newdata = play, times = c(5, 10), cause = 1) 
  
  # riskRegression::Score(list("My Model" = your_fit), data = play, 
    # formula = Hist(time, event) ~ 1, times = c(5,10), cause = 1)

# or, for separate model/times combo... 
# riskRegression::Score(list("My Model" = your_fit5), data = play, 
                      # formula = Hist(time, event) ~ 1, times = 5, cause = 1)

# riskRegression::Score(list("My Model" = your_fit10), data = play, 
                      # formula = Hist(time, event) ~ 1, times = 10, cause = 1)

#----------------------------------------------------------------------

### Packages ---
library(lava)
library(data.table)
library(riskRegression)
library(survival)

### Generate 'play' data ---
# Returns data.table with baseline covariates and time to cvd (time,event) 
# where event has values 0 for right censored 1 for cvd and 2 for death without cvd. 
source("https://raw.githubusercontent.com/JuliaWhitL/juliaphd-stenoT1update/main/simulateStenoT1.R")
set.seed(111)
play <- simulateStenoT1(n=4000)

### Format for analysis ---
# create age category that can be interacted with eGFR as in original Steno model
play$age_cat <- factor(ifelse(play$age < 40, "young", "old"))

# create separate macro- and microAlbuminuria variables as in original Steno model
play$macro_Albuminuria <- ifelse(play$value_Albuminuria=="Macro", 1, 0)
play$micro_Albuminuria <- ifelse(play$value_Albuminuria=="Micro", 1, 0)

# remove intermediate variables and  info past censoring time
play <- play[, -c('time.event.1', 'time.event.2', 'time.event.0', 
                  'uncensored_time', 'uncensored_event')]


## EXAMPLE ---
# fit_bench <- CSC(Hist(time, event) ~ age + sex + diabetes_duration, data = play, cause = 1)
# 
# iv <- riskRegression::Score(list(
#   "benchmark" = fit_bench
#   # , model 2, etc.
#   ),
#    data = play,
#    seed = 8,
#    null.model = TRUE,
#    formula = Hist(time, event) ~ 1,
#    times = c(5,10), 
#    cause = 1,
#    split.method = "cv10",
#    summary = "ipa"
#  )
#  summary(iv, what = 'score')
 
