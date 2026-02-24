## test

### Packages ---
library(lava)
library(data.table)
library(riskRegression)
library(survival)

source(here::here("simulateStenoT1.R"))
newdata <- simulateStenoT1(n=2000)
newdata$macro_Albuminuria <- ifelse(newdata$value_Albuminuria=="Macro", 1, 0)
newdata$micro_Albuminuria <- ifelse(newdata$value_Albuminuria=="Micro", 1, 0)
newdata$age_cat <- factor(ifelse(newdata$age < 40, "young", "old"))


source("mikes-model/03-create-mike-model.R")
source("machines-model/03-create-codex-model.R")

Score(list(mikes_fit, codex_fit), 
  data = newdata,
    seed = 8,
  null.model = TRUE,
    formula = Hist(time, event) ~ 1,
    times = c(5,10), 
    cause = 1,
    split.method = "none",
    summary = "ipa"
  )

