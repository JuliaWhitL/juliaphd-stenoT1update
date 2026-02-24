
source(here::here("mikes-model", "00-model-specs.R"))

fullfits <- readRDS(here::here("mikes-model", "fullfits.rds"))
opt.fit <- readRDS(here::here("mikes-model", "opt-coeffs.rds"))

mikes_fit <- structure(list(fullfits = fullfits, opt.coeff = opt.fit), class = "mike")

predictRisk.mike <- function(object, newdata, times = c(5, 10), cause = 1, ...) {
  
  newdata <- as.data.frame(newdata)
  newdata$age_cat <- factor(ifelse(newdata$age < 40, "young", "old"))
  
  newdata$macro_Albuminuria <- ifelse(newdata$value_Albuminuria=="Macro", 1, 0)
  newdata$micro_Albuminuria <- ifelse(newdata$value_Albuminuria=="Micro", 1, 0)
  newdata$YY <- survival::Surv(newdata$time, newdata$event == 1)
  newdata$sex <- as.numeric(as.character(newdata$sex))
  newdata$value_Smoking <- as.numeric(as.character(newdata$value_Smoking))
  newdata$value_Motion <- as.numeric(as.character(newdata$value_Motion))

  newdata <- newdata[,c("sex", "age", "diabetes_duration", "value_SBP", "value_LDL", 
    "value_HBA1C", "value_Smoking", "value_Motion", "value_Albuminuria", 
    "eGFR", "age_cat", "macro_Albuminuria", "micro_Albuminuria", "YY")]
  
  Zvalid <- do.call(cbind, lapply(object$fullfits, function(f) f(newdata)))
  
  ptrunc <- function(x) {
    pmax(0, pmin(x, 1))
  }
  remat <- cbind(ptrunc(cbind(1, Zvalid) %*% object$opt.coeff$opt.fit5),
        ptrunc(cbind(1, Zvalid) %*% object$opt.coeff$opt.fit10))
  
  if(length(times) == 1 && times == 5) {
    remat[, 1, drop = FALSE]
  } else if (length(times) == 1 && times == 10) {
    remat[, 2, drop = FALSE]
  } else if (all(times == c(5, 10))) {
    remat
  } else {
    stop("Can't predict at times ", times)
  }
  
}

