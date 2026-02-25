library(eventglm) ## remotes::install_github("sachsmc/eventglm")
library(survival)
library(randomForestSRC)
library(CoxBoost)  ## remotes::install_github("binderh/CoxBoost")
library(riskRegression)
library(e1071)
library(splines)
library(glmnet)
library(class)
library(butcher)
## models

## survival
############### Natively ############################
#' @description  Predictions based on those Machine Learning procedures in the library that allow for weights to be specified as an argument of the R function. No bagging occurs. This group of algorithms is denoted as  Native Weights
#' @param dset data set
#' @return  prediction
#' @rdname megalearner-internal
#'
stepwise <- function(dset) {

  cfit <- coxph(YY ~ ., data = dset)
  cfin <- step(cfit, trace = 0)
  function(valid, time = 5) {
    predict(cfin, newdata = valid,
            type = "lp")
  }

}

random.forest <- function(dset, time = 5) {

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL

  cfit <- butcher(rfsrc(Surv(time, status) ~ ., data = dset, save.memory = TRUE))
  function(valid, time = 5) {
    res <- predict(cfit, newdata = valid)
    1 - res$survival[, max(which(res$time.interest <= time))[1]]
  }

}

coxboost <- function(dset) {

  mm <- model.matrix(YY ~ ., data = dset)[, -1]
  cfit <- butcher(CoxBoost(time = dset$YY[, "time"], status = dset$YY[, "status"],
           x = mm))
  function(valid, time = 5) {
    mmv <- model.matrix(YY ~ ., data = valid)[, -1]
    res <- predict(cfit, newdata = mmv,
                   type = "lp")
    res[1, ]
  }

}

## binary

direct.binomial <- function(dset, time = 5 ) {
  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)
  dset$weights <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$time <- dset$status <- NULL

  cfit <- glm(ybin ~ bs(age, degree = 1, knots = c(40), Boundary.knots = c(18, 85)) + 
                bs(diabetes_duration, degree = 1, knots = c(20), Boundary.knots = c(0, 60)) + 
                value_SBP + value_LDL + 
                value_HBA1C + value_Smoking + value_Motion + value_Albuminuria + 
                sex + log(eGFR),
      data = dset, family = "binomial", weights = weights)
  cfin <- step(cfit, trace = 0)
  function(valid, time = 5) {
      predict(cfit, newdata = valid, type = "response")
  }

}

## bagging ipcw

svm <- function(dset, time = 5 ){

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)

  dset <- dset[!is.na(dset$ybin),]
  wtmp <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$samp.wts <- wtmp / sum(wtmp)
  dset$time <- dset$status <- NULL

  svmboot <- lapply(1:10, function(i) {

    dboot <- dset[sample(1:nrow(dset), nrow(dset),
                         replace = TRUE, prob = dset$samp.wts),]
    dboot$samp.wts <- NULL
    e1071::svm(ybin ~ ., data = dboot)

  })


  function(valid, time = 5) {

    vboots <- sapply(svmboot, function(cfit) {
      predict(cfit, newdata = valid)
    })
    rowMeans(vboots)

  }
}




## pseudo obs

pseudo.glm <- function(dset, time = 5 ){

  cfit <- cumincglm(YY ~ bs(age, degree = 1, knots = c(40), Boundary.knots = c(18, 85)) + 
                      bs(diabetes_duration, degree = 1, knots = c(20), Boundary.knots = c(0, 60)) + 
                      value_SBP + value_LDL + 
                      value_HBA1C + value_Smoking + value_Motion + value_Albuminuria + 
                      sex + log(eGFR), data = dset, time = time)
  function(valid, time = 5) {
    predict(cfit, newdata = valid, type = "response")
  }

}


pseudo.glmnet <- function(dset, time = 5 ) {

  cfit <- cumincglm(YY ~  sex * (bs(age, degree = 1, knots = c(40), Boundary.knots = c(18, 85)) + 
                                   bs(diabetes_duration, degree = 1, knots = c(20), Boundary.knots = c(0, 60)) + 
                                   value_SBP + value_LDL + 
                                   value_HBA1C + value_Smoking + value_Motion + value_Albuminuria + 
                                   log(eGFR)), data = dset, time = time, x = TRUE)

  cnfit <- butcher(cv.glmnet(cfit$x, cfit$y))

  function(valid, time = 5) {
    vx <- model.matrix(~ sex * (bs(age, degree = 1, knots = c(40), Boundary.knots = c(18, 85)) + 
                                  bs(diabetes_duration, degree = 1, knots = c(20), Boundary.knots = c(0, 60)) + 
                                  value_SBP + value_LDL + 
                                  value_HBA1C + value_Smoking + value_Motion + value_Albuminuria + 
                                  log(eGFR)), data = valid)
    predict(cnfit, newx = vx)[, 1]

  }

}


