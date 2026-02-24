library(pseudoloss) ## remotes::install_github("sachsmc/pseudoloss")
source("00-model-specs.R")
library(parallel)

Zout <- readRDS("Zout.rds")
fullfits <- readRDS("fullfits.rds")
folds <- readRDS("split.rds")
play <- read.csv("training-sample.csv")
play$X <- NULL
play$event <- as.factor(play$event)
play$YY <- survival::Surv(play$time, play$event == 1)
play$time <- NULL
play$event <- NULL
play$value_Albuminuria <- as.factor(play$value_Albuminuria)
play$age_cat <- as.factor(play$age_cat)
devel <- play
devel$PO <- cumincglm(YY ~ 1, data = devel, time = 5)$y
devel$PO10 <- cumincglm(YY ~ 1, data = devel, time = 10)$y

YYens <- devel$PO[unlist(folds)]
YYens10 <- devel$PO10[unlist(folds)]

Zmat <- do.call(rbind, Zout)

lassofit <- cv.glmnet(Zmat, YYens, standardize = FALSE, alpha = 0)
lfit <- glmnet(Zmat, YYens, standardize = FALSE,
               lambda = lassofit$lambda.1se, alpha = 0)
lassofit10 <- cv.glmnet(Zmat, YYens10, standardize = FALSE, alpha = 0)
lfit10 <- glmnet(Zmat, YYens10, standardize = FALSE,
               lambda = lassofit$lambda.1se, alpha = 0)


saveRDS(list(opt.fit5 = c(lfit$a0, lfit$beta[, 1]), 
             opt.fit10 = c(lfit10$a0, lfit10$beta[, 1])), "opt-coeffs.rds")



