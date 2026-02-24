source("00-model-specs.R")
library(parallel)

play <- read.csv("training-sample.csv")
play$X <- NULL
play$event <- as.factor(play$event)
play$YY <- survival::Surv(play$time, play$event == 1)
play$time <- NULL
play$event <- NULL
play$value_Albuminuria <- as.factor(play$value_Albuminuria)
play$age_cat <- as.factor(play$age_cat)
devel <- play

mymods <- list(stepwise, random.forest,
               direct.binomial, 
               function(xx) pseudo.glm(xx, time = 5 ))

set.seed(420)
ndex <- 1:nrow(devel)
part <- as.factor(sample(1:10, length(ndex), replace = TRUE))
folds <- split(ndex, part)

Zout <- vector(mode = "list", length = 10)
for(j in 1:length(folds)) {

  training <- devel[unlist(folds[-j]), ]
  validation <- devel[folds[[j]], ]

  Zout[[j]] <- do.call(cbind, mclapply(mymods, function(f){

    fhat <- f(training)
    fhat(validation)

  }, mc.cores = length(mymods)))

}


fullfits <- lapply(mymods, function(f) {

  f(devel)

})


saveRDS(fullfits, "fullfits.rds")
saveRDS(Zout, "Zout.rds")
saveRDS(folds, "split.rds")
