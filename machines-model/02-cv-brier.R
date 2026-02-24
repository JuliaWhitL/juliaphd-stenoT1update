## Author: Codex
## Date updated: Feb 23 2026
#----------------------------------------------------------------------
# Cross-validated Brier score estimates for 5- and 10-year CVD risk.

library(data.table)
library(splines)
library(riskRegression)
library(prodlim)
library(survival)
library(glmnet)

set.seed(20260223)

read_training_data <- function(path) {
  dt <- fread(path)
  if (names(dt)[1] %in% c("", "X", "V1")) {
    dt[, (names(dt)[1]) := NULL]
  }

  if ("time" %in% names(dt)) {
    dt[, time := as.numeric(time)]
  }
  if ("event" %in% names(dt)) {
    dt[, event := as.integer(as.character(event))]
  }

  binary_cols <- c(
    "sex",
    "value_Smoking",
    "value_Motion",
    "macro_Albuminuria",
    "micro_Albuminuria"
  )
  for (col in binary_cols) {
    if (col %in% names(dt)) {
      dt[[col]] <- as.integer(as.character(dt[[col]]))
    }
  }

  if ("age_cat" %in% names(dt)) {
    dt[, age_cat := factor(age_cat)]
  }
  if ("value_Albuminuria" %in% names(dt)) {
    dt[, value_Albuminuria := factor(value_Albuminuria)]
  }

  dt
}

build_pseudo_formula <- function() {
  ~ sex + ns(age, 4) + ns(diabetes_duration, 3) + ns(value_SBP, 4) +
    ns(value_LDL, 3) + ns(value_HBA1C, 3) + value_Smoking + value_Motion +
    macro_Albuminuria + micro_Albuminuria + ns(eGFR, 3) + age_cat +
    age_cat:ns(eGFR, 3)
}

make_design_matrix <- function(formula, data, xlevels = NULL, contrasts = NULL) {
  mf <- model.frame(formula, data = data, xlev = xlevels, na.action = na.pass)
  mm <- model.matrix(formula, data = mf, contrasts.arg = contrasts)
  mm <- mm[, -1, drop = FALSE]
  mm
}

fit_pseudo_glmnet <- function(data, times, formula, alpha = 1) {
  terms_obj <- terms(formula)
  xlevels <- .getXlevels(terms_obj, model.frame(formula, data = data))
  contrasts <- attr(model.matrix(formula, data = data), "contrasts")
  x <- make_design_matrix(formula, data, xlevels, contrasts)
  cif_fit <- prodlim(Hist(time, event) ~ 1, data = data)
  pseudo <- prodlim::jackknife(cif_fit, cause = 1, times = times)

  models <- lapply(seq_along(times), function(i) {
    cv.glmnet(
      x = x,
      y = pseudo[, i],
      alpha = alpha,
      family = "gaussian",
      standardize = TRUE
    )
  })

  out <- list(
    models = models,
    times = times,
    formula = formula,
    xlevels = xlevels,
    contrasts = contrasts,
    alpha = alpha
  )
  out$call <- match.call()
  out$call$formula <- formula
  out$call$times <- times
  out$call$alpha <- alpha
  class(out) <- "pseudoGlmnet"
  out
}

predictRisk.pseudoGlmnet <- function(object, newdata, times, cause, ...) {
  if (!identical(cause, 1)) {
    stop("pseudoGlmnet only supports cause = 1.")
  }
  if (length(times) == 0) {
    return(matrix(numeric(0), nrow = nrow(newdata), ncol = 0))
  }
  x <- make_design_matrix(object$formula, newdata, object$xlevels, object$contrasts)
  train_times <- object$times
  preds_train <- sapply(seq_along(train_times), function(i) {
    as.numeric(predict(object$models[[i]], newx = x, s = "lambda.min"))
  })
  colnames(preds_train) <- as.character(train_times)

  interp_one_time <- function(tt) {
    if (tt %in% train_times) {
      preds_train[, as.character(tt)]
    } else {
      apply(preds_train, 1, function(p) {
        approx(train_times, p, xout = tt, rule = 2)$y
      })
    }
  }

  preds_out <- sapply(times, interp_one_time)
  preds_out <- pmin(pmax(preds_out, 0), 1)
  if (is.null(dim(preds_out))) {
    preds_out <- matrix(preds_out, nrow = nrow(newdata), ncol = length(times))
  }
  preds_out
}

main <- function() {
  data_path <- "training-sample.csv"
  model_path <- file.path("models", "best_model.rds")
  results_dir <- "results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  dt <- read_training_data(data_path)
  fit <- readRDS(model_path)

  cv_score <- Score(
    list("BestModel" = fit),
    data = dt,
    formula = Hist(time, event) ~ 1,
    times = c(5, 10),
    cause = 1,
    metrics = "brier",
    split.method = "cv10",
    seed = 20260223,
    null.model = FALSE,
    se.fit = FALSE,
    verbose = 0
  )

  brier_scores <- cv_score$Brier$score[model == "BestModel", .(times, Brier)]
  fwrite(brier_scores, file.path(results_dir, "cv_brier_scores.csv"))
  print(brier_scores)
}

main()
