## Author: Codex
## Date updated: Feb 23 2026
#----------------------------------------------------------------------
# Train and tune a competing-risks model for 5- and 10-year CVD risk.
# Saves the best-performing model (by 5-year CV Brier) for downstream use.

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

build_formulas <- function() {
  list(
    linear = Hist(time, event) ~ sex + age + diabetes_duration + value_SBP +
      value_LDL + value_HBA1C + value_Smoking + value_Motion +
      macro_Albuminuria + micro_Albuminuria + eGFR + age_cat,
    spline = Hist(time, event) ~ sex + ns(age, 4) + ns(diabetes_duration, 3) +
      ns(value_SBP, 4) + ns(value_LDL, 3) + ns(value_HBA1C, 3) +
      value_Smoking + value_Motion + macro_Albuminuria + micro_Albuminuria +
      ns(eGFR, 3) + age_cat,
    spline_interaction = Hist(time, event) ~ sex + ns(age, 4) +
      ns(diabetes_duration, 3) + ns(value_SBP, 4) + ns(value_LDL, 3) +
      ns(value_HBA1C, 3) + value_Smoking + value_Motion +
      macro_Albuminuria + micro_Albuminuria + ns(eGFR, 3) + age_cat +
      age_cat:ns(eGFR, 3)
  )
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

score_candidate <- function(fit, data, seed) {
  sc <- Score(
    list("candidate" = fit),
    data = data,
    formula = Hist(time, event) ~ 1,
    times = 5,
    cause = 1,
    metrics = "brier",
    split.method = "cv5",
    seed = seed,
    null.model = FALSE,
    se.fit = FALSE,
    verbose = 0
  )
  brier <- sc$Brier$score[model == "candidate" & times == 5, Brier]
  brier
}

train_best_model <- function(data) {
  formulas <- build_formulas()
  pseudo_formula <- build_pseudo_formula()
  candidates <- list()

  fit_csc_spline <- CSC(formula = formulas$spline, data = data, cause = 1)
  fit_csc_spline$call$formula <- formulas$spline
  candidates <- append(candidates, list(list(
    name = "CSC_spline",
    type = "CSC",
    fit = fit_csc_spline,
    params = list(formula = "spline")
  )))

  fit_csc_interaction <- CSC(formula = formulas$spline_interaction, data = data, cause = 1)
  fit_csc_interaction$call$formula <- formulas$spline_interaction
  candidates <- append(candidates, list(list(
    name = "CSC_spline_interaction",
    type = "CSC",
    fit = fit_csc_interaction,
    params = list(formula = "spline_interaction")
  )))

  fit_fgr_linear <- FGR(formula = formulas$linear, data = data, cause = 1)
  fit_fgr_linear$call$formula <- formulas$linear
  candidates <- append(candidates, list(list(
    name = "FGR_linear",
    type = "FGR",
    fit = fit_fgr_linear,
    params = list(formula = "linear")
  )))

  glmnet_alpha_grid <- c(1, 0.5, 0)
  for (alpha in glmnet_alpha_grid) {
    fit_pseudo <- fit_pseudo_glmnet(
      data = data,
      times = c(5, 10),
      formula = pseudo_formula,
      alpha = alpha
    )
    candidates <- append(candidates, list(list(
      name = sprintf("PseudoGLMnet_alpha_%s", format(alpha, trim = TRUE)),
      type = "pseudo_glmnet",
      fit = fit_pseudo,
      params = list(alpha = alpha, times = "5,10")
    )))
  }

  results <- rbindlist(lapply(candidates, function(candidate) {
    brier <- score_candidate(candidate$fit, data, seed = 20260223)
    data.table(
      model = candidate$name,
      model_type = candidate$type,
      brier_5y = brier,
      params = paste(unlist(candidate$params), collapse = ";")
    )
  }))
  setorder(results, brier_5y)
  best_name <- results$model[1]
  best_fit <- candidates[[which(vapply(candidates, function(x) x$name, character(1)) == best_name)]][["fit"]]

  list(best_fit = best_fit, tuning = results)
}

main <- function() {
  data_path <- "training-sample.csv"
  models_dir <- "models"
  if (!dir.exists(models_dir)) {
    dir.create(models_dir, recursive = TRUE)
  }

  dt <- read_training_data(data_path)
  trained <- train_best_model(dt)

  best_fit <- trained$best_fit
  tuning <- trained$tuning

  saveRDS(best_fit, file.path(models_dir, "best_model.rds"))
  fwrite(tuning, file.path(models_dir, "model_candidate_results.csv"))

  # Sanity check for predictRisk compatibility
  predictRisk(best_fit, newdata = dt[1:10], times = c(5, 10), cause = 1)
}

main()
