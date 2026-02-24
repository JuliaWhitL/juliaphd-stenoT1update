# Model search summary

We evaluated several competing-risks learners to predict 5- and 10-year CVD risk (event = 1) with non‑CVD death as a competing event (event = 2). The goal was to balance flexibility (capturing non‑linear effects and interactions) with calibrated risk predictions that are compatible with `riskRegression::Score`.

## Models tried and why

- **CSC (cause-specific Cox) with splines and age_cat × eGFR interaction**: This is a strong, interpretable baseline for competing risks that allows non‑linear covariate effects and a clinically motivated interaction. It is directly supported by `riskRegression` and provides calibrated cumulative incidence predictions.
- **CSC with splines (no interaction)**: A slightly simpler variant to test whether the interaction adds predictive value.
- **FGR (Fine–Gray) linear model**: A standard subdistribution hazard model for competing risks; included as a classical alternative.
- **Pseudo‑value + glmnet (elastic net / lasso / ridge)**: A machine‑learning style approach that uses pseudo‑values of the cumulative incidence at 5 and 10 years and fits penalized regression with non‑linear terms and interactions. This provides flexible regularization and can improve generalization when many correlated predictors exist.

## Why the best model was selected

We selected the model with the lowest 5‑year cross‑validated Brier score (cv5 for tuning, then cv10 for final evaluation). The **CSC spline‑interaction** model achieved the best 5‑year Brier score among candidates and maintained strong 10‑year performance. This suggests the interaction between age category and eGFR adds predictive signal and that the spline structure captures non‑linear effects without overfitting. The pseudo‑glmnet models were close but slightly worse on the 5‑year Brier score, so the CSC spline‑interaction model was chosen as the best performer.
