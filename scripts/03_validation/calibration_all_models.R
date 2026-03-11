# ============================================================
# CALIBRATION — All three models, clean slate
# Loads saved models + frozen scores, outputs Excel table
# ============================================================
# Models:
#   1) Berlin 7-protein GLM  → validated on Prague
#   2) Train60 7-protein GLM → validated on Test40
#   3) Train60 Concentration → validated on Test40
# ============================================================

library(pROC)
library(dplyr)
library(readr)
library(ResourceSelection)
library(boot)
library(openxlsx)
set.seed(7)

# ----------------------------------------------------------
# PATHS — adjust if needed
# ----------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
outdir       <- "../proteomics_ml/data/calibration_output/finalx"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# LOAD FROZEN SCORES
# ----------------------------------------------------------
scores_test40 <- read_csv(file.path(artifact_dir, "scores_test40_frozen.csv"))
scores_prague <- read_csv(file.path(artifact_dir, "scores_Prague_frozen.csv"))
scores_test40_con <- read_csv(file.path(artifact_dir, "scores_test40_frozen_concentration.csv"))

cat("Scores loaded:\n")
cat("  Prague n =", nrow(scores_prague), "\n")
cat("  Test40 n =", nrow(scores_test40), "\n")

# ----------------------------------------------------------
# CALIBRATION FUNCTION
# Returns one data.frame row per model
# ----------------------------------------------------------
run_calibration <- function(truth, pred, label) {
  
  truth <- as.integer(truth)
  pred  <- as.numeric(pred)
  
  # AUC + DeLong CI
  roc_obj <- pROC::roc(truth, pred, quiet = TRUE)
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
  auc_ci  <- round(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), 3)
  
  # Brier score + bootstrap CI
  brier_val  <- round(mean((pred - truth)^2), 3)
  brier_boot <- boot::boot(
    data      = data.frame(p = pred, y = truth),
    statistic = function(d, i) mean((d$p[i] - d$y[i])^2),
    R         = 2000
  )
  brier_ci <- round(quantile(brier_boot$t, c(0.025, 0.975)), 3)
  
  # Calibration regression (intercept + slope)
  logit_safe <- function(p) {
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    log(p / (1 - p))
  }
  cal_fit   <- glm(truth ~ logit_safe(pred), family = binomial())
  intercept <- round(coef(cal_fit)[1], 3)
  slope     <- round(coef(cal_fit)[2], 3)
  
  # Hosmer-Lemeshow
  hl_p <- tryCatch(
    round(ResourceSelection::hoslem.test(truth, pred, g = 10)$p.value, 3),
    error = function(e) NA_real_
  )
  
  # Calibration in the large
  mean_pred <- round(mean(pred),  3)
  mean_obs  <- round(mean(truth), 3)
  cal_itl   <- round(mean_pred - mean_obs, 3)
  
  # Label permutation p
  perm_aucs <- replicate(2000, {
    as.numeric(pROC::auc(pROC::roc(sample(truth), pred, quiet = TRUE)))
  })
  perm_p <- round(mean(perm_aucs >= auc_val), 4)
  
  data.frame(
    model            = label,
    n                = length(truth),
    events           = sum(truth),
    AUC              = auc_val,
    AUC_lower95      = auc_ci[1],
    AUC_upper95      = auc_ci[3],
    Brier            = brier_val,
    Brier_lower95    = brier_ci[1],
    Brier_upper95    = brier_ci[2],
    Permutation_p    = perm_p,
    Cal_intercept    = intercept,
    Cal_slope        = slope,
    HL_p             = hl_p,
    Mean_pred        = mean_pred,
    Observed_rate    = mean_obs,
    Cal_in_the_large = cal_itl,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# RUN ALL THREE MODELS
# ----------------------------------------------------------
cat("\nRunning calibration (bootstrap + permutation ~1 min each)...\n")

cal_berlin <- run_calibration(
  truth = scores_prague$group01,
  pred  = scores_prague$pred_berlin_7,
  label = "Berlin 7-protein -> Prague"
)

cal_60_7 <- run_calibration(
  truth = scores_test40$group01,
  pred  = scores_test40$pred_train60_7,
  label = "Train60 7-protein -> Test40"
)

cal_60_conc <- run_calibration(
  truth = scores_test40_con$group01,
  pred  = scores_test40_con$pred_train60_conc,
  label = "Train60 Concentration -> Test40"
)

# ----------------------------------------------------------
# COMBINE AND SAVE
# ----------------------------------------------------------
cal_all <- bind_rows(cal_berlin, cal_60_7, cal_60_conc)

print(cal_all)
row.names(cal_all) <- NULL
write_csv(cal_all, file.path(outdir, "calibration_all_models.csv"))

scores_berlin_train <- read_csv(file.path(artifact_dir, "scores_data_testB_frozen.csv"))

sep_berlin <- scores_berlin_train %>%
  summarise(
    total_n        = n(),
    pred_below_0001 = sum(pred_berlin_7 < 0.001),
    pred_above_0999 = sum(pred_berlin_7 > 0.999),
    n_extreme       = pred_below_0001 + pred_above_0999,
    pct_extreme     = round(100 * n_extreme / total_n, 1)
  )

print(sep_berlin)
write_csv(sep_berlin, file.path(outdir, "separation_berlin_train.csv"))
# Berlin (n=104): 53/104 = 51% extreme

scores_60_train <- read_csv(file.path(artifact_dir, "scores_train60_frozen.csv"))
sep_berlin <- scores_60_train %>%
  summarise(
    total_n        = n(),
    pred_below_0001 = sum(pred_train60_7 < 0.001),
    pred_above_0999 = sum(pred_train60_7 > 0.999),
    n_extreme       = pred_below_0001 + pred_above_0999,
    pct_extreme     = round(100 * n_extreme / total_n, 1)
  )

print(sep_berlin)
write_csv(sep_berlin, file.path(outdir, "separation_60_40_train.csv"))
# Train60 (n=127): 37/127 = 29%


scores_60_train_conc <- read_csv(file.path(artifact_dir, "scores_train60_concentration_frozen.csv"))
sep_berlin <- scores_60_train_conc %>%
  summarise(
    total_n        = n(),
    pred_below_0001 = sum(pred_train60_conc < 0.001),
    pred_above_0999 = sum(pred_train60_conc > 0.999),
    n_extreme       = pred_below_0001 + pred_above_0999,
    pct_extreme     = round(100 * n_extreme / total_n, 1)
  )

print(sep_berlin)
write_csv(sep_berlin, file.path(outdir, "separation_60_40_train_concentration.csv"))
# Train60 (n=127): 6/124 = 5%

# The separation problem is specific to the normalized/transformed protein data used in the Berlin and 
# Train60 7-protein models, and largely disappears when using raw concentrations. This could be worth a 
# sentence in your results or discussion, as it suggests the feature scaling/transformation may be contributing
# to separation on top of the small EPV.

# Fits Ridge on Berlin with 10-fold CV
library(ResourceSelection)
set.seed(7)
# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
data_csv     <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
outdir       <- "../proteomics_ml/data/calibration_output/finalx"
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir,       recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# PREDICTORS
# ----------------------------------------------------------
prots_7 <- c("AHSG_001", "CLEC3B_020", "COMP_023",
             "F9_027", "LRG1_042", "MCAM_047", "MRC1_073..")

# ----------------------------------------------------------
# LOAD DATA + SPLITS
# ----------------------------------------------------------
dfz <- read_csv(data_csv)

data_trainB <- dfz %>% filter(cohort == "Berlin")
data_testP  <- dfz %>% filter(cohort == "Prague")

load("../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r")  # full_meta_train
load("../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r")   # full_meta_test

train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)

cat("Berlin n =", nrow(data_trainB), "\n")
cat("Prague  n =", nrow(data_testP),  "\n")
cat("Train60 n =", nrow(train60),     "\n")
cat("Test40  n =", nrow(test40),      "\n")

dfvsn <- readRDS(file = "../proteomics_ml/data/variance_dfvsn_final.rds")
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_trainn.r")
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_test.r")

# Create VSN‑based train/test using the old variance_final_data_raw
vsn_data_based_test  <- dfvsn[dfvsn$cleaned %in% row.names(data_test), ]
vsn_data_based_train <- dfvsn[dfvsn$cleaned %in% row.names(data_train), ]

# ----------------------------------------------------------
# LOAD EXISTING GLM MODELS (for coefficient comparison)
# ----------------------------------------------------------
model_berlin_7  <- readRDS(file.path(artifact_dir, "model_7panel_Berlin_glm.rds"))
model_train60_7 <- readRDS(file.path(artifact_dir, "model_7panel_Train60_glm.rds"))
model_train60_7_conc <- readRDS(file.path(artifact_dir, "model_conc_Train60_glm.rds"))

# ----------------------------------------------------------
# FIT RIDGE MODEL ON BERLIN
# ----------------------------------------------------------
X_berlin <- as.matrix(data_trainB[, prots_7])
y_berlin  <- data_trainB$group01

cv_ridge <- glmnet::cv.glmnet(
  x      = X_berlin,
  y      = y_berlin,
  family = "binomial",
  alpha  = 0,       # Ridge
  nfolds = 10
)

cat("\nLambda selected (1se):", round(cv_ridge$lambda.1se, 4), "\n")
cat("Lambda selected (min):", round(cv_ridge$lambda.min,  4), "\n")

model_berlin_ridge <- glmnet::glmnet(
  x      = X_berlin,
  y      = y_berlin,
  family = "binomial",
  alpha  = 0,
  lambda = cv_ridge$lambda.1se
)

# ----------------------------------------------------------
# SCORE PRAGUE
# ----------------------------------------------------------
pred_glm   <- as.numeric(predict(model_berlin_7,
                                 newdata = data_testP,
                                 type    = "response"))

pred_ridge <- as.numeric(predict(model_berlin_ridge,
                                 newx  = as.matrix(data_testP[, prots_7]),
                                 type  = "response"))

cat("\nGLM   — Prague extreme predictions:\n")
cat("  pred < 0.001:", sum(pred_glm < 0.001), "| pred > 0.999:", sum(pred_glm > 0.999),
    "| total:", length(pred_glm), "\n")
# pred < 0.001: 29 | pred > 0.999: 17 | total: 108 

cat("Ridge — Prague extreme predictions:\n")
cat("  pred < 0.001:", sum(pred_ridge < 0.001), "| pred > 0.999:", sum(pred_ridge > 0.999),
    "| total:", length(pred_ridge), "\n")
# pred < 0.001: 0 | pred > 0.999: 0 | total: 108 

# ----------------------------------------------------------
# CALIBRATION FUNCTION
# ----------------------------------------------------------
run_calibration <- function(truth, pred, label, R_boot = 2000) {
  
  truth <- as.integer(truth)
  pred  <- as.numeric(pred)
  
  # ----------------------------------------------------------
  # AUC + DeLong CI
  # ----------------------------------------------------------
  roc_obj <- pROC::roc(truth, pred, quiet = TRUE)
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
  auc_ci  <- round(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), 3)
  
  # ----------------------------------------------------------
  # Brier score + bootstrap CI
  # ----------------------------------------------------------
  brier_val  <- round(mean((pred - truth)^2), 3)
  brier_boot <- boot::boot(
    data      = data.frame(p = pred, y = truth),
    statistic = function(d, i) mean((d$p[i] - d$y[i])^2),
    R         = R_boot
  )
  brier_ci <- round(quantile(brier_boot$t, c(0.025, 0.975)), 3)
  
  # ----------------------------------------------------------
  # Calibration regression (intercept + slope) — point estimates
  # ----------------------------------------------------------
  logit_safe <- function(p) {
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    log(p / (1 - p))
  }
  
  df_cal  <- data.frame(y = truth, p = pred)
  cal_fit <- glm(y ~ logit_safe(p), data = df_cal, family = binomial())
  
  intercept <- round(coef(cal_fit)[1], 3)
  slope     <- round(coef(cal_fit)[2], 3)
  
  # ----------------------------------------------------------
  # Bootstrap CI for calibration intercept + slope  [NEW]
  # ----------------------------------------------------------
  cal_boot <- boot::boot(
    data      = df_cal,
    statistic = function(d, i) {
      fit <- tryCatch(
        glm(y ~ logit_safe(p), data = d[i, ], family = binomial()),
        error   = function(e) NULL,
        warning = function(w) suppressWarnings(
          glm(y ~ logit_safe(p), data = d[i, ], family = binomial())
        )
      )
      if (is.null(fit)) return(c(NA_real_, NA_real_))
      coef(fit)
    },
    R = R_boot
  )
  
  # Percentile CIs — gracefully handle failures
  ci_intercept <- tryCatch(
    round(boot::boot.ci(cal_boot, index = 1, type = "perc")$percent[4:5], 3),
    error = function(e) c(NA_real_, NA_real_)
  )
  ci_slope <- tryCatch(
    round(boot::boot.ci(cal_boot, index = 2, type = "perc")$percent[4:5], 3),
    error = function(e) c(NA_real_, NA_real_)
  )
  
  # ----------------------------------------------------------
  # Hosmer-Lemeshow
  # ----------------------------------------------------------
  hl_p <- tryCatch(
    round(ResourceSelection::hoslem.test(truth, pred, g = 10)$p.value, 3),
    error = function(e) NA_real_
  )
  
  # ----------------------------------------------------------
  # Calibration in the large
  # ----------------------------------------------------------
  mean_pred <- round(mean(pred),  3)
  mean_obs  <- round(mean(truth), 3)
  cal_itl   <- round(mean_pred - mean_obs, 3)
  
  # ----------------------------------------------------------
  # Label permutation p
  # ----------------------------------------------------------
  perm_aucs <- replicate(R_boot, {
    as.numeric(pROC::auc(pROC::roc(sample(truth), pred, quiet = TRUE)))
  })
  perm_p <- round(mean(perm_aucs >= auc_val), 4)
  
  # ----------------------------------------------------------
  # Return — all original columns + new CI columns
  # ----------------------------------------------------------
  data.frame(
    model                = label,
    n                    = length(truth),
    events               = sum(truth),
    AUC                  = auc_val,
    AUC_lower95          = auc_ci[1],
    AUC_upper95          = auc_ci[3],
    Brier                = brier_val,
    Brier_lower95        = brier_ci[1],
    Brier_upper95        = brier_ci[2],
    Permutation_p        = perm_p,
    Cal_intercept        = intercept,
    Cal_intercept_lo95   = ci_intercept[1],   # NEW
    Cal_intercept_hi95   = ci_intercept[2],   # NEW
    Cal_slope            = slope,
    Cal_slope_lo95       = ci_slope[1],        # NEW
    Cal_slope_hi95       = ci_slope[2],        # NEW
    HL_p                 = hl_p,
    Mean_pred            = mean_pred,
    Observed_rate        = mean_obs,
    Cal_in_the_large     = cal_itl,
    stringsAsFactors     = FALSE
  )
}

# ============================================================
# The three run_calibration() calls below are UNCHANGED —
# they will now automatically produce the extra CI columns
# because the function signature is backwards-compatible.
# ============================================================

cat("\nRunning calibration (bootstrap + permutation ~2 min each)...\n")

cal_berlin <- run_calibration(
  truth = scores_prague$group01,
  pred  = scores_prague$pred_berlin_7,
  label = "Berlin 7-protein -> Prague"
)

cal_60_7 <- run_calibration(
  truth = scores_test40$group01,
  pred  = scores_test40$pred_train60_7,
  label = "Train60 7-protein -> Test40"
)

cal_60_conc <- run_calibration(
  truth = scores_test40_con$group01,
  pred  = scores_test40_con$pred_train60_conc,
  label = "Train60 Concentration -> Test40"
)

# Combine and save — same as before, now with extra columns
cal_all <- bind_rows(cal_berlin, cal_60_7, cal_60_conc)

print(cal_all)
row.names(cal_all) <- NULL
write_csv(cal_all, file.path(outdir, "calibration_all_models.csv"))

# ============================================================
# QUICK PRINT — formatted slope CIs for easy copy-paste
# into manuscript / response
# ============================================================
cat("\n=== Calibration slope summary (point estimate + 95% bootstrap CI) ===\n\n")
for (i in seq_len(nrow(cal_all))) {
  cat(sprintf("%-40s  slope = %.3f (95%% CI: %.3f – %.3f)\n",
              cal_all$model[i],
              cal_all$Cal_slope[i],
              cal_all$Cal_slope_lo95[i],
              cal_all$Cal_slope_hi95[i]))
  cat(sprintf("%-40s  intercept = %.3f (95%% CI: %.3f – %.3f)\n\n",
              "",
              cal_all$Cal_intercept[i],
              cal_all$Cal_intercept_lo95[i],
              cal_all$Cal_intercept_hi95[i]))
}
# ----------------------------------------------------------
# RUN CALIBRATION — GLM and Ridge on Prague
# ----------------------------------------------------------
cat("\nRunning calibration (bootstrap + permutation ~1 min each)...\n")

cal_glm   <- run_calibration(data_testP$group01, pred_glm,   "Berlin GLM -> Prague")
cal_ridge <- run_calibration(data_testP$group01, pred_ridge, "Berlin Ridge -> Prague")

cal_both <- bind_rows(cal_glm, cal_ridge)
print(cal_both)

row.names(cal_both) <- NULL

write_csv(cal_both, file.path(outdir, "calibration_7Berlin_naive_ridge.csv"))


# ----------------------------------------------------------
# COEFFICIENT COMPARISON TABLE
# ----------------------------------------------------------
coef_tab <- data.frame(
  predictor         = c("(Intercept)", prots_7),
  beta_GLM_berlin   = round(coef(model_berlin_7),  3),
  beta_Ridge        = round(as.numeric(coef(model_berlin_ridge)), 3),
  beta_GLM_train60  = round(coef(model_train60_7), 3),
  ratio_Berlin_Train60 = round(
    abs(coef(model_berlin_7)) / (abs(coef(model_train60_7)) + 0.001), 2
  ),
  stringsAsFactors = FALSE
)
print(coef_tab)

# ----------------------------------------------------------
# SAVE EVERYTHING
# ----------------------------------------------------------

# Model objects
saveRDS(model_berlin_ridge, file.path(artifact_dir, "model_7panel_Berlin_Ridge_glmnet.rds"))
saveRDS(cv_ridge,           file.path(artifact_dir, "model_7panel_Berlin_Ridge_cv_glmnet.rds"))
message("Saved: model_7panel_Berlin_Ridge_glmnet.rds")

# Frozen scores
scores_prague <- data_testP %>%
  select(Patient, cohort, group01) %>%
  mutate(pred_berlin_glm   = pred_glm,
         pred_berlin_ridge = pred_ridge)

write_csv(scores_prague,
          file.path(artifact_dir, "scores_Prague_GLM_and_Ridge_frozen.csv"))
message("Saved: scores_Prague_GLM_and_Ridge_frozen.csv")

# Coefficient table
write_csv(coef_tab,
          file.path(artifact_dir, "coef_Berlin_GLM_Ridge_Train60_comparison.csv"))
message("Saved: coef_Berlin_GLM_Ridge_Train60_comparison.csv")

# Calibration summary
write_csv(cal_both,
          file.path(outdir, "calibration_Berlin_GLM_vs_Ridge.csv"))
message("Saved: calibration_Berlin_GLM_vs_Ridge.csv")


# ============================================================
# EPV + Separation check
# ============================================================

cat("\n=== Events Per Variable (EPV) ===\n\n")

epv_check <- function(label, data, n_pred) {
  n1  <- sum(data$group01 == 1)
  n0  <- sum(data$group01 == 0)
  epv <- round(min(n1, n0) / n_pred, 2)
  cat(label, "\n")
  cat("  Remission (1):", n1, "| Active (0):", n0,
      "| Predictors:", n_pred, "| EPV:", epv, "\n\n")
  data.frame(training_set = label, remission = n1, active = n0,
             n_predictors = n_pred, EPV = epv, stringsAsFactors = FALSE)
}

epv_all <- bind_rows(
  epv_check("Berlin training set (7-protein)",   data_trainB,          7),
  epv_check("Train60 set (7-protein)",           train60,              7),
  epv_check("Train60 set (concentration model)", vsn_data_based_train, 7)
)

cat("=== Separation check — extreme predictions on validation sets ===\n\n")

sep_check <- function(label, pred, n_total) {
  lo <- sum(pred < 0.001)
  hi <- sum(pred > 0.999)
  cat(label, "\n")
  cat("  pred < 0.001:", lo, "| pred > 0.999:", hi, "| total n:", n_total, "\n\n")
  data.frame(model = label, pred_below_0001 = lo, pred_above_0999 = hi,
             total_n = n_total, stringsAsFactors = FALSE)
}

sep_all <- bind_rows(
  sep_check("Berlin GLM   -> Prague",  scores_prague$pred_berlin_glm,        nrow(scores_prague)),
  sep_check("Berlin Ridge -> Prague",  scores_prague$pred_berlin_ridge,       nrow(scores_prague)),
  sep_check("Train60 7p   -> Test40",  scores_test40$pred_train60_7,         nrow(scores_test40)),
  sep_check("Train60 Conc -> Test40",  scores_test40_con$pred_train60_conc,  nrow(scores_test40_con))
)

write_csv(epv_all, file.path(outdir, "epv_summary.csv"))
write_csv(sep_all, file.path(outdir, "separation_check.csv"))

# ============================================================
# Calibration plots (one SVG per model)
# ============================================================

library(ggplot2)

calibration_plot <- function(truth, pred, model_label, n_bins = 5, outdir = NULL) {
  
  truth <- as.integer(truth)
  pred  <- as.numeric(pred)
  
  df <- data.frame(truth = truth, pred = pred)
  
  cal_df <- df %>%
    dplyr::mutate(bin = dplyr::ntile(pred, n_bins)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(pred),
      obs_rate  = mean(truth),
      n         = dplyr::n(),
      obs_lo    = binom.test(sum(truth), dplyr::n())$conf.int[1],
      obs_hi    = binom.test(sum(truth), dplyr::n())$conf.int[2],
      .groups   = "drop"
    ) %>%
    dplyr::arrange(mean_pred)
  
  logit_safe <- function(p) {
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    log(p / (1 - p))
  }
  cal_fit   <- glm(truth ~ logit_safe(pred), data = df, family = binomial())
  intercept <- round(coef(cal_fit)[1], 3)
  slope     <- round(coef(cal_fit)[2], 3)
  brier     <- round(mean((pred - truth)^2), 3)
  hl_p      <- tryCatch(
    round(ResourceSelection::hoslem.test(truth, pred, g = 10)$p.value, 3),
    error = function(e) NA_real_
  )
  
  stats_label <- paste0(
    "Intercept = ", intercept, "\n",
    "Slope     = ", slope,     "\n",
    "Brier     = ", brier,     "\n",
    "H-L p     = ", hl_p
  )
  
  p <- ggplot(df, aes(x = pred, y = truth)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    geom_point(data = cal_df,
               aes(x = mean_pred, y = obs_rate),
               size = 3, inherit.aes = FALSE) +
    geom_errorbar(data = cal_df,
                  aes(x = mean_pred, ymin = obs_lo, ymax = obs_hi),
                  width = 0.02, inherit.aes = FALSE) +
    geom_line(data = cal_df,
              aes(x = mean_pred, y = obs_rate),
              inherit.aes = FALSE) +
    annotate("text", x = 0.02, y = 0.98,
             label = stats_label,
             hjust = 0, vjust = 1, size = 3.2, family = "mono") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 12) +
    labs(title = paste0("Calibration: ", model_label),
         x = "Predicted probability",
         y = "Observed event rate (binwise)")
  
  if (!is.null(outdir)) {
    out_svg <- file.path(outdir,
                         paste0("calibration_", gsub("[ /->]", "_", model_label), ".svg"))
    ggsave(filename = out_svg, plot = p, width = 5, height = 5)
    message("Saved: ", out_svg)
  }
  
  invisible(p)
}

# --- Run for all four models ---
calibration_plot(scores_prague$group01,         scores_prague$pred_berlin_glm,
                 "Berlin GLM -> Prague",         outdir = outdir)

calibration_plot(scores_prague$group01,        scores_prague$pred_berlin_ridge,
                 "Berlin Ridge -> Prague",       outdir = outdir)

calibration_plot(scores_test40$group01,         scores_test40$pred_train60_7,
                 "Train60 7-protein -> Test40",  outdir = outdir)

calibration_plot(scores_test40_con$group01,     scores_test40_con$pred_train60_conc,
                 "Train60 Conc -> Test40",       outdir = outdir)

# Bootstrap CI for calibration intercept & slope
cal_boot_ci <- function(truth, pred, R = 2000, seed = 7) {
  set.seed(seed)
  n <- length(truth)
  df <- data.frame(y = as.integer(truth),
                   p = pmin(pmax(pred, .Machine$double.eps), 1 - .Machine$double.eps))
  logit_safe <- function(p) log(p / (1 - p))
  stat <- function(data, i) {
    d <- data[i, ]
    fit <- glm(y ~ logit_safe(p), data = d, family = binomial())
    coef(fit)
  }
  b <- boot::boot(data = df, statistic = function(d, i) stat(d, i), R = R)
  ci_intercept <- tryCatch(boot::boot.ci(b, index = 1, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  ci_slope     <- tryCatch(boot::boot.ci(b, index = 2, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  list(intercept = c(point = b$t0[1], lower = ci_intercept[1], upper = ci_intercept[2]),
       slope     = c(point = b$t0[2], lower = ci_slope[1],     upper = ci_slope[2]),
       boot_obj   = b)
}
# Example usage:
boot_cal_glm   <- cal_boot_ci(data_testP$group01, pred_glm,   R = 2000)
boot_cal_ridge <- cal_boot_ci(data_testP$group01, pred_ridge, R = 2000)


# Score the Berlin training set with both models (GLM and Ridge)
pred_train_glm   <- as.numeric(predict(model_berlin_7, newdata = data_trainB, type = "response"))
pred_train_ridge <- as.numeric(predict(model_berlin_ridge, newx = as.matrix(data_trainB[, prots_7]), type = "response"))

# Compute separation diagnostics for both
sep_train_glm <- data_trainB %>%
  summarise(
    model = "Berlin GLM (train)",
    total_n = n(),
    pred_below_0001 = sum(pred_train_glm < 0.001),
    pred_above_0999 = sum(pred_train_glm > 0.999),
    n_extreme = pred_below_0001 + pred_above_0999,
    pct_extreme = round(100 * n_extreme / total_n, 1)
  )

sep_train_ridge <- data_trainB %>%
  summarise(
    model = "Berlin Ridge (train)",
    total_n = n(),
    pred_below_0001 = sum(pred_train_ridge < 0.001),
    pred_above_0999 = sum(pred_train_ridge > 0.999),
    n_extreme = pred_below_0001 + pred_above_0999,
    pct_extreme = round(100 * n_extreme / total_n, 1)
  )

sep_train_both <- bind_rows(sep_train_glm, sep_train_ridge)

# Print and save
print(sep_train_both)
write_csv(sep_train_both, file.path(outdir, "separation_berlin_train_by_model.csv"))



