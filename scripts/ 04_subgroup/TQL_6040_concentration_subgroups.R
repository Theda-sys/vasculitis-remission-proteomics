# ============================================================
# TQL 60/40 — CONCENTRATION MODEL: METRICS + SCATTER PLOT
# ============================================================
# Loads the FROZEN saved model and frozen test scores (or falls
# back to re-scoring if the CSV is present but the model is
# needed for the train-set scatter).
# ============================================================

library(pROC)
library(dplyr)
library(readr)
library(boot)
library(ggplot2)

set.seed(7)

# ------------------------------------------------------------------
# PATHS — adjust to match your directory layout
# ------------------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final"
outdir       <- "../proteomics_ml/data/TQL_60_40_concentration/finalx"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Frozen score CSVs written during model training
scores_test40_path  <- file.path(artifact_dir, "scores_test40_frozen_concentration.csv")
scores_train60_path <- file.path(artifact_dir, "scores_train60_concentration_frozen.csv")

# Saved model (used only for the train-set Youden threshold)
model_conc_path <- file.path(artifact_dir, "model_conc_Train60_glm.rds")

# ------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------
logit_safe <- function(p) {
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  log(p / (1 - p))
}

binom_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int[1:2]
}

# Full evaluation: AUC (DeLong + bootstrap), Brier, threshold metrics,
# calibration slope/intercept — returns a one-row data.frame
eval_model_conc <- function(truth01, pred_prob, model_name,
                            outdir, R_boot = 2000, threshold = 0.5) {
  
  truth01   <- as.integer(truth01)
  pred_prob <- as.numeric(pred_prob)
  df_boot   <- data.frame(y = truth01, p = pred_prob)
  
  # ---- AUC ----------------------------------------------------------------
  roc_obj <- pROC::roc(response = truth01, predictor = pred_prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # flip if necessary (very unlikely for a well-trained model)
  if (auc_val < 0.5) {
    warning(model_name, ": AUC < 0.5 — check label coding!")
    roc_obj <- pROC::roc(response = truth01, predictor = 1 - pred_prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
  }
  
  auc_ci_delong <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  
  auc_stat <- function(d, i) {
    di <- d[i, , drop = FALSE]
    r  <- tryCatch(pROC::roc(response = di$y, predictor = di$p, quiet = TRUE),
                   error = function(e) NULL)
    if (is.null(r)) return(NA_real_)
    as.numeric(pROC::auc(r))
  }
  boot_auc <- { set.seed(7); bt <- boot::boot(df_boot, auc_stat, R = R_boot);
  list(ci = quantile(bt$t, c(0.025, 0.975), na.rm = TRUE)) }
  
  # ---- Brier ---------------------------------------------------------------
  brier_val  <- mean((pred_prob - truth01)^2)
  brier_stat <- function(d, i) mean((d$p[i] - d$y[i])^2)
  boot_brier <- { set.seed(8); bt <- boot::boot(df_boot, brier_stat, R = R_boot);
  list(ci = quantile(bt$t, c(0.025, 0.975), na.rm = TRUE)) }
  
  # ---- Threshold metrics at 0.5 -------------------------------------------
  pred_bin <- as.integer(pred_prob >= threshold)
  TP <- sum(pred_bin == 1 & truth01 == 1, na.rm = TRUE)
  TN <- sum(pred_bin == 0 & truth01 == 0, na.rm = TRUE)
  FP <- sum(pred_bin == 1 & truth01 == 0, na.rm = TRUE)
  FN <- sum(pred_bin == 0 & truth01 == 1, na.rm = TRUE)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  
  sens_ci <- if ((TP + FN) > 0) binom_ci(TP, TP + FN) else c(NA, NA)
  spec_ci <- if ((TN + FP) > 0) binom_ci(TN, TN + FP) else c(NA, NA)
  ppv_ci  <- if ((TP + FP) > 0) binom_ci(TP, TP + FP) else c(NA, NA)
  npv_ci  <- if ((TN + FN) > 0) binom_ci(TN, TN + FN) else c(NA, NA)
  
  # ---- Youden threshold ----------------------------------------------------
  roc_coords    <- pROC::coords(roc_obj, "best", best.method = "youden",
                                ret = c("threshold", "sensitivity", "specificity"))
  youden_thresh <- as.numeric(roc_coords["threshold"])
  sens_youden   <- as.numeric(roc_coords["sensitivity"])
  spec_youden   <- as.numeric(roc_coords["specificity"])
  
  # ---- Calibration intercept & slope (bootstrap CIs) ----------------------
  cal_fit       <- glm(truth01 ~ logit_safe(pred_prob), family = binomial())
  cal_intercept <- as.numeric(coef(cal_fit)[1])
  cal_slope     <- as.numeric(coef(cal_fit)[2])
  
  cal_stat <- function(d, i) {
    di  <- d[i, , drop = FALSE]
    fit <- tryCatch(
      withCallingHandlers(
        glm(di$y ~ logit_safe(di$p), family = binomial()),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    as.numeric(coef(fit))
  }
  boot_cal  <- { set.seed(9); boot::boot(df_boot, cal_stat, R = R_boot) }
  ci_cal    <- t(apply(boot_cal$t, 2, function(x)
    if (all(is.na(x))) c(NA, NA) else quantile(x, c(0.025, 0.975), na.rm = TRUE)))
  
  # ---- Cal-in-the-large & Hosmer-Lemeshow ----------------------------------
  citl   <- mean(pred_prob) - mean(truth01)
  hl_p   <- tryCatch(
    ResourceSelection::hoslem.test(truth01, pred_prob, g = 10)$p.value,
    error = function(e) NA_real_
  )
  
  # ---- Assemble one-row summary --------------------------------------------
  one_row <- data.frame(
    model           = model_name,
    n_test          = length(truth01),
    events_test     = sum(truth01 == 1, na.rm = TRUE),
    auc             = auc_val,
    auc_lo_delong   = auc_ci_delong[1],
    auc_hi_delong   = auc_ci_delong[3],
    auc_lo_boot     = as.numeric(boot_auc$ci[1]),
    auc_hi_boot     = as.numeric(boot_auc$ci[2]),
    brier           = brier_val,
    brier_lo_boot   = as.numeric(boot_brier$ci[1]),
    brier_hi_boot   = as.numeric(boot_brier$ci[2]),
    sens_05         = sens,  sens_05_lo = sens_ci[1], sens_05_hi = sens_ci[2],
    spec_05         = spec,  spec_05_lo = spec_ci[1], spec_05_hi = spec_ci[2],
    ppv_05          = ppv,   ppv_05_lo  = ppv_ci[1],  ppv_05_hi  = ppv_ci[2],
    npv_05          = npv,   npv_05_lo  = npv_ci[1],  npv_05_hi  = npv_ci[2],
    youden_thresh   = youden_thresh,
    sens_youden     = sens_youden,
    spec_youden     = spec_youden,
    cal_intercept   = cal_intercept,
    cal_intercept_lo = ci_cal[1, 1], cal_intercept_hi = ci_cal[1, 2],
    cal_slope       = cal_slope,
    cal_slope_lo    = ci_cal[2, 1], cal_slope_hi = ci_cal[2, 2],
    citl            = citl,
    hl_p            = hl_p,
    TP = TP, TN = TN, FP = FP, FN = FN,
    stringsAsFactors = FALSE
  )
  
  # save individual CSV
  write.csv(one_row,
            file = file.path(outdir, paste0("summary_one_row_", model_name, ".csv")),
            row.names = FALSE)
  
  # save diagnostics RDS (for later use / DeLong comparisons)
  saveRDS(list(one_row  = one_row,
               roc      = roc_obj,
               boot_auc = boot_auc,
               boot_brier = boot_brier,
               boot_cal = boot_cal),
          file = file.path(outdir, paste0("diagnostics_", model_name, ".rds")))
  
  message(sprintf(
    "[%s]  AUC=%.3f (%.3f–%.3f)  |  Brier=%.3f  |  Sens=%.3f  Spec=%.3f  PPV=%.3f  NPV=%.3f",
    model_name, auc_val, auc_ci_delong[1], auc_ci_delong[3],
    brier_val, sens, spec, ppv, npv))
  
  return(list(one_row = one_row, roc_obj = roc_obj))
}

# ------------------------------------------------------------------
# SCATTER PLOT HELPER (same style as TQL_models.R)
# ------------------------------------------------------------------
make_scatter_plot <- function(scores_df, pred_col, group_col,
                              title_str, outdir, filename,
                              threshold = 0.5) {
  
  stopifnot(pred_col  %in% colnames(scores_df))
  stopifnot(group_col %in% colnames(scores_df))
  
  set.seed(7)
  plot_df <- scores_df %>%
    mutate(
      risk_score = as.numeric(.data[[pred_col]]),
      outcome    = factor(.data[[group_col]],
                          levels = c(0, 1),
                          labels = c("Active", "Remission")),
      y_jitter   = runif(n(), -0.08, 0.08)
    )
  
  p <- ggplot(plot_df, aes(x = risk_score, y = y_jitter)) +
    geom_hline(yintercept = 0, colour = "gray80") +
    
    # shape AND fill both driven by outcome — no hardcoded pch override
    geom_point(aes(fill = outcome, shape = outcome),
               size = 3, stroke = 0.4, colour = "black") +
    
    # fixed 0.5 threshold
    geom_vline(xintercept = threshold,
               linetype = "dashed", color = "black", linewidth = 0.6) +
    annotate("text",
             x = threshold, y = 0.18,
             label = paste0("Threshold = ", threshold),
             angle = 90, vjust = -0.5, hjust = 0.5, size = 3.5) +
    
    geom_rug(aes(x = risk_score), sides = "b", alpha = 0.3) +
    
    # 21 = filled circle, 24 = filled upward triangle (both respect fill + colour)
    scale_shape_manual(values = c("Active"    = 21,
                                  "Remission" = 24)) +
    scale_fill_manual(values  = c("Active"    = "#ffd973",
                                  "Remission" = "#0000ab")) +
    
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.25, 0.25)) +
    labs(x     = "Predicted probability of remission",
         y     = NULL,
         title = title_str,
         fill  = "True outcome",
         shape = "True outcome") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  out_path <- file.path(outdir, filename)
  ggsave(filename = out_path, plot = p, width = 8, height = 3.0)
  message("Saved scatter: ", out_path)
  invisible(p)
}

# ------------------------------------------------------------------
# ROC PLOT HELPER
# ------------------------------------------------------------------
make_roc_plot <- function(roc_obj, model_label, outdir, filename) {
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  
  out_path <- file.path(outdir, filename)
  # use svg if svglite is available, otherwise png
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(out_path, width = 5.5, height = 5.2)
  } else {
    out_path <- sub("\\.svg$", ".png", out_path)
    grDevices::png(out_path, width = 550, height = 520)
  }
  
  plot(roc_obj, col = "#c70000", lwd = 2.2,
       xlab = "1 – Specificity", ylab = "Sensitivity",
       main = sprintf("ROC — Concentration model (Train60 → Test40)\nn=%d, events=%d",
                      length(roc_obj$response),
                      sum(roc_obj$response == 1)))
  abline(a = 0, b = 1, lty = 3, col = "grey60")
  legend("bottomright",
         legend = sprintf("%s\nAUC = %.3f (%.3f–%.3f)",
                          model_label, auc_val, auc_ci[1], auc_ci[3]),
         col = "#c70000", lwd = 2.2, bty = "n", cex = 0.82)
  dev.off()
  message("Saved ROC plot: ", out_path)
}

# ==================================================================
# MAIN EXECUTION
# ==================================================================

# --- 1) Load frozen TEST scores -------------------------------------------
if (!file.exists(scores_test40_path)) {
  stop("Frozen test scores not found:\n  ", scores_test40_path,
       "\nPlease ensure scores_test40_frozen_concentration.csv exists.")
}
scores_test40_con <- read_csv(scores_test40_path, show_col_types = FALSE)

# Sanity checks
stopifnot("group01"          %in% colnames(scores_test40_con))
stopifnot("pred_train60_conc" %in% colnames(scores_test40_con))

n_test    <- nrow(scores_test40_con)
n_events  <- sum(scores_test40_con$group01 == 1, na.rm = TRUE)
obs_rate  <- mean(scores_test40_con$group01,    na.rm = TRUE)
mean_pred <- mean(scores_test40_con$pred_train60_conc, na.rm = TRUE)

message(sprintf(
  "Test40 loaded: n=%d | events (remission)=%d | obs rate=%.3f | mean pred=%.3f",
  n_test, n_events, obs_rate, mean_pred))

# Plausibility check — flag if predictions look reversed or collapsed
if (mean_pred < 0.1 || mean_pred > 0.9)
  warning("Mean predicted probability is extreme (", round(mean_pred, 3),
          "). Check label / prediction coding!")
if (abs(mean_pred - obs_rate) > 0.15)
  warning("Cal-in-the-large gap is large: mean_pred=", round(mean_pred, 3),
          " vs obs_rate=", round(obs_rate, 3), ". Consider recalibration.")

# --- 2) Evaluate on test set ----------------------------------------------
res_conc <- eval_model_conc(
  truth01    = scores_test40_con$group01,
  pred_prob  = scores_test40_con$pred_train60_conc,
  model_name = "Train60_Conc_Test40",
  outdir     = outdir,
  R_boot     = 2000,
  threshold  = 0.5
)

# Save rounded summary table
summary_rounded <- res_conc$one_row
num_cols <- sapply(summary_rounded, is.numeric)
summary_rounded[, num_cols] <- lapply(summary_rounded[, num_cols],
                                      function(x) round(x, 3))
write.csv(summary_rounded,
          file = file.path(outdir, "supp_table_conc_perform.csv"),
          row.names = FALSE)
message("Saved performance summary: supp_table_conc_perform.csv")

# Pretty print to console
with(res_conc$one_row, {
  cat(sprintf("  n = %d | events (remission) = %d\n", n_test, events_test))
  cat(sprintf("  AUC     : %.3f (DeLong 95%% CI: %.3f – %.3f)\n",
              auc, auc_lo_delong, auc_hi_delong))
  cat(sprintf("  Brier   : %.3f (95%% CI: %.3f – %.3f)\n",
              brier, brier_lo_boot, brier_hi_boot))
  cat(sprintf("  Sens    : %.3f (%.3f – %.3f)  @  0.5 threshold\n",
              sens_05, sens_05_lo, sens_05_hi))
  cat(sprintf("  Spec    : %.3f (%.3f – %.3f)  @  0.5 threshold\n",
              spec_05, spec_05_lo, spec_05_hi))
  cat(sprintf("  PPV     : %.3f (%.3f – %.3f)\n", ppv_05, ppv_05_lo, ppv_05_hi))
  cat(sprintf("  NPV     : %.3f (%.3f – %.3f)\n", npv_05, npv_05_lo, npv_05_hi))
  cat(sprintf("  TP=%d  TN=%d  FP=%d  FN=%d\n", TP, TN, FP, FN))
  cat(sprintf("  Youden threshold: %.3f  (sens=%.3f, spec=%.3f)\n",
              youden_thresh, sens_youden, spec_youden))
  cat(sprintf("  Cal intercept: %.3f (%.3f – %.3f)\n",
              cal_intercept, cal_intercept_lo, cal_intercept_hi))
  cat(sprintf("  Cal slope    : %.3f (%.3f – %.3f)\n",
              cal_slope, cal_slope_lo, cal_slope_hi))
  cat(sprintf("  CITL : %.4f\n", citl))
  cat(sprintf("  H-L p: %.3f\n", hl_p))
})

# --- 3) Scatter plot — TEST set -------------------------------------------
make_scatter_plot(
  scores_df  = scores_test40_con,
  pred_col   = "pred_train60_conc",
  group_col  = "group01",
  title_str  = "Concentration model — predicted probability (test set, 40%)",
  outdir     = outdir,
  filename   = "scatter_conc_test40_outcome.svg"
)

# --- 4) Scatter plot — TRAIN set (requires frozen train scores CSV) --------
if (file.exists(scores_train60_path)) {
  scores_train60_con <- read_csv(scores_train60_path, show_col_types = FALSE)
  stopifnot("group01"           %in% colnames(scores_train60_con))
  stopifnot("pred_train60_conc" %in% colnames(scores_train60_con))
  
  make_scatter_plot(
    scores_df  = scores_train60_con,
    pred_col   = "pred_train60_conc",
    group_col  = "group01",
    title_str  = "Concentration model — predicted probability (training set, 60%)",
    outdir     = outdir,
    filename   = "scatter_conc_train60_outcome.svg"
  )
  
  # Quick sanity: apparent AUC on training set (expected to be higher — check for overfitting)
  roc_train   <- pROC::roc(response  = as.integer(scores_train60_con$group01),
                           predictor = scores_train60_con$pred_train60_conc,
                           quiet = TRUE)
  auc_train   <- as.numeric(pROC::auc(roc_train))
  auc_test    <- res_conc$one_row$auc
  message(sprintf("Apparent train AUC=%.3f vs test AUC=%.3f  (optimism=%.3f)",
                  auc_train, auc_test, auc_train - auc_test))
  if ((auc_train - auc_test) > 0.10)
    warning("Optimism > 0.10 — consider bootstrap or cross-validation correction.")
  
} else {
  message("Train60 scores not found at: ", scores_train60_path,
          "\nSkipping train-set scatter. Create it by adding the following to your training script:\n",
          "  write_csv(vsn_data_based_train %>% select(Patient, cohort, group01, pred_train60_conc),\n",
          "            file.path(outdir, 'scores_train60_concentration_frozen.csv'))")
}

# --- 5) ROC plot ----------------------------------------------------------
make_roc_plot(
  roc_obj     = res_conc$roc_obj,
  model_label = "7-protein concentration model",
  outdir      = outdir,
  filename    = "ROC_conc_train60_test40.svg"
)



# ==================================================================
# ENHANCED CONFUSION MATRIX PLOT ("Uwe-Plot" style)
# ==================================================================
# Actual class on x-axis, predicted class on y-axis.
# Points jittered within each quadrant; quadrant shading makes
# TP/TN (correct) vs FP/FN (errors) immediately visible.
# Marginal counts annotated in each quadrant.
# ==================================================================

make_confusion_scatter <- function(scores_df, pred_col, group_col,
                                   outdir, filename,
                                   threshold = 0.5,
                                   title_str = "Classification result (test set)") {
  
  library(ggplot2)
  library(dplyr)
  
  truth01   <- as.integer(scores_df[[group_col]])
  pred_prob <- as.numeric(scores_df[[pred_col]])
  pred_bin  <- as.integer(pred_prob >= threshold)
  
  set.seed(7)
  plot_df <- data.frame(
    Actual    = truth01,
    Predicted = pred_bin,
    Prob      = pred_prob
  ) %>%
    mutate(
      x_jit = Actual    + runif(n(), -0.28, 0.28),
      y_jit = Predicted + runif(n(), -0.28, 0.28),
      quad_label = case_when(
        Actual == 1 & Predicted == 1 ~ "TP",
        Actual == 0 & Predicted == 0 ~ "TN",
        Actual == 0 & Predicted == 1 ~ "FP",
        Actual == 1 & Predicted == 0 ~ "FN"
      ),
      # fill = PREDICTED class (yellow = predicted Active, blue = predicted Remission)
      fill_class = factor(ifelse(Predicted == 1, "Remission", "Active"),
                          levels = c("Active", "Remission")),
      # shape = TRUE class (circle 21 = Active, triangle 24 = Remission)
      shape_class = factor(ifelse(Actual == 1, "Remission", "Active"),
                           levels = c("Active", "Remission"))
    )
  
  # quadrant count labels
  counts <- plot_df %>%
    group_by(quad_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      x     = case_when(quad_label %in% c("TP", "FN") ~ 1,  TRUE ~ 0),
      y     = case_when(quad_label %in% c("TP", "FP") ~ 1,  TRUE ~ 0),
      label = paste0(quad_label, "\nn = ", n)
    )
  
  # background shading: predicted class colour tint
  # bottom row (predicted Active) = yellow tint, top row (predicted Remission) = blue tint
  shading <- data.frame(
    xmin = c(-0.5, -0.5,  0.5,  0.5),
    xmax = c( 0.5,  0.5,  1.5,  1.5),
    ymin = c(-0.5,  0.5, -0.5,  0.5),
    ymax = c( 0.5,  1.5,  0.5,  1.5),
    fill = c("#FFFFC5", "#7373FF", "#7373FF", "#FFFFC5")  # TN=yellow, FN=blue, FP=blue, TP=yellow
  )
  
  p <- ggplot(plot_df, aes(x = x_jit, y = y_jit)) +
    
    # quadrant backgrounds
    geom_rect(data = shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = shading$fill, inherit.aes = FALSE, alpha = 0.5) +
    
    # dividing lines
    geom_hline(yintercept = 0.5, colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 0.5, colour = "grey40", linewidth = 0.5) +
    
    # points: fill = predicted class, shape = true class
    geom_point(aes(fill = fill_class, shape = shape_class),
               size = 3.8, stroke = 0.5, colour = "black", alpha = 0.88) +
    
    # fill scale: predicted class colours
    scale_fill_manual(
      name   = "Predicted class",
      values = c("Active" = "#ffd973", "Remission" = "#0000ab"),
      guide  = guide_legend(override.aes = list(shape = 21, size = 4))
    ) +
    
    # shape scale: true class — circle for Active, filled triangle for Remission
    scale_shape_manual(
      name   = "True class",
      values = c("Active" = 21, "Remission" = 24),
      guide  = guide_legend(override.aes = list(fill = "grey60", size = 4))
    ) +
    
    # quadrant count labels
    geom_label(data = counts,
               aes(x = x, y = y, label = label),
               inherit.aes = FALSE,
               size = 3.8, fontface = "bold",
               fill = "white", label.size = 0.3, alpha = 0.88,
               vjust = "inward", hjust = "inward") +
    
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Active\n(true)", "Remission\n(true)"),
      limits = c(-0.5, 1.5)
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("Active\n(predicted)", "Remission\n(predicted)"),
      limits = c(-0.5, 1.5)
    ) +
    
    labs(
      title    = title_str,
      subtitle = sprintf("Threshold = %.1f  |  n = %d  |  AUC = %.3f",
                         threshold, nrow(plot_df),
                         as.numeric(pROC::auc(
                           pROC::roc(truth01, pred_prob, quiet = TRUE)))),
      x = "Actual class",
      y = "Predicted class"
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey35"),
      axis.text     = element_text(size = 11),
      axis.title    = element_text(size = 12),
      panel.grid    = element_blank(),
      legend.position = "right",
      legend.title  = element_text(size = 10, face = "bold"),
      legend.text   = element_text(size = 9)
    )
  
  out_path <- file.path(outdir, filename)
  ggsave(filename = out_path, plot = p, width = 5.5, height = 5.0)
  message("Saved confusion scatter: ", out_path)
  invisible(p)
}

# --- Run confusion matrix plot on test set -----------------------------------
# install.packages("ggnewscale")  # run once if not yet installed
make_confusion_scatter(
  scores_df  = scores_test40_con,
  pred_col   = "pred_train60_conc",
  group_col  = "group01",
  outdir     = outdir,
  filename   = "confusion_scatter_conc_test40.svg",
  threshold  = 0.5,
  title_str  = "Concentration model — test set (40%)"
)


# ============================================================
# MODEL COMPARISON — concentration vs clinical markers
# ============================================================

# Load full dataset to get train60 with clinical markers
data_csv <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
train60_file <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test40_file  <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"

dfz <- readr::read_csv(data_csv, show_col_types = FALSE)
load(train60_file)  # full_meta_train
load(test40_file)   # full_meta_test

train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)

# Attach frozen concentration scores to test40
# (scores_test40_con already has group01 and pred_train60_conc)
test40 <- test40 %>%
  left_join(scores_test40_con %>% select(Patient, pred_train60_conc),
            by = "Patient")

# Sanity check
stopifnot(sum(!is.na(test40$pred_train60_conc)) > 0)
cat("Matched concentration scores:", sum(!is.na(test40$pred_train60_conc)),
    "of", nrow(test40), "test patients\n")

# ---- helpers ---------------------------------------------------------------
binom_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int[1:2]
}

compute_metrics <- function(truth, prob, model_name, threshold = 0.5) {
  truth <- as.integer(truth)
  prob  <- as.numeric(prob)
  complete <- !is.na(truth) & !is.na(prob)
  truth <- truth[complete]; prob <- prob[complete]
  n      <- length(truth)
  events <- sum(truth == 1)
  non_ev <- sum(truth == 0)
  if (events > 0 && non_ev > 0) {
    roc_obj <- pROC::roc(truth, prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  } else {
    roc_obj <- NULL; auc_val <- NA_real_; auc_ci <- rep(NA_real_, 3)
  }
  pred_bin <- as.integer(prob >= threshold)
  TP <- sum(pred_bin == 1 & truth == 1)
  TN <- sum(pred_bin == 0 & truth == 0)
  FP <- sum(pred_bin == 1 & truth == 0)
  FN <- sum(pred_bin == 0 & truth == 1)
  sens <- if ((TP+FN) > 0) TP/(TP+FN) else NA_real_
  spec <- if ((TN+FP) > 0) TN/(TN+FP) else NA_real_
  ppv  <- if ((TP+FP) > 0) TP/(TP+FP) else NA_real_
  npv  <- if ((TN+FN) > 0) TN/(TN+FN) else NA_real_
  sens_ci <- if ((TP+FN) > 0) binom.test(TP, TP+FN)$conf.int else c(NA, NA)
  spec_ci <- if ((TN+FP) > 0) binom.test(TN, TN+FP)$conf.int else c(NA, NA)
  ppv_ci  <- if ((TP+FP) > 0) binom.test(TP, TP+FP)$conf.int else c(NA, NA)
  npv_ci  <- if ((TN+FN) > 0) binom.test(TN, TN+FN)$conf.int else c(NA, NA)
  list(
    roc_obj = roc_obj,
    table = data.frame(
      model = model_name, n = n, events = events,
      auc = auc_val, auc_lo = auc_ci[1], auc_hi = auc_ci[3],
      sens = sens, sens_lo = sens_ci[1], sens_hi = sens_ci[2],
      spec = spec, spec_lo = spec_ci[1], spec_hi = spec_ci[2],
      ppv  = ppv,  ppv_lo  = ppv_ci[1],  ppv_hi  = ppv_ci[2],
      npv  = npv,  npv_lo  = npv_ci[1],  npv_hi  = npv_ci[2],
      TP = TP, TN = TN, FP = FP, FN = FN,
      stringsAsFactors = FALSE
    )
  )
}

# ---- fit augmented models on train60 ---------------------------------------
# concentration score as single predictor, augmented with CRP / ANCA
# Note: train60 needs pred_train60_conc for train-set scoring
train60 <- train60 %>%
  left_join(
    readr::read_csv(scores_train60_path, show_col_types = FALSE) %>%
      select(Patient, pred_train60_conc),
    by = "Patient"
  )

mod_conc      <- glm(group01 ~ pred_train60_conc,                           data = train60, family = binomial())

# ---- score test40 ----------------------------------------------------------
test40$prob_conc      <- predict(mod_conc,       newdata = test40, type = "response")


# ---- compute metrics -------------------------------------------------------
res_conc_main  <- compute_metrics(test40$group01, test40$prob_conc,      "7-protein concentration")

metrics_conc_all <- do.call(rbind, lapply(
  list(res_conc_main),
  function(x) x$table
))
rownames(metrics_conc_all) <- NULL
num_cols <- sapply(metrics_conc_all, is.numeric)
metrics_conc_all[, num_cols] <- lapply(metrics_conc_all[, num_cols], round, 3)

write.csv(metrics_conc_all,
          file.path(outdir, "conc_model_comparison_metrics.csv"),
          row.names = FALSE)


