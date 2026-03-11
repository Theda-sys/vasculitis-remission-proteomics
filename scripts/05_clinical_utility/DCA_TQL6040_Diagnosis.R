# ============================================================
# DECISION CURVE ANALYSIS — TQL 60/40 — v7 (FINAL)
# ============================================================
# Clinical framing:
#   Outcome event = REMISSION (group01 == 1).
#   Clinical question: "Is this patient in remission and safe
#   to taper/stop immunosuppression?"
#   False positive = tapering a patient still in active disease
#   = relapse risk (the harm penalised by the NB formula).
# ============================================================

library(dcurves)
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(readr)
library(boot)

set.seed(7)

# ============================================================
# PATHS
# ============================================================
data_csv          <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
model_dir         <- "../proteomics_ml/data/models_and_artifacts/final"
frozen_scores_csv <- file.path(model_dir, "scores_test40_frozen.csv")
cal_output_csv    <- "../proteomics_ml/data/calibration_output/final/calibration_all_models.csv"
train60_file      <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test40_file       <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"
outdir            <- "../proteomics_ml/data/TQL_60_40_DCA/finalx"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# TUNABLE SETTINGS — all defined here, never redefined below
# ============================================================
N_BOOT_DCA       <- 500   # dcurves CI bands (increase to 2000 for final figures)
B_DIFF           <- 1000  # bootstrap reps for manual NB difference CI
B_BOOT           <- 1000  # bootstrap reps for cohort/weighted NB tables
clin_thresh      <- c(0.20, 0.25, 0.30, 0.40)
annot_thresholds <- c(0.20, 0.30, 0.40)
THRESH_FULL      <- seq(0.05, 0.50, by = 0.01)
THRESH_CLINICAL  <- seq(0.15, 0.50, by = 0.01)
THRESH_REPORT    <- c(0.15, 0.20, 0.25, 0.30, 0.40)

model_cols <- c(
  "7-protein panel"        = "#C70500",
  "CRP alone"              = "#89CC48",
  "ANCA alone"             = "#66CCD1",
  "ANCA + CRP"   = "#A260D8",
  "7-protein + CRP"        = "#888888",
  "7-protein + ANCA"       = "#aaaaaa",
  "7-protein + CRP + ANCA" = "#BBC2D5",
  "All"                    = "#D5CD78",
  "None"                   = "black"
)

# ============================================================
# HELPERS — defined once
# ============================================================

# DCA wrapper — detects dcurves version at runtime.
# dcurves >= 0.4.0: bootstrapped CI bands in figures.
# dcurves <  0.4.0: runs without CI bands, prints clear message.
# Quantified uncertainty always available via manual bootstrap
# difference plot regardless of dcurves version.
dca_with_ci <- function(formula, data, thresholds, label, n_boot = N_BOOT_DCA) {
  has_boot_arg <- "bootstraps" %in% names(formals(dcurves::dca))
  if (has_boot_arg) {
    dca(formula = formula, data = data, thresholds = thresholds,
        label = label, bootstraps = n_boot)
  } else {
    message("dcurves < 0.4.0: running without CI bands in figures. ",
            "Update dcurves to >= 0.4.0 to enable. ",
            "Manual bootstrap difference plot provides quantified uncertainty.")
    dca(formula = formula, data = data, thresholds = thresholds, label = label)
  }
}

# Robust coercion of dca object to tibble
dca_to_tibble <- function(dca_obj) {
  tryCatch(as_tibble(dca_obj), error = function(e) as_tibble(as.data.frame(dca_obj)))
}

# ggplot2 theme for all figures
dca_theme_nature <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = base_size),
      plot.title       = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1)
    )
}

# Net benefit (optionally weighted).
# obs = 1 means REMISSION (event of interest throughout this script).
nb_single <- function(obs, pred, thresh, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(obs))
  idx_tp <- (pred >= thresh) & (obs == 1)
  idx_fp <- (pred >= thresh) & (obs == 0)
  TP_w   <- sum(weights[idx_tp])
  FP_w   <- sum(weights[idx_fp])
  N_w    <- sum(weights)
  (TP_w / N_w) - (FP_w / N_w) * (thresh / (1 - thresh))
}

# Bootstrap CI for a single NB estimate
nb_boot_ci <- function(obs, pred, thresh, weights = NULL, B = B_BOOT) {
  df_boot <- data.frame(obs = obs, pred = pred,
                        w = if (is.null(weights)) 1 else weights)
  boot_stat <- function(data, i) {
    nb_single(as.numeric(data[i, "obs"]), as.numeric(data[i, "pred"]), thresh,
              weights = if (all(data[i, "w"] == 1)) NULL else as.numeric(data[i, "w"]))
  }
  b  <- boot::boot(data = df_boot, statistic = boot_stat, R = B)
  ci <- tryCatch(boot::boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  list(nb = b$t0, ci_low = ci[1], ci_high = ci[2])
}

# Paired bootstrap ΔNB (pred1 - pred2)
paired_delta_boot <- function(obs, pred1, pred2, thresh, B = B_BOOT) {
  df_boot <- data.frame(obs = obs, pred1 = pred1, pred2 = pred2)
  boot_stat <- function(data, i) {
    nb_single(data[i, "obs"], data[i, "pred1"], thresh) -
      nb_single(data[i, "obs"], data[i, "pred2"], thresh)
  }
  b  <- boot::boot(data = df_boot, statistic = boot_stat, R = B)
  ci <- tryCatch(boot::boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
  list(delta = b$t0, ci_low = ci[1], ci_high = ci[2], ci_width = ci[2] - ci[1])
}

# Calibration summary — one-row data.frame per model
calibration_summary <- function(obs, pred, label) {
  obs  <- as.integer(obs)
  pred <- as.numeric(pred)
  logit_safe <- function(p) {
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    log(p / (1 - p))
  }
  roc_obj   <- pROC::roc(obs, pred, quiet = TRUE)
  auc_val   <- round(as.numeric(pROC::auc(roc_obj)), 3)
  auc_ci    <- round(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), 3)
  cal_fit   <- glm(obs ~ logit_safe(pred), family = binomial())
  intercept <- round(coef(cal_fit)[1], 3)
  slope     <- round(coef(cal_fit)[2], 3)
  mean_pred <- round(mean(pred), 3)
  mean_obs  <- round(mean(obs),  3)
  cal_itl   <- round(mean_pred - mean_obs, 3)
  data.frame(
    model            = label,
    n                = length(obs),
    events           = sum(obs),
    AUC              = auc_val,
    AUC_lo           = auc_ci[1],
    AUC_hi           = auc_ci[3],
    Cal_intercept    = intercept,
    Cal_slope        = slope,
    slope_flag       = if (slope < 0.7 | slope > 1.3) "CHECK slope far from 1" else "OK",
    Mean_pred        = mean_pred,
    Observed_rate    = mean_obs,
    Cal_in_the_large = cal_itl,
    itl_flag         = if (abs(cal_itl) > 0.05) "CHECK >0.05 NB may be biased" else "OK",
    stringsAsFactors = FALSE
  )
}


# ============================================================
# LOAD DATA
# ============================================================
dfz <- read_csv(data_csv, show_col_types = FALSE)

load(train60_file)   # -> full_meta_train
load(test40_file)    # -> full_meta_test

train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)

n_test    <- nrow(test40)
ev_test   <- sum(test40$group01 == 1)   # remission = event
prev_test <- mean(test40$group01)       # prevalence of remission

cat(sprintf("train60 : n=%d | remission=%d\n", nrow(train60), sum(train60$group01 == 1)))
cat(sprintf("test40  : n=%d | remission=%d | prevalence=%.3f\n\n",
            n_test, ev_test, prev_test))

stopifnot(
  "CRP_Plasma_mg_L" %in% colnames(dfz),
  "zANCA"           %in% colnames(dfz),
  "AAV_group"       %in% colnames(dfz)
)

# ============================================================
# COHORT COMPOSITION AUDIT
# ============================================================
cat(strrep("=", 65), "\n")
cat("COHORT COMPOSITION AUDIT\n")
cat(strrep("=", 65), "\n\n")

if ("cohort" %in% colnames(test40)) {
  cohort_counts <- test40 %>%
    count(cohort, name = "n") %>%
    mutate(
      remission = sapply(cohort, function(cc)
        sum(test40$group01[test40$cohort == cc] == 1)),
      prev_remission = round(remission / n, 3)
    )
  cat("test40 cohort breakdown (event = remission):\n")
  print(cohort_counts)
  cat(sprintf(
    "\ntest40: %d Berlin + %d Prague patients.\n",
    sum(test40$cohort == "Berlin", na.rm = TRUE),
    sum(test40$cohort == "Prague",  na.rm = TRUE)
  ))
  cat("Protein panel selected on Berlin-only data before combined\n")
  cat("dataset was assembled — test40 patients were never used in\n")
  cat("feature selection. The 60/40 DCA is a clean held-out check.\n\n")
  write.csv(cohort_counts,
            file.path(outdir, "test40_cohort_composition.csv"), row.names = FALSE)
} else {
  warning("'cohort' column not found in test40 — cannot verify composition.")
}

# ============================================================
# LOAD FROZEN PROTEIN SCORES
# ============================================================
cat(strrep("=", 65), "\n")
cat("FROZEN SCORE LOADING\n")
cat(strrep("=", 65), "\n\n")

frozen_scores <- read_csv(frozen_scores_csv, show_col_types = FALSE)

test40 <- test40 %>%
  left_join(frozen_scores %>% select(Patient, pred_train60_7), by = "Patient")

n_missing <- sum(is.na(test40$pred_train60_7))
if (n_missing > 0) {
  stop(n_missing, " test40 patients have no frozen score. ",
       "Re-run TQL_models.R on the same patient list.")
}
cat(sprintf("OK: frozen scores matched for all %d test40 patients.\n\n", n_test))

# protein_score = P(remission) from the frozen 7-panel model
# Used directly in DCA C, D, and all manual bootstrap tables.
test40$protein_score <- test40$pred_train60_7

# train60 protein_score — input to augmented wrapper GLMs only
model_train60_7 <- readRDS(file.path(model_dir, "model_7panel_Train60_glm.rds"))
cat("Main model formula:", deparse(formula(model_train60_7)), "\n\n")
train60$protein_score <- predict(model_train60_7, newdata = train60, type = "response")

# ============================================================
# FIT WRAPPER COMPARISON MODELS ON TRAIN60
# Single-predictor logistic models used for DCA A/B (all-model
# supplement). All predict group01 (1 = remission).
# NOT used in DCA C (main figure) — that uses protein_score.
# ============================================================
mod_crp           <- glm(group01 ~ CRP_Plasma_mg_L,                        data = train60, family = binomial())
mod_anca          <- glm(group01 ~ zANCA,                                   data = train60, family = binomial())
mod_prot          <- glm(group01 ~ protein_score,                           data = train60, family = binomial())
mod_anca_crp      <- glm(group01 ~ zANCA + CRP_Plasma_mg_L,             data = train60, family = binomial())
mod_prot_crp      <- glm(group01 ~ protein_score + CRP_Plasma_mg_L,         data = train60, family = binomial())
mod_prot_anca     <- glm(group01 ~ protein_score + zANCA,                   data = train60, family = binomial())
mod_prot_crp_anca <- glm(group01 ~ protein_score + CRP_Plasma_mg_L + zANCA, data = train60, family = binomial())

for (nm in c("mod_crp", "mod_anca", "mod_prot", "mod_anca_crp",
             "mod_prot_crp", "mod_prot_anca", "mod_prot_crp_anca")) {
  if (!isTRUE(get(nm)$converged))
    warning(nm, " did not converge — check for complete separation.")
}

# Score test40 — all wrapper scores predict P(remission)
test40$prob_prot          <- predict(mod_prot,          newdata = test40, type = "response")
test40$prob_crp           <- predict(mod_crp,           newdata = test40, type = "response")
test40$prob_anca          <- predict(mod_anca,          newdata = test40, type = "response")
test40$prob_anca_crp      <- predict(mod_anca_crp,          newdata = test40, type = "response")
test40$prob_prot_crp      <- predict(mod_prot_crp,      newdata = test40, type = "response")
test40$prob_prot_anca     <- predict(mod_prot_anca,     newdata = test40, type = "response")
test40$prob_prot_crp_anca <- predict(mod_prot_crp_anca, newdata = test40, type = "response")

# ============================================================
# SANITY CHECK 1 — FROZEN SCORE CALIBRATION
# Cross-check against calibration_all_models.R output.
# That script is the authoritative source — avoids computing
# the same statistics twice with potentially different code.
# ============================================================

if (file.exists(cal_output_csv)) {
  cal_summary <- read_csv(cal_output_csv, show_col_types = FALSE)
  row_7p <- cal_summary[cal_summary$model == "Train60 7-protein -> Test40", ]
  if (nrow(row_7p) == 0) {
    warning("'Train60 7-protein -> Test40' not found in calibration_all_models.csv.")
  } else {
    cat("Stored calibration (authoritative):\n\n")
    print(row_7p[, c("model", "AUC", "AUC_lower95", "AUC_upper95",
                     "Cal_intercept", "Cal_slope", "HL_p",
                     "Mean_pred", "Observed_rate", "Cal_in_the_large")])
    cat("\n")
    itl   <- row_7p$Cal_in_the_large
    slope <- row_7p$Cal_slope
    cat(sprintf("  Cal-in-the-large: %+.3f  %s\n", itl,
                if (abs(itl) > 0.05) "*** CHECK >0.05 — NB may be biased ***" else "OK"))
    cat(sprintf("  Cal slope       :  %.3f  %s\n", slope,
                if (slope < 0.7 | slope > 1.3) "*** CHECK far from 1 ***" else "OK"))
    cat(sprintf("  HL p-value      :  %.3f  %s\n\n", row_7p$HL_p,
                if (!is.na(row_7p$HL_p) && row_7p$HL_p < 0.05)
                  "*** CHECK p<0.05 — significant lack of fit ***" else "OK"))
    # Cross-check: live mean should match stored
    live_mean <- mean(test40$protein_score, na.rm = TRUE)
    delta_x   <- abs(live_mean - row_7p$Mean_pred)
    cat(sprintf("  Live mean(protein_score)=%.4f | stored=%.4f | delta=%.4f  %s\n\n",
                live_mean, row_7p$Mean_pred, delta_x,
                if (delta_x > 0.005)
                  "*** MISMATCH — check frozen scores are consistent ***"
                else "OK (scores consistent)"))
  }
} else {
  warning("calibration_all_models.csv not found at: ", cal_output_csv,
          "\nRun calibration_all_models.R first.")
  cat("Falling back to live calibration computation:\n")
  live_cal <- calibration_summary(test40$group01, test40$protein_score,
                                  "Train60 7-protein -> Test40 (live, remission=1)")
  print(live_cal)
  cat("\n")
}

# ============================================================
# SANITY CHECK 2 — WRAPPER MODEL CALIBRATION (CRP, ANCA)
# These only affect DCA A/B (supplement all-model figures).
# DCA C (main figure) uses protein_score directly — unaffected.
# ============================================================


wrapper_cal <- bind_rows(
  calibration_summary(test40$group01, test40$prob_crp,  "CRP wrapper  (remission=1)"),
  calibration_summary(test40$group01, test40$prob_anca, "ANCA wrapper (remission=1)")
)
print(wrapper_cal[, c("model", "n", "events", "AUC", "AUC_lo", "AUC_hi",
                      "Cal_intercept", "slope_flag",
                      "Cal_in_the_large", "itl_flag")]) 

write.csv(wrapper_cal, file.path(outdir, "wrapper_calibration_CRP_ANCA.csv"), row.names = FALSE)

for (col in c("prob_crp", "prob_anca")) {
  mn   <- mean(test40[[col]], na.rm = TRUE)
  diff <- mn - prev_test
  cat(sprintf("  %-12s  mean_pred=%.3f  obs(remission)=%.3f  diff=%+.3f  %s\n",
              col, mn, prev_test, diff,
              if (abs(diff) > 0.05) "*** CHECK — NB in DCA A/B may be biased ***" else "OK"))
}
cat("\n  Miscalibration here affects DCA A/B (supplement) only.\n")
cat("  DCA C (main figure) is unaffected.\n\n")

# ============================================================
# SANITY CHECK 3 — AUC ON TEST40

cat("SANITY CHECK 3 — AUC ON TEST40 (event = remission)\n")


auc_cols <- c("protein_score", "prob_crp", "prob_anca", "prob_anca_crp",
              "prob_prot_crp", "prob_prot_anca", "prob_prot_crp_anca")
auc_df <- data.frame(
  model = c("7-protein panel (frozen)", "CRP alone (wrapper)",
            "ANCA alone (wrapper)", "ANCA + CRP", "7-protein + CRP",
            "7-protein + ANCA", "7-protein + CRP + ANCA"),
  AUC = sapply(auc_cols, function(col)
    round(as.numeric(pROC::auc(
      pROC::roc(test40$group01, test40[[col]], quiet = TRUE))), 3)),
  row.names = NULL
)
print(auc_df)

# ============================================================
# SANITY CHECK 4 — SCORE DISTRIBUTIONS
# ============================================================
cat("SANITY CHECK 4 — SCORE RANGES ON TEST40\n")

for (col in c("protein_score", "prob_crp", "prob_anca", "prob_anca_crp",
              "prob_prot_crp", "prob_prot_anca", "prob_prot_crp_anca")) {
  q <- quantile(test40[[col]], c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  cat(sprintf("  %-28s  Q5=%.3f Q25=%.3f Med=%.3f Q75=%.3f Q95=%.3f\n",
              col, q[1], q[2], q[3], q[4], q[5]))
}

# ============================================================
# SANITY CHECK 5 — TREAT-ALL NET BENEFIT
# "Treat all" = taper ALL patients regardless of score.
# Prevalence here = P(remission) = proportion safe to taper.
# ============================================================
cat("SANITY CHECK 5 — TREAT-ALL NB (taper everyone)\n")


treat_all_nb <- function(prev, pt) prev - (1 - prev) * (pt / (1 - pt))
for (pt in c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)) {
  cat(sprintf("  threshold=%.2f -> NB_treat_all=%.4f\n",
              pt, treat_all_nb(prev_test, pt)))
}
cat(sprintf(
  "\n  P(remission)=%.0f%% — treat-all NB is high at low thresholds\n",
  prev_test * 100))
cat("  because most patients are actually in remission.\n")
cat("  Panel advantage over treat-all grows at higher thresholds.\n\n")

# ============================================================
# SCORE DISTRIBUTION PLOT (diagnostic only)
# ============================================================
dist_df <- test40 %>%
  select(group01, protein_score, prob_crp, prob_anca) %>%
  pivot_longer(cols = c(protein_score, prob_crp, prob_anca),
               names_to = "model", values_to = "prob") %>%
  mutate(
    model   = recode(model,
                     protein_score = "7-protein panel (frozen, P(remission))",
                     prob_crp      = "CRP alone (wrapper)",
                     prob_anca     = "ANCA alone (wrapper)"),
    outcome = factor(group01, labels = c("Active (0)", "Remission (1)"))
  )

p_dist <- ggplot(dist_df, aes(x = prob, fill = outcome)) +
  geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Active (0)" = "#c70000", "Remission (1)" = "#2ecc71")) +
  labs(title    = "Score distributions by outcome — DIAGNOSTIC",
       subtitle = "Scores predict P(remission). Not for paper.",
       x = "Predicted probability of remission", y = "Count", fill = "Outcome") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

svglite(file.path(outdir, "DIAGNOSTIC_score_distributions.svg"), width = 5.5, height = 7)
print(p_dist); dev.off()

# ============================================================
# DCA A — FULL RANGE (5–50%), all models
# All scores predict P(remission); group01==1 is the event.
# prob_prot = wrapper GLM (consistent with other wrappers here).
# ============================================================
cat("Running DCA A: full range, all models...\n")

dca_full <- dca_with_ci(
  formula = group01 ~ prob_prot + prob_crp + prob_anca + prob_anca_crp+
    prob_prot_crp + prob_prot_anca + prob_prot_crp_anca,
  data       = test40,
  thresholds = THRESH_FULL,
  label      = list(
    prob_prot          = "7-protein panel",
    prob_crp           = "CRP alone",
    prob_anca          = "ANCA alone",
    prob_anca_crp      = "ANCA + CRP",
    prob_prot_crp      = "7-protein + CRP",
    prob_prot_anca     = "7-protein + ANCA",
    prob_prot_crp_anca = "7-protein + CRP + ANCA"
  )
)

p_full <- dca_full %>%
  plot(smooth = TRUE, rug = FALSE) +
  scale_color_manual(values = model_cols) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  annotate("rect", xmin = 0.20, xmax = 0.50, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.35) +
  geom_vline(xintercept = clin_thresh, colour = "grey70", linetype = "dashed", size = 0.3) +
  annotate("text", x = 0.35, y = 0.02, size = 2.8, colour = "grey40",
           label = "Panel consistently\nsuperior to CRP (>=20%)") +
  labs(
    title    = sprintf("Decision curve analysis — test set (n=%d) | Full range 5–50%%", n_test),
    subtitle = sprintf("Outcome: remission | Prevalence = %d%% | Shaded = clinically relevant range",
                       round(prev_test * 100)),
    x = "Threshold probability", y = "Net benefit"
  ) +
  dca_theme_nature(11)

svglite(file.path(outdir, "DCA_A_full_range_all_models_improved.svg"), width = 8, height = 5.5)
print(p_full); dev.off()


df_full <- dca_to_tibble(dca_full)
write.csv(df_full, file.path(outdir, "DCA_full_range_table.csv"), row.names = FALSE)

# ============================================================
# DCA B — CLINICAL RANGE (15–50%), all models
# ============================================================
dca_clin <- dca_with_ci(
  formula = group01 ~ prob_prot + prob_crp + prob_anca + prob_anca_crp+
    prob_prot_crp + prob_prot_anca + prob_prot_crp_anca,
  data       = test40,
  thresholds = THRESH_CLINICAL,
  label      = list(
    prob_prot          = "7-protein panel",
    prob_crp           = "CRP alone",
    prob_anca          = "ANCA alone",
    prob_anca_crp      = "ANCA + CRP",
    prob_prot_crp      = "7-protein + CRP",
    prob_prot_anca     = "7-protein + ANCA",
    prob_prot_crp_anca = "7-protein + CRP + ANCA"
  )
)

p_clin <- dca_clin %>%
  plot(smooth = TRUE) +
  scale_color_manual(values = model_cols) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  geom_vline(xintercept = clin_thresh, colour = "grey80", linetype = "dashed", size = 0.25) +
  labs(
    title    = sprintf("Decision curve analysis — clinical range 15–50%% (n=%d)", n_test),
    subtitle = sprintf("Outcome: remission | Prevalence = %d%%", round(prev_test * 100)),
    x = "Threshold probability", y = "Net benefit"
  ) +
  dca_theme_nature(11)

svglite(file.path(outdir, "DCA_B_clinical_range_all_models_improved.svg"), width = 8, height = 5.5)
print(p_clin); dev.off()


df_clin <- dca_to_tibble(dca_clin)
write.csv(df_clin, file.path(outdir, "DCA_clinical_range_table.csv"), row.names = FALSE)

# ============================================================
# DCA C — FOCUSED: 7-protein panel vs CRP alone
# MAIN PAPER FIGURE.
# Uses frozen protein_score (P(remission)) directly — not the
# logistic wrapper — for the panel. CRP uses wrapper prob_crp.
# Both predict P(remission); group01==1 is the event.
# ============================================================
dca_vs_crp <- dca_with_ci(
  formula    = group01 ~ protein_score + prob_crp,
  data       = test40,
  thresholds = THRESH_CLINICAL,
  label      = list(protein_score = "7-protein panel", prob_crp = "CRP alone")
)

p_vs_crp <- dca_vs_crp %>%
  plot(smooth = TRUE) +
  scale_color_manual(values = c("7-protein panel" = "#C70500",
                                "CRP alone"       = "#89CC48",
                                "All"             = "#D5CD78",
                                "None"            = "black")) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  annotate("rect", xmin = 0.20, xmax = 0.50, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.35) +
  annotate("text", x = 0.345, y = 0.02, size = 2.8, colour = "grey40",
           label = "Panel consistently\nsuperior to CRP (>=20%)") +
  geom_vline(xintercept = clin_thresh, colour = "grey80", linetype = "dashed", size = 0.25) +
  labs(
    title    = sprintf("Decision curve analysis — 7-protein panel vs CRP | test set (n=%d)", n_test),
    subtitle = sprintf("Remission: n=%d (%.0f%%) | Threshold = minimum certainty to taper",
                       ev_test, prev_test * 100),
    x = "Threshold probability of remission", y = "Net benefit"
  ) +
  dca_theme_nature(11)

svglite(file.path(outdir, "DCA_C_panel_vs_CRP.svg"), width = 6.5, height = 5)
print(p_vs_crp); dev.off()

p_vs_crp <- dca_vs_crp %>%
  plot(smooth = TRUE) +
  scale_color_manual(values = c("7-protein panel" = "#C70500",
                                "CRP alone"       = "#89CC48",
                                "All"             = "#D5CD78",
                                "None"            = "black")) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  annotate("rect", xmin = 0.15, xmax = 0.25, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.40) +
  annotate("text", x = 0.20, y = -0.03, size = 2.5, colour = "grey40",
           label = "Panel > CRP\n(CI excl. zero)") +
  geom_vline(xintercept = clin_thresh, colour = "grey80", linetype = "dashed", size = 0.25) +
  labs(
    title    = sprintf("Decision curve analysis — 7-protein panel vs CRP | test set (n=%d)", n_test),
    subtitle = sprintf("Remission: n=%d (%.0f%%) | Threshold = minimum certainty to taper",
                       ev_test, prev_test * 100),
    x = "Threshold probability of remission", y = "Net benefit"
  ) +
  dca_theme_nature(11)


svglite(file.path(outdir, "DCA_C_panel_vs_CRP_final.svg"), width = 6.5, height = 3.5)
print(p_vs_crp); dev.off()


# Extract NB table and raw difference
df_vs <- dca_to_tibble(dca_vs_crp)
write.csv(df_vs, file.path(outdir, "DCA_vs_crp_tibble.csv"), row.names = FALSE)

nb_wide <- df_vs %>%
  filter(variable %in% c("protein_score", "prob_crp")) %>%
  select(variable, threshold, net_benefit) %>%
  pivot_wider(names_from = variable, values_from = net_benefit) %>%
  arrange(threshold)

if (!all(c("protein_score", "prob_crp") %in% colnames(nb_wide)))
  stop("Pivot failed. Found columns: ", paste(colnames(nb_wide), collapse = ", "))

adv_tbl <- nb_wide %>%
  rename(nb_panel = protein_score, nb_crp = prob_crp) %>%
  mutate(diff = nb_panel - nb_crp)

write.csv(adv_tbl, file.path(outdir, "DCA_panel_vs_CRP_nb_table_simple.csv"), row.names = FALSE)

# Annotated version of Figure C
annot_df <- adv_tbl %>%
  filter(threshold %in% annot_thresholds) %>%
  mutate(ypos = pmax(nb_panel, nb_crp, na.rm = TRUE) + 0.01)

if (nrow(annot_df) > 0) {
  p_vs_crp_annot <- p_vs_crp +
    geom_point(data = annot_df, aes(x = threshold, y = ypos),
               inherit.aes = FALSE, shape = 21, fill = "white", size = 2) +
    geom_text(data = annot_df, aes(x = threshold, y = ypos,
                                   label = sprintf("\u0394NB=%.3f", diff)),
              inherit.aes = FALSE, size = 2.8, vjust = 0)
  svglite(file.path(outdir, "DCA_C_panel_vs_CRP_annotated.svg"), width = 6.5, height = 5)
  print(p_vs_crp_annot); dev.off()
}

# ============================================================
# NET BENEFIT TABLE (from DCA A — all models)
# ============================================================
target_thresholds <- round(seq(0.10, 0.50, by = 0.05), 10)

nb_full <- dca_full %>%
  as_tibble() %>%
  filter(round(threshold, 10) %in% target_thresholds) %>%
  select(variable, threshold, net_benefit) %>%
  pivot_wider(id_cols = threshold, names_from = variable, values_from = net_benefit) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

if (!0.30 %in% round(nb_full$threshold, 2))
  warning("threshold=0.30 missing — check THRESH_FULL.")

write.csv(nb_full, file.path(outdir, "DCA_net_benefit_full.csv"), row.names = FALSE)

nb_advantage <- nb_full %>%
  mutate(
    prot_vs_crp_advantage = prob_prot - prob_crp,
    prot_vs_all_advantage = prob_prot - all
  ) %>%
  select(threshold, prob_prot, prob_crp, prob_anca, all,
         prot_vs_crp_advantage, prot_vs_all_advantage)

write.csv(nb_advantage, file.path(outdir, "DCA_protein_advantage_table.csv"), row.names = FALSE)

# ============================================================
# DCA D — NET INTERVENTIONS AVOIDED (supplement)
# Uses frozen protein_score, event = remission.
# "Treat all" baseline = taper every patient.
# ============================================================
cat("Running DCA D: net interventions avoided...\n")

dca_nia <- dca_with_ci(
  formula    = group01 ~ protein_score,
  data       = test40,
  thresholds = THRESH_CLINICAL,
  label      = list(protein_score = "7-protein panel")
)

p_nia <- dca_nia %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE) +
  scale_color_manual(values = c("7-protein panel" = "#C70500", "All" = "#D5CD78")) +
  coord_cartesian(ylim = c(-5, NA)) +
  labs(
    title    = sprintf("Net interventions avoided — 7-protein panel | test set (n=%d)", n_test),
    subtitle = sprintf("Per 100 patients vs tapering all | Remission prevalence = %d%%",
                       round(prev_test * 100)),
    x = "Threshold probability of remission",
    y = "Net reduction in unnecessary tapers per 100 patients"
  ) +
  dca_theme_nature(11)

svglite(file.path(outdir, "DCA_D_interventions_avoided.svg"), width = 6.5, height = 5)
print(p_nia); dev.off()


# ============================================================
# BOOTSTRAP ΔNB — PANEL MINUS CRP AT EACH THRESHOLD
# Manual bootstrap on frozen scores (not through dcurves).
# ============================================================
set.seed(2023)

thresh_vec  <- sort(unique(adv_tbl$threshold))
nb_diff_mat <- matrix(NA_real_, nrow = B_DIFF, ncol = length(thresh_vec),
                      dimnames = list(NULL, as.character(thresh_vec)))

for (b_i in seq_len(B_DIFF)) {
  idx <- sample(seq_len(nrow(test40)), replace = TRUE)
  dfb <- test40[idx, , drop = FALSE]
  nb_diff_mat[b_i, ] <- sapply(thresh_vec, function(t) {
    nb_single(dfb$group01, dfb$protein_score, t) -
      nb_single(dfb$group01, dfb$prob_crp, t)
  })
}

ci_tbl <- data.frame(
  threshold = thresh_vec,
  diff      = adv_tbl %>% arrange(threshold) %>% pull(diff),
  lower     = apply(nb_diff_mat, 2, quantile, 0.025, na.rm = TRUE),
  upper     = apply(nb_diff_mat, 2, quantile, 0.975, na.rm = TRUE)
)
write.csv(ci_tbl, file.path(outdir, "DCA_panel_minus_CRP_diff_CI.csv"), row.names = FALSE)

p_adv <- ggplot(ci_tbl, aes(x = threshold, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C70500", alpha = 0.18) +
  geom_line(colour = "#C70500", size = 1) +
  geom_point(data = ci_tbl %>% filter(threshold %in% annot_thresholds),
             aes(x = threshold, y = diff), size = 1.8) +
  geom_text(data = ci_tbl %>% filter(threshold %in% annot_thresholds),
            aes(x = threshold, y = diff, label = sprintf("%.3f", diff)),
            vjust = -0.8, size = 2.8) +
  labs(
    title    = "Net-benefit advantage: 7-protein panel \u2212 CRP",
    subtitle = sprintf("Positive = panel better for remission identification | %d bootstrap reps",
                       B_DIFF),
    x = "Threshold probability of remission",
    y = "Net-benefit difference (panel \u2212 CRP)"
  ) +
  coord_cartesian(ylim = c(min(ci_tbl$lower, na.rm = TRUE) - 0.01, NA)) +
  theme_bw(base_size = 11)

svglite(file.path(outdir, "DCA_panel_minus_CRP_advantage_with_bootCI.svg"), width = 6.5, height = 4)
print(p_adv); dev.off()

# ============================================================
# COHORT-SPECIFIC DCA AND BOOTSTRAP NB TABLES
# ============================================================
cohorts <- sort(unique(test40$cohort[!is.na(test40$cohort)]))

for (cc in cohorts) {
  subset_df <- test40 %>% filter(cohort == cc)
  n_c  <- nrow(subset_df)
  ev_c <- sum(subset_df$group01 == 1)
  message(sprintf("Running cohort DCA for '%s' (n=%d, remission=%d)", cc, n_c, ev_c))
  
  dca_cohort <- dca_with_ci(
    formula    = group01 ~ prob_prot + prob_crp + prob_anca +
      prob_prot_crp + prob_prot_anca + prob_prot_crp_anca,
    data       = subset_df,
    thresholds = THRESH_FULL,
    label      = list(
      prob_prot          = "7-protein panel",
      prob_crp           = "CRP alone",
      prob_anca          = "ANCA alone",
      prob_prot_crp      = "7-protein + CRP",
      prob_prot_anca     = "7-protein + ANCA",
      prob_prot_crp_anca = "7-protein + CRP + ANCA"
    )
  )
  
  p_cohort <- dca_cohort %>%
    plot(smooth = TRUE) +
    scale_color_manual(values = model_cols) +
    labs(title = sprintf("DCA — cohort: %s | n=%d | remission=%d", cc, n_c, ev_c)) +
    dca_theme_nature(11)
  
  ggsave(file.path(outdir, sprintf("DCA_cohort_%s.svg", cc)),
         plot = p_cohort, width = 8, height = 5)
  
  # Bootstrap NB table with CIs (using frozen score for panel)
  nb_rows <- lapply(THRESH_REPORT, function(th) {
    panel_res <- nb_boot_ci(subset_df$group01, subset_df$protein_score, th, B = B_BOOT)
    crp_res   <- nb_boot_ci(subset_df$group01, subset_df$prob_crp,       th, B = B_BOOT)
    data.frame(cohort = cc, threshold = th,
               nb_panel      = panel_res$nb,    ci_low_panel  = panel_res$ci_low,
               ci_high_panel = panel_res$ci_high,
               nb_crp        = crp_res$nb,      ci_low_crp    = crp_res$ci_low,
               ci_high_crp   = crp_res$ci_high,
               delta_nb      = panel_res$nb - crp_res$nb,
               stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, nb_rows),
            file.path(outdir, sprintf("NB_table_cohort_%s.csv", cc)), row.names = FALSE)
  message("Saved: ", sprintf("NB_table_cohort_%s.csv", cc))
}

# After the existing cohort CRP loop, add:
for (cc in cohorts) {
  df_c   <- test40 %>% filter(cohort == cc)
  rows_c <- lapply(THRESH_REPORT, function(th) {
    pb <- paired_delta_boot(df_c$group01, df_c$protein_score, df_c$prob_anca,
                            thresh = th, B = B_BOOT)
    data.frame(scope = paste0("cohort_", cc), threshold = th,
               delta_hat = pb$delta, ci_low = pb$ci_low, ci_high = pb$ci_high,
               ci_width = pb$ci_width, stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, rows_c),
            file.path(outdir, sprintf("paired_deltaNB_ANCA_cohort_%s.csv", cc)),
            row.names = FALSE)
}

# ============================================================
# PREVALENCE-WEIGHTED SENSITIVITY ANALYSIS
# ============================================================
prev_berlin <- mean(test40$group01[test40$cohort == "Berlin"], na.rm = TRUE)
prev_prague <- mean(test40$group01[test40$cohort == "Prague"],  na.rm = TRUE)

compute_weights <- function(obs, target_prev) {
  n_pos <- sum(obs == 1); n_neg <- sum(obs == 0)
  if (n_pos == 0 || n_neg == 0) return(rep(1, length(obs)))
  a <- (target_prev * n_neg) / ((1 - target_prev) * n_pos)
  ifelse(obs == 1, a, 1)
}

targets <- data.frame(
  name        = c("pooled_observed", "Berlin_prev", "Prague_prev"),
  target_prev = c(prev_test, prev_berlin, prev_prague),
  stringsAsFactors = FALSE
)

for (ti in seq_len(nrow(targets))) {
  tgt_name <- targets$name[ti]
  tgt_prev <- targets$target_prev[ti]
  wts      <- compute_weights(test40$group01, tgt_prev)
  
  nb_rows <- lapply(THRESH_REPORT, function(th) {
    res_panel <- nb_boot_ci(test40$group01, test40$protein_score, th, weights = wts, B = B_BOOT)
    res_crp   <- nb_boot_ci(test40$group01, test40$prob_crp,      th, weights = wts, B = B_BOOT)
    data.frame(target = tgt_name, target_prev = tgt_prev, threshold = th,
               nb_panel      = res_panel$nb,    ci_low_panel  = res_panel$ci_low,
               ci_high_panel = res_panel$ci_high,
               nb_crp        = res_crp$nb,      ci_low_crp    = res_crp$ci_low,
               ci_high_crp   = res_crp$ci_high,
               delta_nb      = res_panel$nb - res_crp$nb,
               stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, nb_rows),
            file.path(outdir, sprintf("NB_weighted_target_%s.csv", tgt_name)),
            row.names = FALSE)
}

# ============================================================
# PAIRED BOOTSTRAP ΔNB — POOLED AND PER COHORT
# ============================================================

set.seed(1234)

pooled_rows <- lapply(THRESH_REPORT, function(th) {
  pb <- paired_delta_boot(test40$group01, test40$protein_score, test40$prob_crp,
                          thresh = th, B = B_BOOT)
  data.frame(scope = "pooled", threshold = th,
             delta_hat = pb$delta, ci_low = pb$ci_low, ci_high = pb$ci_high,
             ci_width = pb$ci_width, stringsAsFactors = FALSE)
})
pooled_df <- do.call(rbind, pooled_rows)
write.csv(pooled_df, file.path(outdir, "paired_deltaNB_pooled.csv"), row.names = FALSE)

for (cc in cohorts) {
  df_c   <- test40 %>% filter(cohort == cc)
  rows_c <- lapply(THRESH_REPORT, function(th) {
    pb <- paired_delta_boot(df_c$group01, df_c$protein_score, df_c$prob_crp,
                            thresh = th, B = B_BOOT)
    data.frame(scope = paste0("cohort_", cc), threshold = th,
               delta_hat = pb$delta, ci_low = pb$ci_low, ci_high = pb$ci_high,
               ci_width = pb$ci_width, stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, rows_c),
            file.path(outdir, sprintf("paired_deltaNB_cohort_%s.csv", cc)), row.names = FALSE)
}

# ============================================================
# COHORT INTERACTION TEST
# Non-significant = panel effect is consistent across cohorts,
# pooling is defensible.
# ============================================================

int_model <- glm(group01 ~ protein_score * cohort, data = test40, family = binomial())
cat("Interaction term (protein_score:cohortPrague):\n")
print(coef(summary(int_model))["protein_score:cohortPrague", , drop = FALSE])
cat("\nIf p >> 0.05: no evidence of effect modification by cohort.\n")
cat("Pooling is statistically defensible.\n\n")

# ============================================================
# CONSOLE SUMMARY
# ============================================================

print(pooled_df)

for (cc in cohorts) {
  message(sprintf("\nΔNB cohort = %s:", cc))
  print(read.csv(file.path(outdir, sprintf("paired_deltaNB_cohort_%s.csv", cc))))
}


test40 %>% filter(CRPhigher5_YES1 == 1) %>% 
  group_by(group01) %>% 
  summarise(median = median(pred_train60_7), IQR_lo = quantile(pred_train60_7, 0.25), 
            IQR_hi = quantile(pred_train60_7, 0.75))

# ============================================================
# DCA B2 — FOCUSED: 7-protein panel vs CRP alone vs ANCA alone
# Shows all three single-marker comparators in one figure.
# Uses wrapper scores for CRP and ANCA (consistent with DCA A/B).
# Uses frozen protein_score for panel.
# ============================================================

# "7-protein panel"        = "#C70500",
# "CRP alone"              = "#89CC48",
# "ANCA alone"             = "#66CCD1",
# "ANCA + CRP"   = "#A260D8",
# "7-protein + CRP"        = "#888888",
# "7-protein + ANCA"       = "#aaaaaa",
# "7-protein + CRP + ANCA" = "#BBC2D5",
# "All"                    = "#D5CD78",
# "None"                   = "black"


dca_vs_crp_anca <- dca_with_ci(
  formula    = group01 ~ protein_score + prob_crp + prob_anca,
  data       = test40,
  thresholds = THRESH_CLINICAL,
  label      = list(
    protein_score = "7-protein panel",
    prob_crp      = "CRP alone",
    prob_anca     = "ANCA alone"
  )
)

cols_b2 <- c(
  "7-protein panel" = "#C70500",
  "CRP alone"       = "#89CC48",
  "ANCA alone"      = "#66CCD1",
  "All"             = "#D5CD78",
  "None"            = "black"
)

p_vs_crp_anca <- dca_vs_crp_anca %>%
  plot(smooth = TRUE) +
  scale_color_manual(values = cols_b2) +
  scale_linetype_manual(values = c(
    "7-protein panel" = "solid",
    "CRP alone"       = "dashed",
    "ANCA alone"      = "dotdash",
    "All"             = "solid",
    "None"            = "solid"
  )) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  annotate("rect",
           xmin = 0.15, xmax = 0.26,
           ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.40) +
  annotate("text",
           x = 0.205, y = -0.035, size = 2.5, colour = "grey35",
           label = "Panel > CRP\n(CI excl. 0)") +
  geom_vline(xintercept = clin_thresh,
             colour = "grey75", linetype = "dashed", size = 0.25) +
  labs(
    title    = sprintf("Decision curve analysis — panel vs CRP vs ANCA | test set (n=%d)", n_test),
    subtitle = sprintf("Remission: n=%d (%.0f%%) | Shaded = thresholds where panel CI excludes zero vs CRP",
                       ev_test, round(prev_test * 100)),
    x        = "Threshold probability of remission",
    y        = "Net benefit",
    colour   = NULL,
    linetype = NULL
  ) +
  dca_theme_nature(11) +
  theme(legend.position = "right")

svglite(file.path(outdir, "DCA_B2_panel_vs_CRP_vs_ANCA.svg"),
        width = 7, height = 5)
print(p_vs_crp_anca)
dev.off()

# NB table for this comparison
df_b2 <- dca_to_tibble(dca_vs_crp_anca)
write.csv(df_b2, file.path(outdir, "DCA_B2_panel_vs_CRP_vs_ANCA_table.csv"), row.names = FALSE)

# ΔNB panel vs ANCA (paired bootstrap, same method as panel vs CRP)
cat("Computing bootstrap CI for NB difference (panel - ANCA)...\n")
set.seed(7)

nb_diff_anca_mat <- matrix(NA_real_, nrow = B_DIFF, ncol = length(thresh_vec),
                           dimnames = list(NULL, as.character(thresh_vec)))

for (b_i in seq_len(B_DIFF)) {
  idx  <- sample(seq_len(nrow(test40)), replace = TRUE)
  dfb  <- test40[idx, , drop = FALSE]
  nb_diff_anca_mat[b_i, ] <- sapply(thresh_vec, function(t) {
    nb_single(dfb$group01, dfb$protein_score, t) -
      nb_single(dfb$group01, dfb$prob_anca, t)
  })
}

ci_tbl_anca <- data.frame(
  threshold = thresh_vec,
  diff      = sapply(thresh_vec, function(t)
    nb_single(test40$group01, test40$protein_score, t) -
      nb_single(test40$group01, test40$prob_anca, t)),
  lower     = apply(nb_diff_anca_mat, 2, quantile, 0.025, na.rm = TRUE),
  upper     = apply(nb_diff_anca_mat, 2, quantile, 0.975, na.rm = TRUE)
)
write.csv(ci_tbl_anca,
          file.path(outdir, "DCA_panel_minus_ANCA_diff_CI.csv"),
          row.names = FALSE)


# Quick summary at key thresholds
cat("\nΔNB panel vs ANCA at key thresholds:\n")
print(ci_tbl_anca[round(ci_tbl_anca$threshold, 2) %in% c(0.15, 0.20, 0.25, 0.30, 0.40), ])

p_adv_anca <- ggplot(ci_tbl_anca, aes(x = threshold, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C70500", alpha = 0.25) +
  geom_line(colour = "#C70500", size = 1) +
  geom_point(data = ci_tbl_anca %>% filter(threshold %in% annot_thresholds),
             aes(x = threshold, y = diff), size = 1.8) +
  geom_text(data = ci_tbl_anca %>% filter(threshold %in% annot_thresholds),
            aes(x = threshold, y = diff, label = sprintf("%.3f", diff)),
            vjust = -0.8, size = 2.8) +
  labs(
    title    = "Net-benefit advantage: 7-protein panel \u2212 ANCA",
    subtitle = sprintf("Positive = panel better | %d bootstrap reps", B_DIFF),
    x        = "Threshold probability of remission",
    y        = "Net-benefit difference (panel \u2212 ANCA)"
  ) +
  coord_cartesian(
    ylim = c(min(ci_tbl$lower, na.rm = TRUE) - 0.01, 0.5)
  ) +
  # scale_y_continuous(
  #   breaks = c(-0.5, 0.0, 0.5, 1.0, 1.5)
  # ) +
  theme_bw(base_size = 11)

svglite(file.path(outdir, "DCA_panel_minus_ANCA_advantage_with_bootCI.svg"),
        width = 6.5, height = 4)
print(p_adv_anca)


