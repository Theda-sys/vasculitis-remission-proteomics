# ============================================================
# RIDGE BERLIN → PRAGUE — Complete supplement script
# ============================================================
# HONEST FRAMING:
#   Ridge = supplement / sensitivity analysis only.
#   Model retained after observing Prague calibration metrics.
#   Prague is not a fully independent validation for Ridge.
#   GLM remains primary model for all discrimination claims.
# ============================================================

library(pROC)
library(dplyr)
library(readr)
library(glmnet)
library(ResourceSelection)
library(boot)
library(svglite)

set.seed(7)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
data_csv     <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
outdir       <- "~/Documents/MDC/ml_proteomics/26_RevisionNature/PragueTQL/Ridge_full_eval/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prots_7 <- c("AHSG_001", "CLEC3B_020", "COMP_023",
             "F9_027", "LRG1_042", "MCAM_047", "MRC1_073..")

# ----------------------------------------------------------
# LOAD MODELS + DATA
# ----------------------------------------------------------
model_glm   <- readRDS(file.path(artifact_dir, "model_7panel_Berlin_glm.rds"))
model_ridge <- readRDS(file.path(artifact_dir, "model_7panel_Berlin_Ridge_glmnet.rds"))

dfzt        <- read_csv(data_csv)
data_trainB <- dfzt %>% filter(cohort == "Berlin")
data_testP  <- dfzt %>% filter(cohort == "Prague")

cat(sprintf("Berlin n=%d | Prague n=%d | Prague remission=%d | Prague active=%d\n\n",
            nrow(data_trainB), nrow(data_testP),
            sum(data_testP$group01 == 1, na.rm = TRUE),
            sum(data_testP$group01 == 0, na.rm = TRUE)))

# ----------------------------------------------------------
# FROZEN SCORES
# ----------------------------------------------------------
frozen_path <- file.path(artifact_dir, "scores_Prague_GLM_and_Ridge_frozen.csv")
if (file.exists(frozen_path)) {
  frozen <- read_csv(frozen_path)
  if (nrow(frozen) == nrow(data_testP) && all(frozen$Patient == data_testP$Patient)) {
    data_testP$pred_glm   <- frozen$pred_berlin_glm
    data_testP$pred_ridge <- frozen$pred_berlin_ridge
    cat("\u2713 Frozen scores loaded and patient order verified.\n\n")
  } else {
    warning("Patient order mismatch — using live predictions.")
    data_testP$pred_glm   <- as.numeric(predict(model_glm, newdata = data_testP, type = "response"))
    data_testP$pred_ridge <- as.numeric(predict(model_ridge,
                                                newx = as.matrix(data_testP[, prots_7]),
                                                type = "response"))
  }
} else {
  data_testP$pred_glm   <- as.numeric(predict(model_glm, newdata = data_testP, type = "response"))
  data_testP$pred_ridge <- as.numeric(predict(model_ridge,
                                              newx = as.matrix(data_testP[, prots_7]),
                                              type = "response"))
}

# Ridge score in Berlin — needed as input to incremental model fitting
data_trainB$protein_score <- as.numeric(predict(model_ridge,
                                                newx = as.matrix(data_trainB[, prots_7]),
                                                type = "response"))
data_testP$protein_score  <- data_testP$pred_ridge

# ----------------------------------------------------------
# SHARED HELPERS
# ----------------------------------------------------------

# Formats AUC + DeLong CI for plot legend labels
auc_fmt <- function(roc_obj) {
  a  <- as.numeric(pROC::auc(roc_obj))
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  sprintf("AUC %.3f (%.3f\u2013%.3f)", a, ci[1], ci[3])
}

logit_safe <- function(p) {
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  log(p / (1 - p))
}

binom_ci <- function(x, n) {
  if (n == 0 || is.na(x) || is.na(n)) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int[1:2]
}

# ----------------------------------------------------------
# UNIFIED METRICS CORE
# Computes AUC, DeLong CI, and all threshold metrics once.
# Used by every subgroup / stratum / incremental-model call.
# ----------------------------------------------------------
compute_metrics_core <- function(y, p, threshold = 0.5) {
  y <- as.integer(y)
  p <- as.numeric(p)
  
  roc_obj <- pROC::roc(y, p, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  
  pred_bin <- as.integer(p >= threshold)
  TP <- sum(pred_bin == 1 & y == 1)
  TN <- sum(pred_bin == 0 & y == 0)
  FP <- sum(pred_bin == 1 & y == 0)
  FN <- sum(pred_bin == 0 & y == 1)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  
  list(
    roc_obj = roc_obj,
    auc_val = auc_val, auc_ci = auc_ci,
    TP = TP, TN = TN, FP = FP, FN = FN,
    sens = sens, spec = spec, ppv = ppv, npv = npv
  )
}

# Turns a compute_metrics_core result into one tidy data.frame row.
# id_cols: named list of identifier columns prepended (e.g. list(subgroup="Naive"))
metrics_row <- function(m, id_cols = list()) {
  as.data.frame(c(
    id_cols,
    list(
      n          = m$TP + m$TN + m$FP + m$FN,
      events     = m$TP + m$FN,
      non_events = m$TN + m$FP,
      AUC        = round(m$auc_val,   3),
      AUC_lo     = round(m$auc_ci[1], 3),
      AUC_hi     = round(m$auc_ci[3], 3),
      Sens       = round(m$sens, 3),
      Sens_lo    = round(binom_ci(m$TP, m$TP + m$FN)[1], 3),
      Sens_hi    = round(binom_ci(m$TP, m$TP + m$FN)[2], 3),
      Spec       = round(m$spec, 3),
      Spec_lo    = round(binom_ci(m$TN, m$TN + m$FP)[1], 3),
      Spec_hi    = round(binom_ci(m$TN, m$TN + m$FP)[2], 3),
      PPV        = round(m$ppv,  3),
      PPV_lo     = round(binom_ci(m$TP, m$TP + m$FP)[1], 3),
      PPV_hi     = round(binom_ci(m$TP, m$TP + m$FP)[2], 3),
      NPV        = round(m$npv,  3),
      NPV_lo     = round(binom_ci(m$TN, m$TN + m$FN)[1], 3),
      NPV_hi     = round(binom_ci(m$TN, m$TN + m$FN)[2], 3),
      TP = m$TP, TN = m$TN, FP = m$FP, FN = m$FN
    )
  ), stringsAsFactors = FALSE)
}

# ============================================================
# 1) OVERALL DISCRIMINATION + CALIBRATION + THRESHOLD METRICS
# ============================================================
cat("=== 1) OVERALL DISCRIMINATION + CALIBRATION ===\n")

truth      <- as.integer(data_testP$group01)
pred_ridge <- as.numeric(data_testP$pred_ridge)
pred_glm   <- as.numeric(data_testP$pred_glm)
n          <- length(truth)
threshold  <- 0.5

roc_ridge_overall <- pROC::roc(truth, pred_ridge, quiet = TRUE)
roc_glm_overall   <- pROC::roc(truth, pred_glm,   quiet = TRUE)
auc_ridge <- as.numeric(pROC::auc(roc_ridge_overall))
auc_glm   <- as.numeric(pROC::auc(roc_glm_overall))
ci_ridge  <- as.numeric(pROC::ci.auc(roc_ridge_overall, method = "delong"))
ci_glm    <- as.numeric(pROC::ci.auc(roc_glm_overall,   method = "delong"))
delong_glm_ridge <- pROC::roc.test(roc_glm_overall, roc_ridge_overall, method = "delong")

cal_ridge <- coef(glm(truth ~ logit_safe(pred_ridge), family = binomial()))
cal_glm   <- coef(glm(truth ~ logit_safe(pred_glm),   family = binomial()))
hl_ridge  <- tryCatch(ResourceSelection::hoslem.test(truth, pred_ridge, g = 10)$p.value, error = function(e) NA)
hl_glm    <- tryCatch(ResourceSelection::hoslem.test(truth, pred_glm,   g = 10)$p.value, error = function(e) NA)
brier_ridge <- mean((pred_ridge - truth)^2)
brier_glm   <- mean((pred_glm   - truth)^2)

brier_ci_ridge <- { b <- boot::boot(data.frame(p = pred_ridge, y = truth),
                                    function(d, i) mean((d$p[i] - d$y[i])^2), R = 2000)
round(quantile(b$t, c(.025, .975)), 3) }
brier_ci_glm   <- { b <- boot::boot(data.frame(p = pred_glm, y = truth),
                                    function(d, i) mean((d$p[i] - d$y[i])^2), R = 2000)
round(quantile(b$t, c(.025, .975)), 3) }

cat(sprintf("  GLM   AUC=%.3f (%.3f-%.3f) | slope=%.3f | Brier=%.3f | HL p=%.4f\n",
            auc_glm,   ci_glm[1],   ci_glm[3],   cal_glm[2],   brier_glm,   hl_glm))
cat(sprintf("  Ridge AUC=%.3f (%.3f-%.3f) | slope=%.3f | Brier=%.3f | HL p=%.4f\n",
            auc_ridge, ci_ridge[1], ci_ridge[3], cal_ridge[2], brier_ridge, hl_ridge))
cat(sprintf("  DeLong p (GLM vs Ridge): %.4f\n\n", delong_glm_ridge$p.value))

# Threshold metrics — Ridge only (GLM miscalibrated)
m_ov <- compute_metrics_core(truth, pred_ridge, threshold)

cat("  Threshold metrics at 0.5 (Ridge only — GLM miscalibrated):\n")
cat(sprintf("  n=%d | TP=%d TN=%d FP=%d FN=%d\n", n, m_ov$TP, m_ov$TN, m_ov$FP, m_ov$FN))
cat(sprintf("  Sensitivity : %.3f (%.3f-%.3f)\n", m_ov$sens, binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[2]))
cat(sprintf("  Specificity : %.3f (%.3f-%.3f)\n", m_ov$spec, binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[2]))
cat(sprintf("  PPV         : %.3f (%.3f-%.3f)\n", m_ov$ppv,  binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[2]))
cat(sprintf("  NPV         : %.3f (%.3f-%.3f)\n\n", m_ov$npv, binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[2]))

# Save overall table
overall_tab <- data.frame(
  model         = c("Berlin GLM", "Berlin Ridge"),
  n             = c(n, n),
  events        = c(sum(truth), sum(truth)),
  AUC           = round(c(auc_glm, auc_ridge), 3),
  AUC_lo        = round(c(ci_glm[1],    ci_ridge[1]),    3),
  AUC_hi        = round(c(ci_glm[3],    ci_ridge[3]),    3),
  DeLong_p      = c(NA, round(delong_glm_ridge$p.value, 4)),
  Brier         = round(c(brier_glm,       brier_ridge),       3),
  Brier_lo      = round(c(brier_ci_glm[1], brier_ci_ridge[1]), 3),
  Brier_hi      = round(c(brier_ci_glm[2], brier_ci_ridge[2]), 3),
  Cal_intercept = round(c(cal_glm[1],   cal_ridge[1]),   3),
  Cal_slope     = round(c(cal_glm[2],   cal_ridge[2]),   3),
  HL_p          = round(c(hl_glm,       hl_ridge),       4),
  Sens    = c(NA, round(m_ov$sens, 3)),
  Sens_lo = c(NA, round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[1], 3)),
  Sens_hi = c(NA, round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[2], 3)),
  Spec    = c(NA, round(m_ov$spec, 3)),
  Spec_lo = c(NA, round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[1], 3)),
  Spec_hi = c(NA, round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[2], 3)),
  PPV     = c(NA, round(m_ov$ppv,  3)),
  PPV_lo  = c(NA, round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[1], 3)),
  PPV_hi  = c(NA, round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[2], 3)),
  NPV     = c(NA, round(m_ov$npv,  3)),
  NPV_lo  = c(NA, round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[1], 3)),
  NPV_hi  = c(NA, round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[2], 3)),
  stringsAsFactors = FALSE
)
write_csv(overall_tab, file.path(outdir, "overall_GLM_vs_Ridge_Prague.csv"))

# ============================================================
# 2) TREATMENT SUBGROUPS (Ridge)
# ============================================================
cat("=== 2) TREATMENT SUBGROUPS (Ridge) ===\n")

data_testP$treatment_group <- as.character(
  data_testP$ACTIVE_treatment_naivenaive_GCGC_GCCYCRTXGCCYCRTX_remissionrem
)
print(table(data_testP$treatment_group, data_testP$group01, useNA = "ifany"))

# Subgroup function — filters to active_group + remission, then calls core
compute_subgroup <- function(df, active_group, group_label,
                             score_col = "pred_ridge", threshold = 0.5) {
  sub <- df[df$treatment_group %in% c(active_group, "remission") &
              !is.na(df$treatment_group), , drop = FALSE]
  
  y          <- as.integer(sub$group01)
  p          <- as.numeric(sub[[score_col]])
  events     <- sum(y == 1, na.rm = TRUE)
  non_events <- sum(y == 0, na.rm = TRUE)
  
  cat(sprintf("\n  %s vs remission: n=%d | active=%d | remission=%d\n",
              group_label, nrow(sub), non_events, events))
  
  if (events == 0 || non_events == 0) {
    cat("  SKIPPED (one class missing)\n")
    return(NULL)
  }
  
  m <- compute_metrics_core(y, p, threshold)
  
  cat(sprintf("    AUC=%.3f (%.3f-%.3f)  Sens=%.2f  Spec=%.2f  PPV=%.2f  NPV=%.2f\n",
              m$auc_val, m$auc_ci[1], m$auc_ci[3],
              m$sens, m$spec, m$ppv, m$npv))
  
  list(
    roc_obj = m$roc_obj,
    table   = metrics_row(m, id_cols = list(subgroup = group_label))
  )
}

res_naive    <- compute_subgroup(data_testP, "naive",    "Na\u00efve")
res_GC       <- compute_subgroup(data_testP, "GC",       "GC")
res_GCCYCRTX <- compute_subgroup(data_testP, "GCCYCRTX", "GC+CYC+RTX")

write_csv(
  bind_rows(res_naive$table, res_GC$table, res_GCCYCRTX$table),
  file.path(outdir, "subgroup_treatment_Ridge_Prague.csv")
)

# DeLong between subgroups
run_delong <- function(r1, r2, lab) {
  if (is.null(r1) || is.null(r2)) { cat(sprintf("  %s: SKIPPED\n", lab)); return(NULL) }
  dl <- pROC::roc.test(r1, r2, method = "delong")
  cat(sprintf("  %-30s AUC1=%.3f  AUC2=%.3f  p=%.4f\n", lab,
              as.numeric(pROC::auc(r1)), as.numeric(pROC::auc(r2)), dl$p.value))
  dl
}

cat("\n  DeLong tests (naive as reference):\n")
dl_naive_gc       <- run_delong(res_naive$roc_obj, res_GC$roc_obj,       "Naive vs GC")
dl_naive_gccycrtx <- run_delong(res_naive$roc_obj, res_GCCYCRTX$roc_obj, "Naive vs GC+CYC+RTX")
dl_gc_gccycrtx    <- run_delong(res_GC$roc_obj,    res_GCCYCRTX$roc_obj, "GC vs GC+CYC+RTX")

# ROC plot — treatment subgroups
cols_sub <- c("naive" = "black", "GC" = "#ff009c", "GCCYCRTX" = "#00a86b")

svglite(file.path(outdir, "ROC_subgroups_Ridge_Prague.svg"), width = 5.5, height = 5.2)
plot(res_naive$roc_obj, col = cols_sub["naive"], lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "Treatment subgroups \u2014 Berlin Ridge \u2192 Prague [SUPPLEMENT]",
     cex.main = 0.9)
plot(res_GC$roc_obj,       col = cols_sub["GC"],       lwd = 2.2, add = TRUE)
plot(res_GCCYCRTX$roc_obj, col = cols_sub["GCCYCRTX"], lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")

leg_entries <- character(0); leg_cols <- character(0); leg_lwd <- numeric(0); leg_lty <- numeric(0)
if (!is.null(res_naive$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("Na\u00efve \u2014 7-prot panel\n", auc_fmt(res_naive$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["naive"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
if (!is.null(res_GC$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("GC \u2014 7-prot panel\n", auc_fmt(res_GC$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["GC"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
if (!is.null(res_GCCYCRTX$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("GC+CYC+RTX \u2014 7-prot panel\n", auc_fmt(res_GCCYCRTX$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["GCCYCRTX"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
legend("bottomright", legend = leg_entries, col = leg_cols, lwd = leg_lwd, lty = leg_lty,
       bty = "n", cex = 0.65, y.intersp = 1.6)
dev.off()

# ============================================================
# 3) INCREMENTAL MODELS (Ridge protein score)
# ============================================================
cat("\n=== 3) INCREMENTAL MODELS (Ridge protein score) ===\n")

# Fit all models on Berlin training data
mod_prot          <- glm(group01 ~ protein_score,                           data = data_trainB, family = binomial())
mod_prot_crp      <- glm(group01 ~ protein_score + CRP_Plasma_mg_L,         data = data_trainB, family = binomial())
mod_prot_anca     <- glm(group01 ~ protein_score + zANCA,                   data = data_trainB, family = binomial())
mod_prot_crp_anca <- glm(group01 ~ protein_score + CRP_Plasma_mg_L + zANCA, data = data_trainB, family = binomial())
mod_crp           <- glm(group01 ~ CRP_Plasma_mg_L,                         data = data_trainB, family = binomial())
mod_anca          <- glm(group01 ~ zANCA,                                   data = data_trainB, family = binomial())

# Score Prague test set
data_testP$prob_prot          <- predict(mod_prot,          newdata = data_testP, type = "response")
data_testP$prob_prot_crp      <- predict(mod_prot_crp,      newdata = data_testP, type = "response")
data_testP$prob_prot_anca     <- predict(mod_prot_anca,     newdata = data_testP, type = "response")
data_testP$prob_prot_crp_anca <- predict(mod_prot_crp_anca, newdata = data_testP, type = "response")
data_testP$prob_crp           <- predict(mod_crp,           newdata = data_testP, type = "response")
data_testP$prob_anca          <- predict(mod_anca,          newdata = data_testP, type = "response")

# Compute metrics for all six models via the unified core
inc_models <- list(
  list(col = "prob_prot",          label = "Ridge 7-protein"),
  list(col = "prob_crp",           label = "CRP alone"),
  list(col = "prob_anca",          label = "ANCA alone"),
  list(col = "prob_prot_crp",      label = "7-protein + CRP"),
  list(col = "prob_prot_anca",     label = "7-protein + ANCA"),
  list(col = "prob_prot_crp_anca", label = "7-protein + CRP + ANCA")
)

inc_results <- lapply(inc_models, function(md) {
  m <- compute_metrics_core(truth, data_testP[[md$col]])
  list(roc_obj = m$roc_obj,
       table   = metrics_row(m, id_cols = list(model = md$label)))
})
names(inc_results) <- sapply(inc_models, `[[`, "label")

inc_table <- bind_rows(lapply(inc_results, `[[`, "table"))
print(inc_table[, c("model", "n", "events",
                    "AUC",  "AUC_lo",  "AUC_hi",
                    "Sens", "Sens_lo", "Sens_hi",
                    "Spec", "Spec_lo", "Spec_hi",
                    "PPV",  "PPV_lo",  "PPV_hi",
                    "NPV",  "NPV_lo",  "NPV_hi")])
write_csv(inc_table, file.path(outdir, "incremental_models_Ridge_Prague.csv"))
# Convenience ROC objects for DeLong + plotting
roc_prot          <- inc_results[["Ridge 7-protein"]]$roc_obj
roc_prot_crp      <- inc_results[["7-protein + CRP"]]$roc_obj
roc_prot_anca     <- inc_results[["7-protein + ANCA"]]$roc_obj
roc_prot_crp_anca <- inc_results[["7-protein + CRP + ANCA"]]$roc_obj
roc_crp           <- inc_results[["CRP alone"]]$roc_obj
roc_anca          <- inc_results[["ANCA alone"]]$roc_obj

auc_prot          <- as.numeric(pROC::auc(roc_prot))
auc_prot_crp      <- as.numeric(pROC::auc(roc_prot_crp))
auc_prot_anca     <- as.numeric(pROC::auc(roc_prot_anca))
auc_prot_crp_anca <- as.numeric(pROC::auc(roc_prot_crp_anca))

# DeLong tests — Ridge protein as reference
del_prot_vs_crp       <- pROC::roc.test(roc_prot, roc_crp,           method = "delong")
del_prot_vs_protCRP   <- pROC::roc.test(roc_prot, roc_prot_crp,      method = "delong")
del_prot_vs_protANCA  <- pROC::roc.test(roc_prot, roc_prot_anca,     method = "delong")
del_prot_vs_full      <- pROC::roc.test(roc_prot, roc_prot_crp_anca, method = "delong")

cat("\n=== DeLong: Ridge protein vs other models (Prague) ===\n")
cat(sprintf("  Ridge protein vs CRP alone:        AUC_prot=%.3f | AUC_CRP=%.3f  | p=%.4f\n",
            auc_prot, as.numeric(pROC::auc(roc_crp)), del_prot_vs_crp$p.value))
cat(sprintf("  Ridge protein vs protein+CRP:      AUC_prot=%.3f | AUC+CRP=%.3f  | p=%.4f\n",
            auc_prot, auc_prot_crp,      del_prot_vs_protCRP$p.value))
cat(sprintf("  Ridge protein vs protein+ANCA:     AUC_prot=%.3f | AUC+ANCA=%.3f | p=%.4f\n",
            auc_prot, auc_prot_anca,     del_prot_vs_protANCA$p.value))
cat(sprintf("  Ridge protein vs protein+CRP+ANCA: AUC_prot=%.3f | AUC_full=%.3f | p=%.4f\n",
            auc_prot, auc_prot_crp_anca, del_prot_vs_full$p.value))

writeLines(c(
  "DeLong incremental \u2014 Prague, Ridge protein score [SUPPLEMENT]",
  sprintf("  Ridge protein vs CRP alone:        \u0394AUC=%+.4f p=%.4f",
          as.numeric(pROC::auc(roc_crp)) - auc_prot, del_prot_vs_crp$p.value),
  sprintf("  Ridge protein vs protein+CRP:      \u0394AUC=%+.4f p=%.4f",
          auc_prot_crp      - auc_prot, del_prot_vs_protCRP$p.value),
  sprintf("  Ridge protein vs protein+ANCA:     \u0394AUC=%+.4f p=%.4f",
          auc_prot_anca     - auc_prot, del_prot_vs_protANCA$p.value),
  sprintf("  Ridge protein vs protein+CRP+ANCA: \u0394AUC=%+.4f p=%.4f",
          auc_prot_crp_anca - auc_prot, del_prot_vs_full$p.value)
), file.path(outdir, "delong_incremental_Ridge_Prague.txt"))

# ROC plot — incremental models
cols_inc <- c(
  "Ridge 7-protein"        = "mediumpurple2",
  "7-protein + CRP"        = "#d95f02",
  "7-protein + ANCA"       = "#e7298a",
  "7-protein + CRP + ANCA" = "#66a61e"
)

svglite(file.path(outdir, "ROC_incremental_Ridge_Prague.svg"), width = 6, height = 5.2)
plot(roc_prot,          col = cols_inc["Ridge 7-protein"],        lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "Incremental value \u2014 Berlin Ridge \u2192 Prague [SUPPLEMENT]",
     cex.main = 0.9)
plot(roc_prot_crp,      col = cols_inc["7-protein + CRP"],        lwd = 2.2, add = TRUE)
plot(roc_prot_anca,     col = cols_inc["7-protein + ANCA"],       lwd = 2.2, add = TRUE)
plot(roc_prot_crp_anca, col = cols_inc["7-protein + CRP + ANCA"], lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")
legend("bottomright",
       legend = c(
         paste0("Ridge 7-protein\n",    auc_fmt(roc_prot)),
         paste0("Ridge + CRP\n",        auc_fmt(roc_prot_crp)),
         paste0("Ridge + ANCA\n",       auc_fmt(roc_prot_anca)),
         paste0("Ridge + CRP + ANCA\n", auc_fmt(roc_prot_crp_anca))
       ),
       col = cols_inc, lwd = 2.2, lty = 1,
       bty = "n", cex = 0.65, y.intersp = 1.6)
dev.off()

# ============================================================
# 4) eGFR STRATUM ANALYSIS
# ============================================================
cat("\n=== 4) eGFR STRATUM ANALYSIS ===\n")

data_testP$eGFR_group <- ifelse(data_testP$eGFR_2009_higher45_YES1 == 1, ">45", "\u226445")

cat("eGFR group \u00d7 outcome:\n")
print(table(data_testP$eGFR_group, data_testP$group01))

# Stratum function — same pattern as compute_subgroup, adds min-events guard
compute_metrics_stratum <- function(df, stratum_label,
                                    score_col = "pred_ridge", threshold = 0.5) {
  y          <- as.integer(df$group01)
  p          <- as.numeric(df[[score_col]])
  events     <- sum(y == 1, na.rm = TRUE)
  non_events <- sum(y == 0, na.rm = TRUE)
  
  if (events < 5 || non_events < 5) {
    cat(sprintf("  %s: SKIPPED (events=%d, non-events=%d)\n",
                stratum_label, events, non_events))
    return(NULL)
  }
  
  m <- compute_metrics_core(y, p, threshold)
  
  cat(sprintf("  %s (n=%d): AUC=%.3f (%.3f-%.3f)  Sens=%.2f  Spec=%.2f  PPV=%.2f  NPV=%.2f\n",
              stratum_label, length(y),
              m$auc_val, m$auc_ci[1], m$auc_ci[3],
              m$sens, m$spec, m$ppv, m$npv))
  
  list(
    roc_obj = m$roc_obj,
    table   = metrics_row(m, id_cols = list(stratum = stratum_label))
  )
}

res_egfr_low  <- compute_metrics_stratum(subset(data_testP, eGFR_group == "\u226445"), "\u226445")
res_egfr_high <- compute_metrics_stratum(subset(data_testP, eGFR_group == ">45"),     ">45")

egfr_tab <- bind_rows(res_egfr_low$table, res_egfr_high$table)
write_csv(egfr_tab, file.path(outdir, "eGFR_stratum_Ridge_Prague.csv"))

# DeLong between strata
if (!is.null(res_egfr_low) && !is.null(res_egfr_high)) {
  dl_egfr <- pROC::roc.test(res_egfr_low$roc_obj, res_egfr_high$roc_obj, method = "delong")
  cat(sprintf("\n  DeLong \u226445 vs >45: p=%.4f\n", dl_egfr$p.value))
  capture.output(dl_egfr, file = file.path(outdir, "delong_eGFR_Ridge_Prague.txt"))
}

# ROC plot — eGFR strata
svglite(file.path(outdir, "ROC_eGFR_Ridge_Prague.svg"), width = 5.5, height = 5.2)
plot(res_egfr_low$roc_obj,
     col = "#34495E", lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "eGFR strata \u2014 Berlin Ridge \u2192 Prague [SUPPLEMENT]",
     cex.main = 0.9)
plot(res_egfr_high$roc_obj, col = "#D35400", lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")
legend("bottomright",
       legend = c(
         paste0("eGFR \u226445 mL/min/1.73 m\u00b2 (n=", res_egfr_low$table$n,  ")\n", auc_fmt(res_egfr_low$roc_obj)),
         paste0("eGFR >45 mL/min/1.73 m\u00b2 (n=",      res_egfr_high$table$n, ")\n", auc_fmt(res_egfr_high$roc_obj))
       ),
       col = c("#34495E", "#D35400"), lwd = 2.2, lty = 1,
       bty = "n", cex = 0.65, y.intersp = 1.6)
dev.off()

# ============================================================
# 5) REMISSION ON vs OFF TREATMENT (descriptive — no model)
# ============================================================
cat("\n=== 5) REMISSION ON vs OFF TREATMENT ===\n")

df_rem <- dfzt %>%
  filter(diseasestatus == "remission") %>%
  filter(!is.na(Remission_OFF_treatment_last_6month_YES1))

cat("Cohort \u00d7 treatment status:\n")
print(table(df_rem$cohort, df_rem$Remission_OFF_treatment_last_6month_YES1, useNA = "ifany"))

desc_list <- lapply(prots_7, function(v) {
  sub   <- df_rem[!is.na(df_rem[[v]]), ]
  g_off <- sub[sub$Remission_OFF_treatment_last_6month_YES1 == "1", ][[v]]
  g_on  <- sub[sub$Remission_OFF_treatment_last_6month_YES1 == "0", ][[v]]
  if (length(g_off) < 3 || length(g_on) < 3) return(NULL)
  p_w <- wilcox.test(g_off, g_on)$p.value
  d   <- (mean(g_off) - mean(g_on)) /
    sqrt(((length(g_off) - 1) * var(g_off) + (length(g_on) - 1) * var(g_on)) /
           (length(g_off) + length(g_on) - 2))
  data.frame(
    protein  = v,
    n_OFF    = length(g_off), med_OFF = round(median(g_off), 3), IQR_OFF = round(IQR(g_off), 3),
    n_ON     = length(g_on),  med_ON  = round(median(g_on),  3), IQR_ON  = round(IQR(g_on),  3),
    p_wilcox = round(p_w, 4), cohens_d = round(d, 3)
  )
})
desc_rem <- bind_rows(desc_list)
print(desc_rem)
write_csv(desc_rem, file.path(outdir, "Remission_ON_OFF_descriptives.csv"))

# Linear models adjusting for cohort
lm_tab <- bind_rows(lapply(prots_7, function(v) {
  sub <- df_rem[!is.na(df_rem[[v]]) & !is.na(df_rem$cohort), ]
  if (nrow(sub) < 10) return(NULL)
  fit <- lm(sub[[v]] ~ Remission_OFF_treatment_last_6month_YES1 + cohort, data = sub)
  s   <- summary(fit)
  data.frame(protein  = v,
             term     = rownames(coef(s)),
             estimate = round(coef(s)[, "Estimate"],    3),
             p        = round(coef(s)[, "Pr(>|t|)"],    4))
}))
write_csv(lm_tab, file.path(outdir, "Remission_ON_OFF_linear_models.csv"))

# ============================================================
# 6) SUGGESTED RESULTS TEXT
# ============================================================
cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SUGGESTED RESULTS TEXT (fill in [X] with actual numbers)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(strwrap(paste0(
  "Supplementary sensitivity analysis \u2014 Ridge-penalised model. ",
  "To confirm that the calibration failure of the primary GLM in Prague was ",
  "attributable to quasi-complete separation rather than the protein panel itself, ",
  "we fitted a Ridge-penalised logistic regression (L2 penalty; lambda selected by ",
  "10-fold cross-validation on Berlin data only) and applied it without modification ",
  "to the Prague cohort. ",
  sprintf("The Ridge model achieved AUC=%.3f (95%% CI %.3f\u2013%.3f) in Prague, ",
          auc_ridge, ci_ridge[1], ci_ridge[3]),
  sprintf("indistinguishable from the primary GLM (AUC=%.3f; DeLong p=%.3f), ",
          auc_glm, delong_glm_ridge$p.value),
  sprintf("but with substantially improved calibration (slope=%.2f, ", cal_ridge[2]),
  sprintf("Brier=%.3f (95%% CI %.3f\u2013%.3f), Hosmer\u2013Lemeshow p=%.3f). ",
          brier_ridge, brier_ci_ridge[1], brier_ci_ridge[2], hl_ridge),
  sprintf("At a decision threshold of 0.5, the Ridge model yielded sensitivity=%.2f (95%% CI %.2f\u2013%.2f), ",
          m_ov$sens, binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[2]),
  sprintf("specificity=%.2f (95%% CI %.2f\u2013%.2f), ",
          m_ov$spec, binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[2]),
  sprintf("PPV=%.2f (95%% CI %.2f\u2013%.2f), ",
          m_ov$ppv, binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[2]),
  sprintf("and NPV=%.2f (95%% CI %.2f\u2013%.2f). ",
          m_ov$npv, binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[2]),
  "Adding ANCA titres or CRP to the Ridge protein score did not significantly ",
  sprintf("improve discrimination (protein+ANCA: AUC=%.3f, DeLong p=%.3f; ",
          auc_prot_anca, del_prot_vs_protANCA$p.value),
  sprintf("protein+CRP: AUC=%.3f, DeLong p=%.3f). ",
          auc_prot_crp, del_prot_vs_protCRP$p.value),
  "Discrimination was consistent across treatment subgroups ",
  sprintf("(Na\u00efve: AUC=%.3f [%.3f\u2013%.3f]; ",
          res_naive$table$AUC, res_naive$table$AUC_lo, res_naive$table$AUC_hi),
  sprintf("GC: AUC=%.3f [%.3f\u2013%.3f]; ",
          res_GC$table$AUC, res_GC$table$AUC_lo, res_GC$table$AUC_hi),
  sprintf("GC+CYC+RTX: AUC=%.3f [%.3f\u2013%.3f]; all pairwise DeLong p>0.05) ",
          res_GCCYCRTX$table$AUC, res_GCCYCRTX$table$AUC_lo, res_GCCYCRTX$table$AUC_hi),
  sprintf("and across eGFR strata (\u226445: AUC=%.3f [%.3f\u2013%.3f]; ",
          res_egfr_low$table$AUC, res_egfr_low$table$AUC_lo, res_egfr_low$table$AUC_hi),
  sprintf(">45: AUC=%.3f [%.3f\u2013%.3f]; DeLong p=%.3f). ",
          res_egfr_high$table$AUC, res_egfr_high$table$AUC_lo, res_egfr_high$table$AUC_hi,
          dl_egfr$p.value),
  "We note that the Ridge model was retained after inspecting its Prague calibration; ",
  "Prague therefore cannot be considered a fully independent validation for this model, ",
  "and these results should be interpreted as supportive rather than confirmatory."
), width = 90), sep = "\n")

message("\n\u2713 Ridge Berlin \u2192 Prague analysis complete.")
message("  Outputs in: ", outdir)

