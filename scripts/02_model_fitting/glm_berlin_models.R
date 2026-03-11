# ============================================================
# GLM BERLIN → PRAGUE — Complete analysis script
# ============================================================
# PRIMARY MODEL: 7-protein GLM trained on Berlin cohort,
# validated externally on Prague cohort.
# ============================================================

library(pROC)
library(dplyr)
library(readr)
library(ResourceSelection)
library(boot)
library(svglite)

set.seed(7)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
data_csv     <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
outdir       <- "../proteomics_ml/data/PragueTQL/final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prots_7 <- c("AHSG_001", "CLEC3B_020", "COMP_023",
             "F9_027", "LRG1_042", "MCAM_047", "MRC1_073..")

# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
dfzt        <- read_csv(data_csv)
data_trainB <- dfzt %>% filter(cohort == "Berlin")
data_testP  <- dfzt %>% filter(cohort == "Prague")

message("Berlin n=", nrow(data_trainB),
        "; Prague n=", nrow(data_testP),
        " (events in Prague = ", sum(data_testP$group01 == 1, na.rm = TRUE), ")")

# ----------------------------------------------------------
# LOAD FROZEN MODEL & SCORE PRAGUE
# ----------------------------------------------------------
model_7_berlin <- readRDS(file.path(artifact_dir, "model_7panel_Berlin_glm.rds"))

data_testP$pred_berlin <- predict(model_7_berlin,
                                  newdata = data_testP,
                                  type    = "response")

cat("NAs in pred_berlin:", sum(is.na(data_testP$pred_berlin)), "\n")

# ----------------------------------------------------------
# SHARED HELPERS
# ----------------------------------------------------------

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

# Core: AUC + threshold metrics for any truth/prob pair
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

# Turns compute_metrics_core output into one tidy data.frame row
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

truth     <- as.integer(data_testP$group01)
pred_glm  <- as.numeric(data_testP$pred_berlin)
n         <- length(truth)
threshold <- 0.5

roc_glm_overall <- pROC::roc(truth, pred_glm, quiet = TRUE)
auc_glm  <- as.numeric(pROC::auc(roc_glm_overall))
ci_glm   <- as.numeric(pROC::ci.auc(roc_glm_overall, method = "delong"))

cal_glm  <- coef(glm(truth ~ logit_safe(pred_glm), family = binomial()))
hl_glm   <- tryCatch(ResourceSelection::hoslem.test(truth, pred_glm, g = 10)$p.value,
                     error = function(e) NA)
brier_glm <- mean((pred_glm - truth)^2)
brier_ci_glm <- { b <- boot::boot(data.frame(p = pred_glm, y = truth),
                                  function(d, i) mean((d$p[i] - d$y[i])^2), R = 2000)
round(quantile(b$t, c(.025, .975)), 3) }

cat(sprintf("  GLM AUC=%.3f (%.3f-%.3f) | slope=%.3f | intercept=%.3f | Brier=%.3f (%.3f-%.3f) | HL p=%.4f\n",
            auc_glm, ci_glm[1], ci_glm[3],
            cal_glm[2], cal_glm[1],
            brier_glm, brier_ci_glm[1], brier_ci_glm[2],
            hl_glm))

# Threshold metrics
m_ov <- compute_metrics_core(truth, pred_glm, threshold)

cat(sprintf("\n  Threshold metrics at %.1f:\n", threshold))
cat(sprintf("  n=%d | TP=%d TN=%d FP=%d FN=%d\n", n, m_ov$TP, m_ov$TN, m_ov$FP, m_ov$FN))
cat(sprintf("  Sensitivity : %.3f (%.3f-%.3f)\n", m_ov$sens,
            binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[2]))
cat(sprintf("  Specificity : %.3f (%.3f-%.3f)\n", m_ov$spec,
            binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[2]))
cat(sprintf("  PPV         : %.3f (%.3f-%.3f)\n", m_ov$ppv,
            binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[1], binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[2]))
cat(sprintf("  NPV         : %.3f (%.3f-%.3f)\n\n", m_ov$npv,
            binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[1], binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[2]))

# Save overall table
overall_tab <- data.frame(
  model         = "Berlin GLM",
  n             = n,
  events        = sum(truth),
  AUC           = round(auc_glm,          3),
  AUC_lo        = round(ci_glm[1],        3),
  AUC_hi        = round(ci_glm[3],        3),
  Brier         = round(brier_glm,        3),
  Brier_lo      = round(brier_ci_glm[1],  3),
  Brier_hi      = round(brier_ci_glm[2],  3),
  Cal_intercept = round(cal_glm[1],       3),
  Cal_slope     = round(cal_glm[2],       3),
  HL_p          = round(hl_glm,           4),
  Sens    = round(m_ov$sens, 3),
  Sens_lo = round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[1], 3),
  Sens_hi = round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FN)[2], 3),
  Spec    = round(m_ov$spec, 3),
  Spec_lo = round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[1], 3),
  Spec_hi = round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FP)[2], 3),
  PPV     = round(m_ov$ppv,  3),
  PPV_lo  = round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[1], 3),
  PPV_hi  = round(binom_ci(m_ov$TP, m_ov$TP + m_ov$FP)[2], 3),
  NPV     = round(m_ov$npv,  3),
  NPV_lo  = round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[1], 3),
  NPV_hi  = round(binom_ci(m_ov$TN, m_ov$TN + m_ov$FN)[2], 3),
  TP = m_ov$TP, TN = m_ov$TN, FP = m_ov$FP, FN = m_ov$FN,
  stringsAsFactors = FALSE
)
write_csv(overall_tab, file.path(outdir, "overall_GLM_Berlin_Prague.csv"))

# ============================================================
# 2) TREATMENT SUBGROUPS (GLM)
# ============================================================
cat("=== 2) TREATMENT SUBGROUPS (GLM) ===\n")

data_testP$treatment_group <- as.character(
  data_testP$ACTIVE_treatment_naivenaive_GCGC_GCCYCRTXGCCYCRTX_remissionrem
)

cat("Treatment subgroup x outcome counts (Prague):\n")
print(table(data_testP$treatment_group, data_testP$group01, useNA = "ifany"))

compute_subgroup <- function(df, active_group, group_label,
                             score_col = "pred_berlin", threshold = 0.5) {
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

subgroup_tab <- bind_rows(res_naive$table, res_GC$table, res_GCCYCRTX$table)
write_csv(subgroup_tab, file.path(outdir, "subgroup_treatment_GLM_Prague.csv"))

# DeLong between subgroups
run_delong <- function(r1, r2, lab, boot_n = 2000) {
  if (is.null(r1) || is.null(r2)) {
    cat(sprintf("  %s: SKIPPED\n", lab)); return(NULL)
  }
  
  # DeLong test for p-value
  dl    <- pROC::roc.test(r1, r2, method = "delong")
  delta <- as.numeric(pROC::auc(r1)) - as.numeric(pROC::auc(r2))
  
  # Manual bootstrap CI for the AUC difference
  # r1 and r2 must come from the same cases (paired) or different (unpaired)
  # We extract the original data from the roc objects
  cases1    <- r1$cases;    controls1 <- r1$controls
  cases2    <- r2$cases;    controls2 <- r2$controls
  
  set.seed(7)
  boot_deltas <- replicate(boot_n, {
    # resample cases and controls separately (stratified)
    c1_b <- sample(cases1,    length(cases1),    replace = TRUE)
    ct1_b <- sample(controls1, length(controls1), replace = TRUE)
    c2_b <- sample(cases2,    length(cases2),    replace = TRUE)
    ct2_b <- sample(controls2, length(controls2), replace = TRUE)
    
    auc1_b <- pROC::auc(pROC::roc(
      response  = c(rep(1, length(c1_b)),  rep(0, length(ct1_b))),
      predictor = c(c1_b, ct1_b), quiet = TRUE))
    auc2_b <- pROC::auc(pROC::roc(
      response  = c(rep(1, length(c2_b)),  rep(0, length(ct2_b))),
      predictor = c(c2_b, ct2_b), quiet = TRUE))
    
    as.numeric(auc1_b) - as.numeric(auc2_b)
  })
  
  ci_lo <- quantile(boot_deltas, 0.025)
  ci_hi <- quantile(boot_deltas, 0.975)
  
  cat(sprintf("  %-35s AUC1=%.3f  AUC2=%.3f  ΔAUC=%+.3f (95%% CI %+.3f – %+.3f)  p=%.4f\n",
              lab,
              as.numeric(pROC::auc(r1)),
              as.numeric(pROC::auc(r2)),
              delta, ci_lo, ci_hi, dl$p.value))
  
  list(D     = as.numeric(dl$statistic),
       df    = as.numeric(dl$parameter),
       p     = dl$p.value,
       delta = delta,
       ci_lo = as.numeric(ci_lo),
       ci_hi = as.numeric(ci_hi),
       auc1  = as.numeric(pROC::auc(r1)),
       auc2  = as.numeric(pROC::auc(r2)))
}
cat("\n  DeLong tests (naive as reference):\n")
dl_naive_gc       <- run_delong(res_naive$roc_obj, res_GC$roc_obj,       "Naive vs GC")
dl_naive_gccycrtx <- run_delong(res_naive$roc_obj, res_GCCYCRTX$roc_obj, "Naive vs GC+CYC+RTX")
dl_gc_gccycrtx    <- run_delong(res_GC$roc_obj,    res_GCCYCRTX$roc_obj, "GC vs GC+CYC+RTX")



# ROC plot — treatment subgroups
cols_sub <- c("naive" = "black", "GC" = "#ff009c", "GCCYCRTX" = "#00a86b")

svglite(file.path(outdir, "ROC_subgroups_GLM_Prague.svg"), width = 5.5, height = 5.2)
plot(res_naive$roc_obj, col = cols_sub["naive"], lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "Treatment subgroups \u2014 Berlin GLM \u2192 Prague",
     cex.main = 0.9)
plot(res_GC$roc_obj,       col = cols_sub["GC"],       lwd = 2.2, add = TRUE)
plot(res_GCCYCRTX$roc_obj, col = cols_sub["GCCYCRTX"], lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")

leg_entries <- character(0); leg_cols <- character(0)
leg_lwd <- numeric(0);       leg_lty <- numeric(0)
if (!is.null(res_naive$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("Na\u00efve\n", auc_fmt(res_naive$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["naive"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
if (!is.null(res_GC$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("GC\n", auc_fmt(res_GC$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["GC"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
if (!is.null(res_GCCYCRTX$roc_obj)) {
  leg_entries <- c(leg_entries, paste0("GC+CYC+RTX\n", auc_fmt(res_GCCYCRTX$roc_obj)))
  leg_cols <- c(leg_cols, cols_sub["GCCYCRTX"]); leg_lwd <- c(leg_lwd, 2.2); leg_lty <- c(leg_lty, 1)
}
legend("bottomright", legend = leg_entries, col = leg_cols, lwd = leg_lwd, lty = leg_lty,
       bty = "n", cex = 0.65, y.intersp = 1.6)
dev.off()

# ============================================================
# 3) INCREMENTAL MODELS (GLM protein score vs CRP / ANCA)
# ============================================================
cat("\n=== 3) INCREMENTAL MODELS (GLM protein score) ===\n")

# Derive protein_score from the frozen Berlin GLM
data_trainB$protein_score <- as.numeric(predict(model_7_berlin,
                                                newdata = data_trainB,
                                                type    = "response"))
data_testP$protein_score  <- as.numeric(predict(model_7_berlin,
                                                newdata = data_testP,
                                                type    = "response"))

# Fit augmented models on Berlin training data
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

# Compute full metrics for all six models
inc_models <- list(
  list(col = "prob_prot",          label = "GLM 7-protein"),
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

cat("\n  Incremental model metrics (Prague):\n")
print(inc_table[, c("model", "n", "events",
                    "AUC",  "AUC_lo",  "AUC_hi",
                    "Sens", "Sens_lo", "Sens_hi",
                    "Spec", "Spec_lo", "Spec_hi",
                    "PPV",  "PPV_lo",  "PPV_hi",
                    "NPV",  "NPV_lo",  "NPV_hi")])

write_csv(inc_table, file.path(outdir, "incremental_models_GLM_Prague.csv"))

# Convenience ROC objects for DeLong + plotting
roc_prot          <- inc_results[["GLM 7-protein"]]$roc_obj
roc_prot_crp      <- inc_results[["7-protein + CRP"]]$roc_obj
roc_prot_anca     <- inc_results[["7-protein + ANCA"]]$roc_obj
roc_prot_crp_anca <- inc_results[["7-protein + CRP + ANCA"]]$roc_obj
roc_crp           <- inc_results[["CRP alone"]]$roc_obj
roc_anca          <- inc_results[["ANCA alone"]]$roc_obj

auc_prot          <- as.numeric(pROC::auc(roc_prot))
auc_prot_crp      <- as.numeric(pROC::auc(roc_prot_crp))
auc_prot_anca     <- as.numeric(pROC::auc(roc_prot_anca))
auc_prot_crp_anca <- as.numeric(pROC::auc(roc_prot_crp_anca))
auc_anca <- as.numeric(pROC::auc(roc_anca))

# DeLong tests — GLM protein as reference
# del_prot_vs_crp       <- pROC::roc.test(roc_prot, roc_crp,           method = "delong")
# del_prot_vs_protCRP   <- pROC::roc.test(roc_prot, roc_prot_crp,      method = "delong")
# del_prot_vs_protANCA  <- pROC::roc.test(roc_prot, roc_prot_anca,     method = "delong")
# del_prot_vs_full      <- pROC::roc.test(roc_prot, roc_prot_crp_anca, method = "delong")

dl_prot_vs_crp       <- run_delong(roc_prot, roc_crp,           "GLM protein vs CRP alone")
dl_prot_vs_anca       <- run_delong(roc_prot, roc_anca,           "GLM protein vs ANCA alone")
dl_prot_vs_protCRP   <- run_delong(roc_prot, roc_prot_crp,      "GLM protein vs protein+CRP")
dl_prot_vs_protANCA  <- run_delong(roc_prot, roc_prot_anca,     "GLM protein vs protein+ANCA")
dl_prot_vs_full      <- run_delong(roc_prot, roc_prot_crp_anca, "GLM protein vs protein+CRP+ANCA")


# ROC plot — incremental models
cols_inc <- c(
  "GLM 7-protein"          = "#C70500",
  "7-protein + ANCA"       = "#6094ED"
)

svglite(file.path(outdir, "ROC_incremental_GLM_Prague_pro_ANCA.svg"), width = 6, height = 5.2)
plot(roc_prot,          col = cols_inc["GLM 7-protein"],          lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "Incremental value \u2014 Berlin GLM \u2192 Prague",
     cex.main = 0.9)
plot(roc_prot_anca,     col = cols_inc["7-protein + ANCA"],       lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")
legend("bottomright",
       legend = c(
         paste0("GLM 7-protein\n",     auc_fmt(roc_prot)),
         paste0("GLM + ANCA\n",        auc_fmt(roc_prot_anca))
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

compute_metrics_stratum <- function(df, stratum_label,
                                    score_col = "pred_berlin", threshold = 0.5) {
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
write_csv(egfr_tab, file.path(outdir, "eGFR_stratum_GLM_Prague.csv"))

# DeLong between strata
if (!is.null(res_egfr_low) && !is.null(res_egfr_high)) {
  dl_egfr <- run_delong(res_egfr_low$roc_obj, res_egfr_high$roc_obj, "eGFR ≤45 vs >45")
}

# ROC plot — eGFR strata
svglite(file.path(outdir, "ROC_eGFR_GLM_Prague.svg"), width = 5.5, height = 5.2)
plot(res_egfr_low$roc_obj,
     col = "#34495E", lwd = 2.2,
     xlab = "1 \u2013 Specificity", ylab = "Sensitivity",
     main = "eGFR strata \u2014 Berlin GLM \u2192 Prague",
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
write_csv(desc_rem, file.path(outdir, "Remission_ON_OFF_descriptives_GLM.csv"))

# Linear models adjusting for cohort
lm_tab <- bind_rows(lapply(prots_7, function(v) {
  sub <- df_rem[!is.na(df_rem[[v]]) & !is.na(df_rem$cohort), ]
  if (nrow(sub) < 10) return(NULL)
  fit <- lm(sub[[v]] ~ Remission_OFF_treatment_last_6month_YES1 + cohort, data = sub)
  s   <- summary(fit)
  data.frame(protein  = v,
             term     = rownames(coef(s)),
             estimate = round(coef(s)[, "Estimate"], 3),
             p        = round(coef(s)[, "Pr(>|t|)"], 4))
}))
write_csv(lm_tab, file.path(outdir, "Remission_ON_OFF_linear_models_GLM.csv"))


delong_to_row <- function(dl, comparison, type) {
  if (is.null(dl)) return(NULL)
  data.frame(
    comparison = comparison,
    type       = type,
    auc1       = as.numeric(dl$auc1)[1],
    auc2       = as.numeric(dl$auc2)[1],
    delta      = as.numeric(dl$delta)[1],
    ci_lo      = as.numeric(dl$ci_lo)[1],
    ci_hi      = as.numeric(dl$ci_hi)[1],
    D          = as.numeric(dl$D)[1],
    df         = as.numeric(dl$df)[1],
    p          = as.numeric(dl$p)[1],
    stringsAsFactors = FALSE
  )
}

delong_all <- bind_rows(
  delong_to_row(dl_naive_gc,        "Naive vs GC",           "treatment"),
  delong_to_row(dl_naive_gccycrtx,  "Naive vs GC+CYC+RTX",  "treatment"),
  delong_to_row(dl_gc_gccycrtx,     "GC vs GC+CYC+RTX",     "treatment"),
  delong_to_row(dl_egfr,            "eGFR ≤45 vs >45",       "eGFR"),
  delong_to_row(dl_prot_vs_anca,     "Protein vs ANCA alone",  "incremental"),
  delong_to_row(dl_prot_vs_crp,     "Protein vs CRP alone",  "incremental"),
  delong_to_row(dl_prot_vs_protCRP, "Protein vs +CRP",       "incremental"),
  delong_to_row(dl_prot_vs_protANCA,"Protein vs +ANCA",      "incremental"),
  delong_to_row(dl_prot_vs_full,    "Protein vs +CRP+ANCA",  "incremental")
)

write_csv(delong_all, file.path(outdir, "delong_all_comparisons_GLM_Prague.csv"))

