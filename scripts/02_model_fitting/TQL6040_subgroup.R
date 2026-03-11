# ============================================================
# TQL 60/40 SPLIT — MODEL COMPARISON ANALYSIS
# 7-protein panel vs CRP, ANCA, protein+CRP, protein+ANCA,
# protein+CRP+ANCA; subgroup analyses (MPO/PR3)
# ============================================================
# Requirements: model_7panel_Train60_glm.rds must already be
# saved (from TQL_models.R). train60 and test40 patient lists
# must be loadable from their .r files.
# ============================================================

library(pROC)
library(ggplot2)
library(svglite)
library(dplyr)
library(tidyr)
library(readr)

set.seed(7)  # report in Methods

# PATHS — adjust to your directory layout
# ----------------------------------------------------------
data_csv      <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
model_dir     <- "../proteomics_ml/data/models_and_artifacts/final"
train60_file  <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test40_file   <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"
outdir        <- "../proteomics_ml/data/TQL_60_40_ModelComparisons/finalx"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
dfz <- read_csv(data_csv)

load(train60_file)   # provides full_meta_train
load(test40_file)    # provides full_meta_test

train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)

cat("Train60 n =", nrow(train60), "| events =", sum(train60$group01 == 1), "\n")
cat("Test40  n =", nrow(test40),  "| events =", sum(test40$group01  == 1), "\n")

# ----------------------------------------------------------
# PREDICTORS
# ----------------------------------------------------------
prots_7 <- c("AHSG_001","CLEC3B_020","COMP_023",
             "F9_027","LRG1_042","MCAM_047","MRC1_073..")

# Sanity-check columns exist
stopifnot(all(prots_7 %in% colnames(dfz)))
stopifnot("CRP_Plasma_mg_L" %in% colnames(dfz))
stopifnot("zANCA"           %in% colnames(dfz))
stopifnot("AAV_group"       %in% colnames(dfz))

# ----------------------------------------------------------
# LOAD FROZEN 7-PANEL TRAIN60 MODEL
# ----------------------------------------------------------
model_train60_7 <- readRDS(file.path(model_dir, "model_7panel_Train60_glm.rds"))

# ── Derive protein_score (predicted probability from 7-panel)
# in BOTH train60 and test40 — needed as input to incremental models
train60$protein_score <- predict(model_train60_7,
                                 newdata = train60,
                                 type    = "response")
test40$protein_score  <- predict(model_train60_7,
                                 newdata = test40,
                                 type    = "response")

# ============================================================
# SECTION 1: FIT ALL COMPARISON MODELS ON TRAIN60
# ============================================================
# Note: protein_score is the frozen 7-panel score used as a
# single linear predictor, then augmented with CRP / zANCA.
# This mirrors exactly how the Berlin analysis was done.

# 1a) Standalone reference models (single predictor, trained on train60)
mod_crp  <- glm(group01 ~ CRP_Plasma_mg_L, data = train60, family = binomial())
mod_anca <- glm(group01 ~ zANCA,           data = train60, family = binomial())

# 1b) 7-protein panel alone (re-expressed as a single score input)
mod_prot <- glm(group01 ~ protein_score,   data = train60, family = binomial())

# 1c) Augmented models (protein_score + clinical marker)
mod_prot_crp      <- glm(group01 ~ protein_score + CRP_Plasma_mg_L,
                         data = train60, family = binomial())
mod_prot_anca     <- glm(group01 ~ protein_score + zANCA,
                         data = train60, family = binomial())
mod_prot_crp_anca <- glm(group01 ~ protein_score + CRP_Plasma_mg_L + zANCA,
                         data = train60, family = binomial())

# 1d) CRP + ANCA combined (no protein panel) — standalone clinical comparator
mod_crp_anca      <- glm(group01 ~ CRP_Plasma_mg_L + zANCA,
                         data = train60, family = binomial())

# Convergence warnings
for (nm in c("mod_crp","mod_anca","mod_prot",
             "mod_prot_crp","mod_prot_anca","mod_prot_crp_anca","mod_crp_anca")) {
  m <- get(nm)
  if (!isTRUE(m$converged)) warning(nm, " did not converge — check for separation.")
}

# ============================================================
# SECTION 2: PREDICT ON TEST40
# ============================================================
test40$prob_prot          <- predict(mod_prot,          newdata = test40, type = "response")
test40$prob_crp           <- predict(mod_crp,           newdata = test40, type = "response")
test40$prob_anca          <- predict(mod_anca,          newdata = test40, type = "response")
test40$prob_prot_crp      <- predict(mod_prot_crp,      newdata = test40, type = "response")
test40$prob_prot_anca     <- predict(mod_prot_anca,     newdata = test40, type = "response")
test40$prob_prot_crp_anca <- predict(mod_prot_crp_anca, newdata = test40, type = "response")
test40$prob_crp_anca      <- predict(mod_crp_anca,      newdata = test40, type = "response")

# ============================================================
# SECTION 3: HELPERS
# ============================================================
binom_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int[1:2]
}

# Full metrics function (AUC + threshold metrics)
compute_metrics <- function(truth, prob, model_name, threshold = 0.5) {
  truth <- as.integer(truth)
  prob  <- as.numeric(prob)
  complete <- !is.na(truth) & !is.na(prob)
  truth <- truth[complete]; prob <- prob[complete]
  n      <- length(truth)
  events <- sum(truth == 1)
  non_ev <- sum(truth == 0)
  
  # AUC
  if (events > 0 && non_ev > 0) {
    roc_obj <- pROC::roc(truth, prob, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  } else {
    roc_obj <- NULL; auc_val <- NA_real_; auc_ci <- rep(NA_real_, 3)
  }
  
  # Threshold metrics
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
      model   = model_name, n = n, events = events,
      auc     = auc_val,  auc_lo = auc_ci[1], auc_hi = auc_ci[3],
      sens    = sens,     sens_lo = sens_ci[1], sens_hi = sens_ci[2],
      spec    = spec,     spec_lo = spec_ci[1], spec_hi = spec_ci[2],
      ppv     = ppv,      ppv_lo  = ppv_ci[1],  ppv_hi  = ppv_ci[2],
      npv     = npv,      npv_lo  = npv_ci[1],  npv_hi  = npv_ci[2],
      # --- added ---
      TP = TP, TN = TN, FP = FP, FN = FN,
      stringsAsFactors = FALSE
    )
  )
}
# ============================================================
# SECTION 4: AUC COMPARISON TABLE (test40, all models)
# ============================================================
res_prot          <- compute_metrics(test40$group01, test40$prob_prot,          "7-protein panel")
res_crp           <- compute_metrics(test40$group01, test40$prob_crp,           "CRP alone")
res_anca          <- compute_metrics(test40$group01, test40$prob_anca,          "ANCA alone")
res_crp_anca      <- compute_metrics(test40$group01, test40$prob_crp_anca,      "CRP + ANCA")
res_prot_crp      <- compute_metrics(test40$group01, test40$prob_prot_crp,      "7-protein + CRP")
res_prot_anca     <- compute_metrics(test40$group01, test40$prob_prot_anca,     "7-protein + ANCA")
res_prot_crp_anca <- compute_metrics(test40$group01, test40$prob_prot_crp_anca, "7-protein + CRP + ANCA")

metrics_all <- do.call(rbind, lapply(
  list(res_prot, res_crp, res_anca, res_crp_anca,
       res_prot_crp, res_prot_anca, res_prot_crp_anca),
  function(x) x$table
))
rownames(metrics_all) <- NULL

# Round
num_cols <- sapply(metrics_all, is.numeric)
metrics_all[, num_cols] <- lapply(metrics_all[, num_cols], round, 3)

cat("\n=== Model AUC comparison (test40) ===\n")
print(metrics_all[, c("model","n","events","auc","auc_lo","auc_hi","sens","spec","ppv","npv")])

write.csv(metrics_all,
          file.path(outdir, "TQL60_40_model_comparison_metrics.csv"),
          row.names = FALSE)

# ============================================================
# SECTION 5: DELONG TESTS (7-protein panel as reference)
# ============================================================


run_delong <- function(r1, r2, lab, boot_n = 2000) {
  if (is.null(r1) || is.null(r2)) {
    cat(sprintf("  %s: SKIPPED\n", lab)); return(NULL)
  }
  dl    <- pROC::roc.test(r1, r2, method = "delong")
  delta <- as.numeric(pROC::auc(r1)) - as.numeric(pROC::auc(r2))
  
  set.seed(7)
  boot_deltas <- replicate(boot_n, {
    c1_b  <- sample(r1$cases,    length(r1$cases),    replace = TRUE)
    ct1_b <- sample(r1$controls, length(r1$controls), replace = TRUE)
    c2_b  <- sample(r2$cases,    length(r2$cases),    replace = TRUE)
    ct2_b <- sample(r2$controls, length(r2$controls), replace = TRUE)
    as.numeric(pROC::auc(pROC::roc(
      response  = c(rep(1, length(c1_b)),  rep(0, length(ct1_b))),
      predictor = c(c1_b, ct1_b), quiet = TRUE))) -
      as.numeric(pROC::auc(pROC::roc(
        response  = c(rep(1, length(c2_b)),  rep(0, length(ct2_b))),
        predictor = c(c2_b, ct2_b), quiet = TRUE)))
  })
  
  ci_lo <- quantile(boot_deltas, 0.025)
  ci_hi <- quantile(boot_deltas, 0.975)
  
  cat(sprintf("  %-40s AUC1=%.3f  AUC2=%.3f  ΔAUC=%+.3f (95%% CI %+.3f – %+.3f)  p=%.4f\n",
              lab, as.numeric(pROC::auc(r1)), as.numeric(pROC::auc(r2)),
              delta, ci_lo, ci_hi, dl$p.value))
  
  list(D     = as.numeric(dl$statistic)[1],
       df    = as.numeric(dl$parameter)[1],
       p     = as.numeric(dl$p.value)[1],
       delta = as.numeric(delta)[1],
       ci_lo = as.numeric(ci_lo)[1],
       ci_hi = as.numeric(ci_hi)[1],
       auc1  = as.numeric(pROC::auc(r1))[1],
       auc2  = as.numeric(pROC::auc(r2))[1])
}

delong_to_row <- function(dl, comparison, type) {
  if (is.null(dl)) return(NULL)
  data.frame(
    comparison = comparison, type = type,
    auc1  = dl$auc1,  auc2  = dl$auc2,
    delta = dl$delta, ci_lo = dl$ci_lo, ci_hi = dl$ci_hi,
    D = dl$D, df = dl$df, p = dl$p,
    stringsAsFactors = FALSE
  )
}

cat("\n=== DeLong tests: 7-protein panel vs comparators (test40) ===\n")
dl_prot_vs_crp       <- run_delong(res_prot$roc_obj, res_crp$roc_obj,
                                   "7-prot vs CRP alone")
dl_prot_vs_anca      <- run_delong(res_prot$roc_obj, res_anca$roc_obj,
                                   "7-prot vs ANCA alone")
dl_prot_vs_crp_anca  <- run_delong(res_prot$roc_obj, res_crp_anca$roc_obj,
                                   "7-prot vs CRP+ANCA")
dl_prot_vs_protCRP   <- run_delong(res_prot$roc_obj, res_prot_crp$roc_obj,
                                   "7-prot vs 7-prot+CRP")
dl_prot_vs_protANCA  <- run_delong(res_prot$roc_obj, res_prot_anca$roc_obj,
                                   "7-prot vs 7-prot+ANCA")
dl_prot_vs_full      <- run_delong(res_prot$roc_obj, res_prot_crp_anca$roc_obj,
                                   "7-prot vs 7-prot+CRP+ANCA")

delong_all_6040 <- bind_rows(
  delong_to_row(dl_prot_vs_crp,      "7-prot vs CRP alone",      "incremental"),
  delong_to_row(dl_prot_vs_anca,     "7-prot vs ANCA alone",     "incremental"),
  delong_to_row(dl_prot_vs_crp_anca, "7-prot vs CRP+ANCA",       "incremental"),
  delong_to_row(dl_prot_vs_protCRP,  "7-prot vs 7-prot+CRP",     "incremental"),
  delong_to_row(dl_prot_vs_protANCA, "7-prot vs 7-prot+ANCA",    "incremental"),
  delong_to_row(dl_prot_vs_full,     "7-prot vs 7-prot+CRP+ANCA","incremental")
)

write_csv(delong_all_6040,
          file.path(outdir, "delong_all_comparisons_TQL60_40.csv"))

# ============================================================
# SECTION 6: ROC CURVES — all models (test40)
# ============================================================
# Colour palette (consistent with Prague script style)
cols_models <- c(
  "7-protein panel"        = "#C70500",
  "CRP alone"              = "#89CC48",
  "ANCA alone"             = "#66CCD1",
  "ANCA + CRP"   = "#A260D8",
  "7-protein + CRP"        = "#888888",
  "7-protein + ANCA"       = "#aaaaaa",
  "7-protein + CRP + ANCA" = "#BBC2D5"
)

# Build legend labels with AUC + 95% CI
make_label <- function(name, res) {
  t <- res$table
  sprintf("%s\nAUC %.3f (%.3f–%.3f)", name, t$auc, t$auc_lo, t$auc_hi)
}

leg_labels <- c(
  make_label("7-protein panel",        res_prot),
  make_label("CRP alone",              res_crp),
  make_label("ANCA alone",             res_anca),
  make_label("CRP + ANCA",             res_crp_anca),
  make_label("7-protein + CRP",        res_prot_crp),
  make_label("7-protein + ANCA",       res_prot_anca),
  make_label("7-protein + CRP + ANCA", res_prot_crp_anca)
)

svg_file <- file.path(outdir, "ROC_TQL60_40_model_comparison.svg")
svglite(svg_file, width = 6, height = 5.5)

plot(res_prot$roc_obj,
     col = cols_models["7-protein panel"], lwd = 2.2,
     xlab = "1 – Specificity", ylab = "Sensitivity",
     main = "TQL 60/40 — model comparison (test40)",
     cex.main = 0.95)
plot(res_crp$roc_obj,       col = cols_models["CRP alone"],              lwd = 1.5, lty = 2, add = TRUE)
plot(res_anca$roc_obj,      col = cols_models["ANCA alone"],             lwd = 1.5, lty = 2, add = TRUE)
plot(res_crp_anca$roc_obj,  col = cols_models["CRP + ANCA"],             lwd = 1.5, lty = 2, add = TRUE)
plot(res_prot_crp$roc_obj,      col = cols_models["7-protein + CRP"],        lwd = 2.2, add = TRUE)
plot(res_prot_anca$roc_obj,     col = cols_models["7-protein + ANCA"],       lwd = 2.2, add = TRUE)
plot(res_prot_crp_anca$roc_obj, col = cols_models["7-protein + CRP + ANCA"], lwd = 2.2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")

legend("bottomright",
       legend = leg_labels,
       col    = cols_models,
       lwd    = c(2.2, 1.5, 1.5, 1.5, 2.2, 2.2, 2.2),
       lty    = c(1, 2, 2, 2, 1, 1, 1),
       bty    = "n",
       cex    = 0.72,
       y.intersp = 1.5)

dev.off()
message("ROC plot saved: ", svg_file)

# ============================================================
# SECTION 7: MPO / PR3 SUBGROUP ANALYSIS
# ============================================================
test40$ANCA_group <- as.character(test40$AAV_group)

cat("\n=== MPO/PR3 subgroup counts (test40) ===\n")
print(table(test40$ANCA_group, test40$group01, useNA = "ifany"))

# Subgroup metrics function
compute_subgroup_metrics <- function(df, sg_value, score_col,
                                     model_label, threshold = 0.5) {
  sub <- df[df$ANCA_group == sg_value & !is.na(df$ANCA_group), , drop = FALSE]
  if (nrow(sub) == 0) {
    warning("Subgroup ", sg_value, " for model ", model_label, " is empty.")
    return(NULL)
  }
  res <- compute_metrics(sub$group01, sub[[score_col]], model_label, threshold)
  res$table$subgroup <- sg_value
  res$table
}

sg_models <- list(
  list(col = "prob_prot",          label = "7-protein panel"),
  list(col = "prob_crp",           label = "CRP alone"),
  list(col = "prob_anca",          label = "ANCA alone"),
  list(col = "prob_crp_anca",      label = "CRP + ANCA"),
  list(col = "prob_prot_crp",      label = "7-protein + CRP"),
  list(col = "prob_prot_anca",     label = "7-protein + ANCA"),
  list(col = "prob_prot_crp_anca", label = "7-protein + CRP + ANCA")
)

sg_rows <- list()
for (md in sg_models) {
  for (sg in c("MPO", "PR3")) {
    r <- compute_subgroup_metrics(test40, sg, md$col, md$label)
    if (!is.null(r)) sg_rows[[length(sg_rows) + 1]] <- r
  }
}

subgroup_table <- do.call(rbind, sg_rows)
rownames(subgroup_table) <- NULL
subgroup_table <- subgroup_table[, c("subgroup","model","n","events",
                                     "auc","auc_lo","auc_hi",
                                     "sens","sens_lo","sens_hi",
                                     "spec","spec_lo","spec_hi",
                                     "ppv","ppv_lo","ppv_hi",
                                     "npv","npv_lo","npv_hi")]

num_cols <- sapply(subgroup_table, is.numeric)
subgroup_table[, num_cols] <- lapply(subgroup_table[, num_cols], round, 3)

cat("\n=== MPO/PR3 subgroup metrics — all models ===\n")
print(subgroup_table[, c("subgroup","model","n","events","auc","auc_lo","auc_hi")])

write.csv(subgroup_table,
          file.path(outdir, "TQL60_40_MPO_PR3_subgroup_metrics.csv"),
          row.names = FALSE)

# ============================================================
# SECTION 8: FOREST-STYLE AUC PLOT — MPO vs PR3
# ============================================================
# Show AUC ± 95% CI for each model, split by ANCA subgroup

sg_plot_data <- subgroup_table %>%
  filter(model %in% c("7-protein panel",
                      "7-protein + CRP",
                      "7-protein + ANCA",
                      "7-protein + CRP + ANCA")) %>%
  mutate(
    model = factor(model, levels = c("7-protein panel",
                                     "7-protein + CRP",
                                     "7-protein + ANCA",
                                     "7-protein + CRP + ANCA")),
    subgroup = factor(subgroup, levels = c("MPO", "PR3"))
  )

cols_sg <- c(
  "7-protein panel"        = "mediumpurple2",
  "7-protein + CRP"        = "#d95f02",
  "7-protein + ANCA"       = "#e7298a",
  "7-protein + CRP + ANCA" = "#66a61e"
)

p_forest <- ggplot(sg_plot_data,
                   aes(x = auc, y = model, colour = model,
                       xmin = auc_lo, xmax = auc_hi)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.3f (%.3f–%.3f)", auc, auc_lo, auc_hi)),
            hjust = -0.08, size = 2.7, colour = "grey30") +
  facet_wrap(~ subgroup, ncol = 2) +
  scale_colour_manual(values = cols_sg) +
  scale_x_continuous(limits = c(0.55, 1.08),
                     breaks = seq(0.6, 1.0, 0.1)) +
  labs(title   = "TQL 60/40 test set — AUC by ANCA subgroup",
       subtitle = "Horizontal bars = 95% DeLong CI",
       x       = "AUC (95% CI)",
       y       = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

svg_sg <- file.path(outdir, "AUC_forest_MPO_PR3_TQL60_40.svg")
svglite(svg_sg, width = 8, height = 3.8)
print(p_forest)
dev.off()
message("Subgroup forest plot saved: ", svg_sg)

# ============================================================
# SECTION 9: ROC CURVES — 7-PROTEIN PANEL, MPO vs PR3
# (mirrors Prague_subgroup_tests.R style)
# ============================================================
build_sg_roc <- function(df, sg, score_col) {
  sub <- df[df$ANCA_group == sg & !is.na(df$ANCA_group), ]
  if (sum(sub$group01 == 1) == 0 || sum(sub$group01 == 0) == 0) return(NULL)
  pROC::roc(sub$group01, sub[[score_col]], quiet = TRUE)
}

roc_mpo_prot      <- build_sg_roc(test40, "MPO", "prob_prot")
roc_pr3_prot      <- build_sg_roc(test40, "PR3", "prob_prot")

roc_mpo_prot_crp  <- build_sg_roc(test40, "MPO", "prob_prot_crp")
roc_pr3_prot_crp  <- build_sg_roc(test40, "PR3", "prob_prot_crp")

roc_mpo_prot_anca <- build_sg_roc(test40, "MPO", "prob_prot_anca")
roc_pr3_prot_anca <- build_sg_roc(test40, "PR3", "prob_prot_anca")

# DeLong MPO vs PR3 within each model
cat("\n=== DeLong: MPO vs PR3 within each model (test40) ===\n")
mpo_vs_pr3_delong <- list()
for (nm in c("prot","prot_crp","prot_anca")) {
  roc_m <- get(paste0("roc_mpo_", nm))
  roc_p <- get(paste0("roc_pr3_", nm))
  if (!is.null(roc_m) && !is.null(roc_p)) {
    d <- tryCatch(pROC::roc.test(roc_m, roc_p, method = "delong"),
                  error = function(e) NULL)
    if (!is.null(d)) {
      cat(sprintf("%-28s MPO=%.3f | PR3=%.3f | p=%.4f\n",
                  nm, as.numeric(pROC::auc(roc_m)),
                  as.numeric(pROC::auc(roc_p)), d$p.value))
      mpo_vs_pr3_delong[[nm]] <- d
    }
  }
}

# Colour convention: MPO = #e7298a (pink), PR3 = #1b7837 (green)
# solid = 7-protein panel, dashed = augmented
if (!is.null(roc_mpo_prot) && !is.null(roc_pr3_prot)) {
  
  svg_sg_roc <- file.path(outdir, "ROC_MPO_PR3_7prot_augmented_TQL60_40.svg")
  svglite(svg_sg_roc, width = 5.5, height = 5.2)
  
  auc_fmt <- function(roc_obj) {
    a  <- as.numeric(pROC::auc(roc_obj))
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    sprintf("AUC %.3f (%.3f–%.3f)", a, ci[1], ci[3])
  }
  
  plot(roc_mpo_prot, col = "#e7298a", lwd = 2.2,
       xlab = "1 – Specificity", ylab = "Sensitivity",
       main = "TQL 60/40 — MPO vs PR3 subgroup (test40)",
       cex.main = 0.9)
  plot(roc_pr3_prot, col = "#1b7837", lwd = 2.2, add = TRUE)
  
  if (!is.null(roc_mpo_prot_crp))
    plot(roc_mpo_prot_crp, col = "#e7298a", lwd = 1.5, lty = 2, add = TRUE)
  if (!is.null(roc_pr3_prot_crp))
    plot(roc_pr3_prot_crp, col = "#1b7837", lwd = 1.5, lty = 2, add = TRUE)
  if (!is.null(roc_mpo_prot_anca))
    plot(roc_mpo_prot_anca, col = "#e7298a", lwd = 1.5, lty = 3, add = TRUE)
  if (!is.null(roc_pr3_prot_anca))
    plot(roc_pr3_prot_anca, col = "#1b7837", lwd = 1.5, lty = 3, add = TRUE)
  
  abline(a = 0, b = 1, lty = 3, col = "grey60")
  
  legend_entries <- character(0); legend_cols <- character(0)
  legend_lwd <- numeric(0); legend_lty <- numeric(0)
  
  if (!is.null(roc_mpo_prot)) {
    legend_entries <- c(legend_entries, paste0("MPO — 7-prot panel\n", auc_fmt(roc_mpo_prot)))
    legend_cols <- c(legend_cols, "#e7298a"); legend_lwd <- c(legend_lwd, 2.2); legend_lty <- c(legend_lty, 1)
  }
  if (!is.null(roc_pr3_prot)) {
    legend_entries <- c(legend_entries, paste0("PR3 — 7-prot panel\n", auc_fmt(roc_pr3_prot)))
    legend_cols <- c(legend_cols, "#1b7837"); legend_lwd <- c(legend_lwd, 2.2); legend_lty <- c(legend_lty, 1)
  }
  if (!is.null(roc_mpo_prot_crp)) {
    legend_entries <- c(legend_entries, paste0("MPO — +CRP\n", auc_fmt(roc_mpo_prot_crp)))
    legend_cols <- c(legend_cols, "#e7298a"); legend_lwd <- c(legend_lwd, 1.5); legend_lty <- c(legend_lty, 2)
  }
  if (!is.null(roc_pr3_prot_crp)) {
    legend_entries <- c(legend_entries, paste0("PR3 — +CRP\n", auc_fmt(roc_pr3_prot_crp)))
    legend_cols <- c(legend_cols, "#1b7837"); legend_lwd <- c(legend_lwd, 1.5); legend_lty <- c(legend_lty, 2)
  }
  if (!is.null(roc_mpo_prot_anca)) {
    legend_entries <- c(legend_entries, paste0("MPO — +ANCA\n", auc_fmt(roc_mpo_prot_anca)))
    legend_cols <- c(legend_cols, "#e7298a"); legend_lwd <- c(legend_lwd, 1.5); legend_lty <- c(legend_lty, 3)
  }
  if (!is.null(roc_pr3_prot_anca)) {
    legend_entries <- c(legend_entries, paste0("PR3 — +ANCA\n", auc_fmt(roc_pr3_prot_anca)))
    legend_cols <- c(legend_cols, "#1b7837"); legend_lwd <- c(legend_lwd, 1.5); legend_lty <- c(legend_lty, 3)
  }
  
  legend("bottomright",
         legend    = legend_entries,
         col       = legend_cols,
         lwd       = legend_lwd,
         lty       = legend_lty,
         bty       = "n",
         cex       = 0.65,
         y.intersp = 1.6)
  
  dev.off()
  message("MPO/PR3 ROC plot saved: ", svg_sg_roc)
}


# Colour convention: MPO = #e7298a (pink), PR3 = #1b7837 (green)
# solid = 7-protein panel, dashed = augmented
if (!is.null(roc_mpo_prot) && !is.null(roc_pr3_prot)) {
  
  svg_sg_roc <- file.path(outdir, "ROC_MPO_PR3_7prot_TQL60_40.svg")
  svglite(svg_sg_roc, width = 5.5, height = 5.2)
  
  auc_fmt <- function(roc_obj) {
    a  <- as.numeric(pROC::auc(roc_obj))
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    sprintf("AUC %.3f (%.3f–%.3f)", a, ci[1], ci[3])
  }
  
  plot(roc_mpo_prot, col = "#e7298a", lwd = 2.2,
       xlab = "1 – Specificity", ylab = "Sensitivity",
       main = "TQL 60/40 — MPO vs PR3 subgroup (test40)",
       cex.main = 0.9)
  plot(roc_pr3_prot, col = "#1b7837", lwd = 2.2, add = TRUE)
  
  abline(a = 0, b = 1, lty = 3, col = "grey60")
  
  legend_entries <- character(0); legend_cols <- character(0)
  legend_lwd <- numeric(0); legend_lty <- numeric(0)
  
  if (!is.null(roc_mpo_prot)) {
    legend_entries <- c(legend_entries, paste0("MPO — 7-prot panel\n", auc_fmt(roc_mpo_prot)))
    legend_cols <- c(legend_cols, "#e7298a"); legend_lwd <- c(legend_lwd, 2.2); legend_lty <- c(legend_lty, 1)
  }
  if (!is.null(roc_pr3_prot)) {
    legend_entries <- c(legend_entries, paste0("PR3 — 7-prot panel\n", auc_fmt(roc_pr3_prot)))
    legend_cols <- c(legend_cols, "#1b7837"); legend_lwd <- c(legend_lwd, 2.2); legend_lty <- c(legend_lty, 1)
  }
  
  legend("bottomright",
         legend    = legend_entries,
         col       = legend_cols,
         lwd       = legend_lwd,
         lty       = legend_lty,
         bty       = "n",
         cex       = 0.65,
         y.intersp = 1.6)
  
  dev.off()
  message("MPO/PR3 ROC plot saved: ", svg_sg_roc)
}


# ============================================================
# SECTION 11: SUMMARY TABLE — results text helper
# ============================================================
# Print a clean table suitable for copy-paste into results

summary_auc <- metrics_all[, c("model","n","events","auc","auc_lo","auc_hi")]
cat("\n", strrep("=", 65), "\n")
cat("RESULTS TABLE — TQL 60/40 test set AUC comparison\n")
cat(strrep("=", 65), "\n")
cat(sprintf("%-30s %4s %6s  %6s (%s–%s)\n",
            "Model", "n", "Events", "AUC", "95%CI lo", "95%CI hi"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(summary_auc))) {
  cat(sprintf("%-30s %4d %6d  %.3f  (%.3f–%.3f)\n",
              summary_auc$model[i],
              summary_auc$n[i],
              summary_auc$events[i],
              summary_auc$auc[i],
              summary_auc$auc_lo[i],
              summary_auc$auc_hi[i]))
}
cat(strrep("=", 65), "\n\n")


# ==================================================================
# CONFUSION MATRIX PLOT — TQL 60/40, 7-protein panel (test40)
# Append this block to the end of TQL6040_subgroup_tests.R
# Requires: test40 with column prob_prot and group01 (already in env)
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
      fill_class  = factor(ifelse(Predicted == 1, "Remission", "Active"),
                           levels = c("Active", "Remission")),
      # shape = TRUE class (circle = Active, triangle = Remission)
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
  
  # background shading by predicted class
  # bottom row (predicted Active) = yellow tint, top row (predicted Remission) = blue tint
  shading <- data.frame(
    xmin = c(-0.5, -0.5,  0.5,  0.5),
    xmax = c( 0.5,  0.5,  1.5,  1.5),
    ymin = c(-0.5,  0.5, -0.5,  0.5),
    ymax = c( 0.5,  1.5,  0.5,  1.5),
    fill = c("#FFFFC5", "#7373FF", "#7373FF", "#FFFFC5")  # TN, FN, FP, TP
  )
  
  p <- ggplot(plot_df, aes(x = x_jit, y = y_jit)) +
    
    geom_rect(data = shading,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = shading$fill, inherit.aes = FALSE, alpha = 0.5) +
    
    geom_hline(yintercept = 0.5, colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 0.5, colour = "grey40", linewidth = 0.5) +
    
    # fill = predicted class, shape = true class
    geom_point(aes(fill = fill_class, shape = shape_class),
               size = 3.8, stroke = 0.5, colour = "black", alpha = 0.88) +
    
    scale_fill_manual(
      name   = "Predicted class",
      values = c("Active" = "#ffd973", "Remission" = "#0000ab"),
      guide  = guide_legend(override.aes = list(shape = 21, size = 4))
    ) +
    scale_shape_manual(
      name   = "True class",
      values = c("Active" = 21, "Remission" = 24),
      guide  = guide_legend(override.aes = list(fill = "grey60", size = 4))
    ) +
    
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

# --- Run for 7-protein panel on test40 --------------------------------------
make_confusion_scatter(
  scores_df  = test40,
  pred_col   = "prob_prot",
  group_col  = "group01",
  outdir     = outdir,
  filename   = "confusion_scatter_7prot_test40.svg",
  threshold  = 0.5,
  title_str  = "7-protein panel — test set (40%)"
)


