# ============================================================
# OVERFITTING TESTS — All three models
# ============================================================
# Models tested:
#   1) Berlin 7-protein GLM  → validated on Prague
#      METHOD: Label-permutation on Prague AUC only
#              (can't refit across cohorts — different populations)
#
#   2) Train60 7-protein GLM → validated on Test40
#      METHOD: Full refit-permutation: shuffle who is train/test,
#              refit model each time, measure AUC gap
#
#   3) Train60 Concentration → validated on Test40
#      METHOD: Same full refit-permutation as Model 2
#
# HONEST CAVEAT:
#   - Berlin→Prague gap is NOT a clean overfitting test.
#     Different cohorts = generalizability, not just overfitting.
#     Label permutation tells you if Prague AUC > chance, not if
#     the model memorized Berlin.
#   - Train60→Test40 IS a proper overfitting test (same population,
#     random split). The refit-permutation is the gold standard here.
# ============================================================

library(pROC)
library(dplyr)
library(readr)
library(ggplot2)
library(svglite)

set.seed(7)

# ----------------------------------------------------------
# PATHS — adjust if needed
# ----------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
data_csv     <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
outdir       <- "../proteomics_ml/data/overfitting_output/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 2000  # increase to 5000 for publication; 2000 is fine for checks

# ----------------------------------------------------------
# LOAD FROZEN SCORES (already computed at model-save time)
# ----------------------------------------------------------
scores_prague      <- read_csv(file.path(artifact_dir, "scores_Prague_frozen.csv"))
scores_test40      <- read_csv(file.path(artifact_dir, "scores_test40_frozen.csv"))
scores_test40_con  <- read_csv(file.path(artifact_dir, "scores_test40_frozen_concentration.csv"))

# ----------------------------------------------------------
# LOAD MODELS
# ----------------------------------------------------------
model_berlin_7       <- readRDS(file.path(artifact_dir, "model_7panel_Berlin_glm.rds"))
model_train60_7      <- readRDS(file.path(artifact_dir, "model_7panel_Train60_glm.rds"))
model_train60_conc   <- readRDS(file.path(artifact_dir, "model_conc_Train60_glm.rds"))

# ----------------------------------------------------------
# LOAD RAW DATA (needed for refit-permutation on Train60 models)
# ----------------------------------------------------------
dfz <- read_csv(data_csv)

load("../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r")  # full_meta_train
load("../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r")   # full_meta_test

train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)

# For concentration model — needs VSN data
dfvsn <- readRDS(file = "../proteomics_ml/data/variance_dfvsn_final.rds")
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_trainn.r")  # data_train
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_test.r")    # data_test

vsn_train <- dfvsn %>% filter(Patient %in% full_meta_train$Patient)
vsn_test  <- dfvsn %>% filter(Patient %in% full_meta_test$Patient)


# Predictors
prots_7    <- c("AHSG_001", "CLEC3B_020", "COMP_023",
                "F9_027", "LRG1_042", "MCAM_047", "MRC1_073..")
formula_7  <- as.formula(paste("group01 ~", paste(prots_7, collapse = " + ")))

formula_conc <- group01 ~ AHSG.1 + CLEC3B.20 + COMP.23. + F9.27 +
  LRG1.42 + MCAM.47 + MRC1.73

# ----------------------------------------------------------
# HELPER: safe AUC (returns NA on failure, e.g. separation)
# ----------------------------------------------------------
safe_auc <- function(truth, pred) {
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(truth, pred, quiet = TRUE))),
    error = function(e) NA_real_
  )
}

# ----------------------------------------------------------
# HELPER: permutation plot + save
# ----------------------------------------------------------
plot_perm <- function(perm_vals, obs_val, p_val, title_str, fname) {
  df_plot <- data.frame(x = perm_vals)
  p <- ggplot(df_plot, aes(x = x)) +
    geom_histogram(bins = 50, fill = "grey70", color = "white") +
    geom_vline(xintercept = obs_val, color = "#c70000", linewidth = 1.2) +
    annotate("text",
             x     = obs_val,
             y     = Inf,
             vjust = 2,
             hjust = -0.1,
             label = sprintf("Observed = %.3f\np = %.4f", obs_val, p_val),
             color = "#c70000", size = 3.5) +
    labs(title = title_str, x = "Value under permutation", y = "Count") +
    theme_bw(base_size = 12)
  
  ggsave(file.path(outdir, fname), plot = p, width = 5, height = 4)
  invisible(p)
}

# ============================================================
# MODEL 1: Berlin → Prague
# METHOD: Label-permutation AUC on Prague
#
# What this tests: Is Prague AUC > chance?
# What this does NOT test: Whether Berlin overfitted to Berlin data.
#   For that you'd need Berlin internal cross-validation (separate analysis).
#   The cross-cohort gap is confounded by cohort differences (MPO/PR3 ratio,
#   treatment protocols, lab batch effects), so a simple gap test is misleading.
# ============================================================
cat("\n=== MODEL 1: Berlin GLM -> Prague (label permutation) ===\n")

truth_prague <- as.integer(scores_prague$group01)
pred_prague  <- as.numeric(scores_prague$pred_berlin_7)

auc_berlin_prague <- safe_auc(truth_prague, pred_prague)
cat(sprintf("  Prague AUC: %.3f\n", auc_berlin_prague))

# In-sample AUC on Berlin training set
berlin_train_data <- dfz %>% filter(cohort == "Berlin")
pred_berlin_train <- as.numeric(predict(model_berlin_7, 
                                        newdata = berlin_train_data, 
                                        type = "response"))
auc_berlin_train <- safe_auc(as.integer(berlin_train_data$group01), pred_berlin_train)
gap_berlin <- auc_berlin_train - auc_berlin_prague
cat(sprintf("  Berlin train AUC: %.3f\n", auc_berlin_train))
cat(sprintf("  Train-to-Prague gap: %.3f\n", gap_berlin))
cat("  NOTE: This gap reflects BOTH overfitting AND cohort differences.\n")
cat("        Do not interpret as pure overfitting.\n")
# this is nonsens test! 

# ============================================================
# MODEL 2: Train60 7-protein → Test40
# METHOD: Full refit-permutation (gold standard for same-population split)
#
# Logic: Randomly reassign who is in train vs. test, refit the GLM each time,
#        compute AUC gap. If your observed gap is not larger than the
#        distribution of permuted gaps, there's no strong overfitting signal.
#
# HONEST WARNING: With n=127 train and 7 predictors, you likely have
#   separation in many permuted folds too — NA AUCs will be dropped.
#   Report how many permutations were usable.
# ============================================================
cat("\n=== MODEL 2: Train60 7-protein -> Test40 (refit permutation) ===\n")

full60_40   <- rbind(train60, test40)
n_full      <- nrow(full60_40)
n_train_60  <- nrow(train60)

truth_test40_7 <- as.integer(scores_test40$group01)
pred_test40_7  <- as.numeric(scores_test40$pred_train60_7)

pred_train60_7_insample <- as.numeric(predict(model_train60_7,
                                              newdata = train60,
                                              type = "response"))
auc_train60_7_train <- safe_auc(as.integer(train60$group01), pred_train60_7_insample)
auc_train60_7_test  <- safe_auc(truth_test40_7, pred_test40_7)
gap_obs_60_7        <- auc_train60_7_train - auc_train60_7_test

cat(sprintf("  Train AUC: %.3f\n", auc_train60_7_train))
# Train AUC: 0.975
cat(sprintf("  Test  AUC: %.3f\n", auc_train60_7_test))
# Test  AUC: 0.974
cat(sprintf("  Observed gap: %.3f\n", gap_obs_60_7))
# Observed gap: 0.001

perm_gaps_60_7 <- replicate(N_PERM, {
  idx        <- sample(seq_len(n_full), size = n_train_60, replace = FALSE)
  p_train    <- full60_40[idx,  , drop = FALSE]
  p_test     <- full60_40[-idx, , drop = FALSE]
  
  mod <- tryCatch(
    glm(formula_7, data = p_train, family = binomial()),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      glm(formula_7, data = p_train, family = binomial())
    )
  )
  if (is.null(mod)) return(NA_real_)
  
  auc_tr <- safe_auc(as.integer(p_train$group01),
                     predict(mod, newdata = p_train, type = "response"))
  auc_te <- safe_auc(as.integer(p_test$group01),
                     predict(mod, newdata = p_test,  type = "response"))
  auc_tr - auc_te
})

n_usable_60_7  <- sum(!is.na(perm_gaps_60_7))
perm_gaps_60_7 <- perm_gaps_60_7[!is.na(perm_gaps_60_7)]
p_val_60_7     <- (sum(perm_gaps_60_7 >= gap_obs_60_7) + 1) / (length(perm_gaps_60_7) + 1)

cat(sprintf("  Permutation usable: %d / %d (dropped due to separation/convergence)\n",
            n_usable_60_7, N_PERM))
# Permutation usable: 2000 / 2000 (dropped due to separation/convergence)
cat(sprintf("  Permutation p (gap): %.4f\n", p_val_60_7))
# Permutation p (gap): 0.7446
cat(sprintf("  Verdict: %s\n",
            ifelse(p_val_60_7 < 0.05,
                   "Gap significantly larger than chance — overfitting signal.",
                   "Gap within chance variation — no strong overfitting.")))
# Verdict: Gap within chance variation — no strong overfitting

if (n_usable_60_7 < N_PERM * 0.5) {
  cat("  !! WARNING: >50% of permutations dropped due to separation.\n")
  cat("     This model has serious separation issues in the 7-protein normalized data.\n")
  cat("     Interpret permutation p with caution — it may be biased.\n")
}

plot_perm(perm_gaps_60_7, gap_obs_60_7, p_val_60_7,
          sprintf("Train60 7-prot -> Test40: AUC gap\n(refit permutation, n_usable=%d)", n_usable_60_7),
          "overfit_Train60_7prot_Test40_refit_perm.svg")

# ============================================================
# MODEL 3: Train60 Concentration → Test40
# METHOD: Full refit-permutation (same as Model 2)
#
# This model uses raw concentrations (VSN data) so separation is
# much less of a problem (~5% extreme predictions vs 29% for 7-protein).
# Permutation should be much cleaner here.
# ============================================================
cat("\n=== MODEL 3: Train60 Concentration -> Test40 (refit permutation) ===\n")

full_vsn    <- rbind(vsn_train, vsn_test)
n_full_vsn  <- nrow(full_vsn)
n_train_vsn <- nrow(vsn_train)

truth_test40_conc <- as.integer(scores_test40_con$group01)
pred_test40_conc  <- as.numeric(scores_test40_con$pred_train60_conc)

pred_vsn_train_insample <- as.numeric(predict(model_train60_conc,
                                              newdata = vsn_train,
                                              type = "response"))
auc_vsn_train <- safe_auc(as.integer(vsn_train$group01), pred_vsn_train_insample)
auc_vsn_test  <- safe_auc(truth_test40_conc, pred_test40_conc)
gap_obs_conc  <- auc_vsn_train - auc_vsn_test

cat(sprintf("  Train AUC: %.3f\n", auc_vsn_train))
# Train AUC: 0.947
cat(sprintf("  Test  AUC: %.3f\n", auc_vsn_test))
# Test  AUC: 0.908
cat(sprintf("  Observed gap: %.3f\n", gap_obs_conc))
# Observed gap: 0.039

perm_gaps_conc <- replicate(N_PERM, {
  idx        <- sample(seq_len(n_full_vsn), size = n_train_vsn, replace = FALSE)
  p_train    <- full_vsn[idx,  , drop = FALSE]
  p_test     <- full_vsn[-idx, , drop = FALSE]
  
  mod <- tryCatch(
    glm(formula_conc, data = p_train, family = binomial()),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      glm(formula_conc, data = p_train, family = binomial())
    )
  )
  if (is.null(mod)) return(NA_real_)
  
  auc_tr <- safe_auc(as.integer(p_train$group01),
                     predict(mod, newdata = p_train, type = "response"))
  auc_te <- safe_auc(as.integer(p_test$group01),
                     predict(mod, newdata = p_test,  type = "response"))
  auc_tr - auc_te
})

n_usable_conc  <- sum(!is.na(perm_gaps_conc))
perm_gaps_conc <- perm_gaps_conc[!is.na(perm_gaps_conc)]
p_val_conc     <- (sum(perm_gaps_conc >= gap_obs_conc) + 1) / (length(perm_gaps_conc) + 1)

cat(sprintf("  Permutation usable: %d / %d\n", n_usable_conc, N_PERM))
# Permutation usable: 2000 / 2000
cat(sprintf("  Permutation p (gap): %.4f\n", p_val_conc))
# Permutation p (gap): 0.3188
cat(sprintf("  Verdict: %s\n",
            ifelse(p_val_conc < 0.05,
                   "Gap significantly larger than chance — overfitting signal.",
                   "Gap within chance variation — no strong overfitting.")))
# Verdict: Gap within chance variation — no strong overfitting.

plot_perm(perm_gaps_conc, gap_obs_conc, p_val_conc,
          sprintf("Train60 Conc -> Test40: AUC gap\n(refit permutation, n_usable=%d)", n_usable_conc),
          "overfit_Train60_Conc_Test40_refit_perm.svg")

# ============================================================
# SUMMARY TABLE
# ============================================================
summary_overfit <- data.frame(
  model = c(
    "Train60 7-protein -> Test40",
    "Train60 Concentration -> Test40"
  ),
  test_type = c(
    "Refit permutation (train-test AUC gap)",
    "Refit permutation (train-test AUC gap)"
  ),
  train_AUC = c(
    round(auc_train60_7_train, 3),
    round(auc_vsn_train, 3)
  ),
  test_AUC = c(
    round(auc_train60_7_test, 3),
    round(auc_vsn_test, 3)
  ),
  observed_gap = c(
    round(gap_obs_60_7, 3),
    round(gap_obs_conc, 3)
  ),
  permutation_p = c(
    round(p_val_60_7, 4),
    round(p_val_conc, 4)
  ),
  n_perm_usable = c(
    n_usable_60_7,
    n_usable_conc
  ),
  interpretation = c(
    ifelse(p_val_60_7 < 0.05, "Overfitting signal", "No strong overfitting"),
    ifelse(p_val_conc < 0.05, "Overfitting signal", "No strong overfitting")
  ),
  stringsAsFactors = FALSE
)

print(summary_overfit)
write_csv(summary_overfit, file.path(outdir, "overfitting_summary_all_models.csv"))

