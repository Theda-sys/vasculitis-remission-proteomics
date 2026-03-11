### ============================================================
### GLOBAL GLM-7: FULL DIAGNOSTICS + SEPARATION CHECKS
### Adapted from recalibration_prm.R — same analysis, global data
### ============================================================

### 0) PATHS & CONFIG ----------------------------------------------------------
## !! UPDATE THESE THREE PATHS !!
train_r   <- "../Uwes_wishes/data/PRM_semi_pure_lassofeatures80iger_split.csv"
test_r    <- "../Uwes_wishes/data/PRM_semi_pure_lassofeatures20iger_split.csv"


outdir_base <- "../proteomics_ml/data/glm_model_diagnostics/finalx"
dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)

# Bootstrap / permutation settings (match PRM script exactly)
BOOT_BRIER_R <- 2000
BOOT_AUC_R   <- 2000
BOOT_CAL_R   <- 2000
B_OPTIMISM_B <- 2000
N_PERM       <- 2000

# Seeds — report these in Methods
SEED_MAIN  <- 7
SEED_PERM  <- 2025
SEED_OPTIM <- 42
set.seed(123)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

### 1) PACKAGES ----------------------------------------------------------------
library(pROC); library(boot); library(caret); library(dplyr); library(ggplot2)
library(tidyr); library(ResourceSelection); library(readr); library(tibble)
library(jsonlite); library(svglite); library(logistf); library(glmnet); library(scales)

writeLines(capture.output(sessionInfo()),
           file.path(outdir_base, "sessionInfo.txt"))

### 2) HELPERS (identical to PRM script) ---------------------------------------
safe_rds <- function(obj, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, file)
}

to01 <- function(x, positive = c("1","remission","rem")) {
  if (is.numeric(x) && all(unique(na.omit(x)) %in% c(0,1))) return(as.integer(x))
  x_chr <- tolower(as.character(x))
  out_int <- as.integer(ifelse(x_chr %in% tolower(positive), 1L, 0L))
  if (!all(out_int %in% c(0L,1L))) stop("Non-binary after to01 conversion.")
  out_int
}

bootstrap_ci <- function(df, stat_fun, R = 200, seed = 123) {
  set.seed(seed)
  b <- tryCatch(boot::boot(data = df, statistic = stat_fun, R = R),
                error = function(e) e)
  if (inherits(b, "error")) stop("bootstrap failed: ", b$message)
  ci <- as.numeric(quantile(b$t, c(0.025, 0.975)))
  list(boot_obj = b, ci = ci)
}

logit_safe <- function(p) {
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  log(p / (1 - p))
}

compute_auc_simple <- function(y, p) {
  roc_obj <- pROC::roc(response = y, predictor = p, quiet = TRUE)
  as.numeric(pROC::auc(roc_obj))
}

auc_delong <- function(y, p) {
  roc_obj <- tryCatch(pROC::roc(response = y, predictor = p, quiet = TRUE),
                      error = function(e) NULL)
  if (is.null(roc_obj)) return(list(auc = NA_real_, lo = NA_real_, hi = NA_real_, roc = NULL))
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
                 error = function(e) rep(NA_real_, 3))
  list(auc = auc_val, lo = ci[1], hi = ci[3], roc = roc_obj)
}

### 3) LOAD DATA ---------------------------------------------------------------
train.data <- read.csv(train_r, sep = ";", dec = ",",
                      stringsAsFactors = FALSE, check.names = FALSE)
test.data  <- read.csv(test_r,  sep = ";", dec = ",",
                      stringsAsFactors = FALSE, check.names = FALSE)

if ("Experimental_Group" %in% colnames(train.data)) {
  rownames(train.data) <- train.data$Experimental_Group
  train.data$Experimental_Group <- NULL
}
if ("Experimental_Group" %in% colnames(test.data)) {
  rownames(test.data) <- test.data$Experimental_Group
  test.data$Experimental_Group <- NULL
}

if (!("group" %in% colnames(train.data))) stop("prm_train missing 'group'")
if (!("group" %in% colnames(test.data)))  stop("prm_test missing 'group'")

write.csv(table(train.data$group), file = file.path(outdir_base, "train_group_table.csv"))
write.csv(table(test.data$group),  file = file.path(outdir_base, "test_group_table.csv"))

### 4) SELECT FEATURES ---------------------------------------------------------
# Try to load saved feature list; fall back to hard-coded 7 proteins
if (file.exists(final7_rds)) {
  features7_global <- readRDS(final7_rds)
  cat("Loaded features from RDS:", final7_rds, "\n")
} else {
  # These are the 7 Boruta-confirmed proteins used in the Berlin/PRM models
  features7_global <- c("MCAM_047","AHSG_001","LRG1_042",
                        "CLEC3B_020","COMP_023","MRC1_073","F9_027")
  cat("Using hard-coded 7-protein feature list.\n")
}

features7_global <- intersect(features7_global, colnames(train.data))
if (length(features7_global) != 7)
  stop("Expected 7 features in train.data, found: ", length(features7_global),
       "\nMissing: ", paste(setdiff(features7_global, colnames(train.data)), collapse = ", "))

train.data$group01 <- to01(train.data$group, positive = c("1","remission","rem"))
test.data$group01  <- to01(test.data$group,  positive = c("1","remission","rem"))
stopifnot(all(train.data$group01 %in% c(0,1)))
stopifnot(all(test.data$group01  %in% c(0,1)))

write.csv(table(train.data$group, train.data$group01),
          file = file.path(outdir_base, "train_group_mapping.csv"))
write.csv(table(test.data$group, test.data$group01),
          file = file.path(outdir_base, "test_group_mapping.csv"))

train7 <- train.data[, c(features7_global, "group01"), drop = FALSE]
test7  <- test.data[,  c(features7_global, "group01"), drop = FALSE]

cat(sprintf("Train n=%d  |  Test n=%d\n", nrow(train7), nrow(test7)))
cat(sprintf("Train events (remission=1): %d  |  Test events: %d\n",
            sum(train7$group01), sum(test7$group01)))

### 5) FIT GLM-7 ON TRAIN ONLY -------------------------------------------------
form7 <- as.formula(paste("group01 ~", paste(features7_global, collapse = " + ")))
glm7  <- glm(form7, data = train7, family = binomial())
safe_rds(glm7, file.path(outdir_base, "glm7_global_train_fit.rds"))
write.csv(coef(summary(glm7)), file = file.path(outdir_base, "glm7_coefficients.csv"))

if (!isTRUE(glm7$converged)) warning("GLM-7 did not converge — separation very likely.")

p7_test <- predict(glm7, newdata = test7, type = "response")
y7_test <- test7$group01
stopifnot(length(p7_test) == length(y7_test), all(!is.na(p7_test)))

### 6) EVALUATION FUNCTION (identical to PRM script) --------------------------
eval_binary_probs <- function(y01, p, n_bins = 5,
                              R_brier    = BOOT_BRIER_R,
                              R_auc_boot = BOOT_AUC_R,
                              R_cal      = BOOT_CAL_R,
                              seed       = SEED_MAIN) {
  stopifnot(all(y01 %in% c(0,1)), length(y01) == length(p))
  
  roc_obj       <- pROC::roc(response = y01, predictor = p, quiet = TRUE)
  auc_val       <- as.numeric(pROC::auc(roc_obj))
  auc_ci_delong <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  
  df_auc   <- data.frame(y = y01, p = p)
  auc_stat <- function(d, i) as.numeric(pROC::roc(d$y[i], d$p[i], quiet = TRUE)$auc)
  boot_auc_res <- bootstrap_ci(df_auc, auc_stat, R = R_auc_boot, seed = seed)
  
  brier_val <- mean((p - y01)^2)
  df_brier  <- data.frame(p = p, y = y01)
  brier_stat <- function(d, i) mean((d$p[i] - d$y[i])^2)
  boot_brier_res <- bootstrap_ci(df_brier, brier_stat, R = R_brier, seed = seed + 1)
  
  cal_df  <- tibble(pred = p, obs = y01) %>% mutate(bin = ntile(pred, n_bins))
  cal_sum <- cal_df %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(pred), obs_rate = mean(obs), n = n(), .groups = "drop") %>%
    arrange(mean_pred)
  
  boot_ci_bin <- function(obs_vec, R = R_cal, seed_bin = 1) {
    if (length(obs_vec) <= 1) return(c(NA_real_, NA_real_))
    set.seed(seed_bin)
    b <- boot::boot(data = data.frame(obs = obs_vec),
                    statistic = function(d, i) mean(d$obs[i]), R = R)
    as.numeric(quantile(b$t, c(0.025, 0.975)))
  }
  
  cis <- lapply(seq_len(nrow(cal_sum)), function(i)
    boot_ci_bin(cal_df$obs[cal_df$bin == cal_sum$bin[i]], R = R_cal, seed_bin = seed + i))
  if (length(cis)) {
    cis <- do.call(rbind, cis)
    cal_sum$ci_low  <- cis[,1]
    cal_sum$ci_high <- cis[,2]
  } else {
    cal_sum$ci_low  <- numeric(nrow(cal_sum))
    cal_sum$ci_high <- numeric(nrow(cal_sum))
  }
  
  # H-L test: suppress the "did not allow for requested bins" warning,
  # which occurs when predictions are near 0/1 (separation) or n is very small.
  # Record the warning text so it can be reported rather than silently dropped.
  hl_warn <- NULL
  hl_obj  <- tryCatch(
    withCallingHandlers(
      ResourceSelection::hoslem.test(x = y01, y = p, g = n_bins),
      warning = function(w) {
        hl_warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (!is.null(hl_warn)) {
    message("H-L warning (n=", length(y01), ", g=", n_bins, "): ", hl_warn)
    message("  -> H-L p-value unreliable; report as not estimable for this test set.")
    if (!is.null(hl_obj)) hl_obj$reliable <- FALSE else
      hl_obj <- list(reliable = FALSE, note = hl_warn)
  } else if (!is.null(hl_obj)) {
    hl_obj$reliable <- TRUE
  }
  
  # Calibration slope + intercept — point estimates
  cal_slope_intercept <- tryCatch({
    fit <- glm(y01 ~ logit_safe(p), family = binomial())
    coef(fit)
  }, error = function(e) c(NA_real_, NA_real_))
  
  # Bootstrap CI for calibration intercept + slope  [NEW]
  df_cal <- data.frame(y = y01, p = p)
  
  cal_boot_stat <- function(d, i) {
    tryCatch(
      coef(glm(y ~ logit_safe(p), data = d[i, ], family = binomial())),
      error   = function(e) c(NA_real_, NA_real_),
      warning = function(w) suppressWarnings(
        coef(glm(y ~ logit_safe(p), data = d[i, ], family = binomial()))
      )
    )
  }
  
  boot_cal_res <- tryCatch(
    boot::boot(data = df_cal, statistic = cal_boot_stat, R = R_cal),
    error = function(e) NULL
  )
  
  # Percentile CIs — graceful fallback to NA if bootstrap fails
  if (!is.null(boot_cal_res)) {
    ci_cal_intercept <- tryCatch(
      as.numeric(boot::boot.ci(boot_cal_res, index = 1, type = "perc")$percent[4:5]),
      error = function(e) c(NA_real_, NA_real_)
    )
    ci_cal_slope <- tryCatch(
      as.numeric(boot::boot.ci(boot_cal_res, index = 2, type = "perc")$percent[4:5]),
      error = function(e) c(NA_real_, NA_real_)
    )
  } else {
    ci_cal_intercept <- c(NA_real_, NA_real_)
    ci_cal_slope     <- c(NA_real_, NA_real_)
  }
  
  list(roc_obj              = roc_obj,
       auc                  = auc_val,
       auc_ci_delong        = auc_ci_delong,
       auc_ci_boot          = boot_auc_res$ci,
       brier                = brier_val,
       brier_ci             = boot_brier_res$ci,
       calibration          = cal_sum,
       hosmer_lemeshow      = hl_obj,
       cal_slope_intercept  = cal_slope_intercept,
       cal_intercept_ci     = ci_cal_intercept,   # NEW: c(lo, hi)
       cal_slope_ci         = ci_cal_slope,        # NEW: c(lo, hi)
       bootstrap_auc_obj    = boot_auc_res$boot_obj,
       bootstrap_brier_obj  = boot_brier_res$boot_obj,
       bootstrap_cal_obj    = boot_cal_res)         # NEW: full boot object
  
}

### 7) RUN DIAGNOSTICS ON TEST SET --------------------------------------------
res7 <- eval_binary_probs(y7_test, p7_test, n_bins = 5,
                          R_brier = BOOT_BRIER_R, R_auc_boot = BOOT_AUC_R,
                          R_cal = BOOT_CAL_R, seed = SEED_MAIN)

cat("\n--- PRM GLM-7 TEST SET ---\n")
cat(sprintf("AUC        = %.4f (DeLong 95%% CI: %.4f–%.4f)\n",
            res7$auc, res7$auc_ci_delong[1], res7$auc_ci_delong[3]))
cat(sprintf("AUC        = %.4f (Bootstrap 95%% CI: %.4f–%.4f)\n",
            res7$auc, res7$auc_ci_boot[1], res7$auc_ci_boot[2]))
cat(sprintf("Brier      = %.5f (Bootstrap 95%% CI: %.5f–%.5f)\n",
            res7$brier, res7$brier_ci[1], res7$brier_ci[2]))
cat(sprintf("Cal slope  = %.3f  (ideal = 1)\n",  res7$cal_slope_intercept[2]))
cat(sprintf("Cal intcpt = %.3f  (ideal = 0)\n",  res7$cal_slope_intercept[1]))
if (!is.null(res7$hosmer_lemeshow))
  cat(sprintf("H-L p      = %.4f\n", res7$hosmer_lemeshow$p.value))

write.csv(data.frame(
  model              = "PRM GLM-7",
  n_test             = length(y7_test),
  events_test        = sum(y7_test),
  AUC                = round(res7$auc, 4),
  AUC_delong_lo      = round(res7$auc_ci_delong[1], 4),
  AUC_delong_hi      = round(res7$auc_ci_delong[3], 4),
  AUC_boot_lo        = round(res7$auc_ci_boot[1], 4),
  AUC_boot_hi        = round(res7$auc_ci_boot[2], 4),
  Brier              = round(res7$brier, 5),
  Brier_lo           = round(res7$brier_ci[1], 5),
  Brier_hi           = round(res7$brier_ci[2], 5),
  Cal_intercept      = round(res7$cal_slope_intercept[1], 4),
  Cal_intercept_lo95 = round(res7$cal_intercept_ci[1], 4),   # NEW
  Cal_intercept_hi95 = round(res7$cal_intercept_ci[2], 4),   # NEW
  Cal_slope          = round(res7$cal_slope_intercept[2], 4),
  Cal_slope_lo95     = round(res7$cal_slope_ci[1], 4),        # NEW
  Cal_slope_hi95     = round(res7$cal_slope_ci[2], 4),        # NEW
  HL_p = if (is.null(res7$hosmer_lemeshow)) NA
  else if (!isTRUE(res7$hosmer_lemeshow$reliable)) "not estimable"
  else round(res7$hosmer_lemeshow$p.value, 4)
), file = file.path(outdir_base, "PRM_GLM7_test_metrics_summary.csv"),
row.names = FALSE)

# ── Quick formatted print for copy-paste into manuscript ─────────────────────
cat(sprintf("\n=== Calibration CI summary (bootstrap 95%%) ===\n"))
cat(sprintf("Intercept : %.3f  (95%% CI: %.3f \u2013 %.3f)\n",
            res7$cal_slope_intercept[1],
            res7$cal_intercept_ci[1],
            res7$cal_intercept_ci[2]))
cat(sprintf("Slope     : %.3f  (95%% CI: %.3f \u2013 %.3f)\n",
            res7$cal_slope_intercept[2],
            res7$cal_slope_ci[1],
            res7$cal_slope_ci[2]))

write.csv(res7$calibration,
          file.path(outdir_base, "PRM_GLM7_calibration_5bins_test.csv"),
          row.names = FALSE)
safe_rds(res7, file.path(outdir_base, "PRM_GLM7_test_eval.rds"))

# ── calibration stats for annotation ─────────────────────────────────────────
cal_intercept <- round(res7$cal_slope_intercept[1], 3)
cal_slope     <- round(res7$cal_slope_intercept[2], 3)
brier_val     <- round(res7$brier, 3)
hl_p_val <- if (is.null(res7$hosmer_lemeshow)) {
  "not computed"
} else if (!isTRUE(res7$hosmer_lemeshow$reliable)) {
  "not estimable (insufficient bin variation — separation)"
} else {
  formatC(res7$hosmer_lemeshow$p.value, digits = 3, format = "f")
}

stats_label <- paste0(
  "Intercept = ", cal_intercept, "\n",
  "Slope     = ", cal_slope,     "\n",
  "Brier     = ", brier_val,     "\n",
  "H-L p     = ", hl_p_val
)

# ── calibration plot ─────────────────────────────────────────────────────────
PRM_cal7 <- ggplot(res7$calibration, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.02) +
  annotate("text",
           x = 0.02, y = 0.98,
           label  = stats_label,
           hjust  = 0, vjust = 1,
           size   = 3.2,
           family = "mono") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  labs(title = "Calibration: PRM GLM-7 (20% test set, 5 bins)",
       x     = "Predicted probability",
       y     = "Observed event rate (binwise)")

svglite(file.path(outdir_base, "PRM_GLM7_test_calibration_plot.svg"), width = 4, height = 4)
print(PRM_cal7)
dev.off()

# ── calibration-in-the-large ─────────────────────────────────────────────────
mean_pred      <- mean(p7_test)
mean_obs       <- mean(y7_test)
cal_large_diff <- mean_pred - mean_obs
writeLines(c(
  paste0("Mean predicted probability: ", round(mean_pred, 3)),
  paste0("Observed event rate:        ", round(mean_obs,  3)),
  paste0("Absolute difference (pred - obs): ", round(cal_large_diff, 3))
), file.path(outdir_base, "PRM_GLM7_calibration_in_the_large.txt"))

### 8) LOO AUC SENSITIVITY ----------------------------------------------------
probs    <- p7_test
truth    <- y7_test
auc_orig <- compute_auc_simple(truth, probs)

aucs_loo <- sapply(seq_along(truth), function(i) {
  idx <- setdiff(seq_along(truth), i)
  compute_auc_simple(truth[idx], probs[idx])
})

write.csv(data.frame(loo_index = seq_along(truth), loo_auc = aucs_loo),
          file = file.path(outdir_base, "PRM_GLM7_test_LOO_AUCs.csv"), row.names = FALSE)
writeLines(sprintf("Original AUC = %.3f", auc_orig),
           file.path(outdir_base, "PRM_GLM7_test_LOO_AUC_original.txt"))

cat(sprintf("\nLOO AUC range: %.3f – %.3f  (original: %.3f)\n",
            min(aucs_loo), max(aucs_loo), auc_orig))

### 9) PERMUTATION TEST — 80/20 TEST SET --------------------------------------
set.seed(SEED_PERM)
perm_auc_test <- replicate(N_PERM, {
  yperm   <- sample(truth, replace = FALSE)
  roc_obj <- tryCatch(pROC::roc(yperm, probs, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
})
perm_auc_test <- perm_auc_test[!is.na(perm_auc_test)]
n_eff_perm    <- length(perm_auc_test)
p_perm_test   <- (sum(perm_auc_test >= auc_orig) + 1) / (n_eff_perm + 1)
p_perm_text   <- if (sum(perm_auc_test >= auc_orig) == 0)
  sprintf("p < 1/%d (p < %.6f)", n_eff_perm + 1, 1/(n_eff_perm + 1)) else
    sprintf("p = %.6g", p_perm_test)

cat(sprintf("Permutation test (test set): %s\n", p_perm_text))

safe_rds(list(perm_auc_test = perm_auc_test, p_value = p_perm_test, p_text = p_perm_text),
         file.path(outdir_base, "PRM_GLM7_test_label_permutation.rds"))
writeLines(c(paste0("Permutation test AUC (test set): ", p_perm_text),
             paste0("Effective permutations (non-NA): ", n_eff_perm)),
           file.path(outdir_base, "PRM_GLM7_test_label_permutation.txt"))

### 10) PERMUTATION TEST — FULL COHORT ----------------------------------------
full_dat   <- rbind(cbind(dataset = "train", train7), cbind(dataset = "test", test7))
full_probs <- as.numeric(predict(glm7, newdata = full_dat, type = "response"))
full_truth <- full_dat$group01
auc_full   <- compute_auc_simple(full_truth, full_probs)

set.seed(SEED_PERM + 1)
perm_auc_full <- replicate(N_PERM, {
  yperm   <- sample(full_truth, replace = FALSE)
  roc_obj <- tryCatch(pROC::roc(yperm, full_probs, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
})
perm_auc_full    <- perm_auc_full[!is.na(perm_auc_full)]
n_eff_full       <- length(perm_auc_full)
p_perm_full      <- (sum(perm_auc_full >= auc_full) + 1) / (n_eff_full + 1)
p_perm_full_text <- if (sum(perm_auc_full >= auc_full) == 0)
  sprintf("p < 1/%d (p < %.6f)", n_eff_full + 1, 1/(n_eff_full + 1)) else
    sprintf("p = %.6g", p_perm_full)

cat(sprintf("Full-cohort AUC = %.3f | Permutation test: %s\n", auc_full, p_perm_full_text))

safe_rds(list(perm_auc_full = perm_auc_full, p_value = p_perm_full, p_text = p_perm_full_text),
         file.path(outdir_base, "PRM_GLM7_full_label_permutation.rds"))
writeLines(c(sprintf("Full-cohort AUC = %.3f", auc_full),
             paste0("Permutation test AUC (full cohort): ", p_perm_full_text),
             paste0("Effective permutations (non-NA): ", n_eff_full)),
           file.path(outdir_base, "PRM_GLM7_full_label_permutation.txt"))

### 11) BOOTSTRAP OPTIMISM ON TRAIN -------------------------------------------
app_train_probs <- as.numeric(predict(glm7, newdata = train7, type = "response"))
app_train_auc   <- compute_auc_simple(train7$group01, app_train_probs)
set.seed(SEED_OPTIM)
opt_diffs <- numeric(B_OPTIMISM_B)

for (b in seq_len(B_OPTIMISM_B)) {
  idx        <- sample(seq_len(nrow(train7)), replace = TRUE)
  boot_train <- train7[idx, , drop = FALSE]
  boot_mod   <- tryCatch(glm(form7, data = boot_train, family = binomial()),
                         error = function(e) NULL)
  if (is.null(boot_mod)) { opt_diffs[b] <- NA; next }
  auc_boot_train    <- compute_auc_simple(boot_train$group01,
                                          as.numeric(predict(boot_mod, newdata = boot_train, type = "response")))
  auc_on_orig_train <- compute_auc_simple(train7$group01,
                                          as.numeric(predict(boot_mod, newdata = train7, type = "response")))
  opt_diffs[b] <- auc_boot_train - auc_on_orig_train
}

opt_est <- mean(opt_diffs, na.rm = TRUE)
cat(sprintf("\nApparent train AUC : %.4f\nOptimism estimate  : %.4f\nCorrected AUC      : %.4f\n",
            app_train_auc, opt_est, app_train_auc - opt_est))

writeLines(c(
  paste0("Apparent train AUC (glm7): ",     round(app_train_auc, 4)),
  paste0("Estimated optimism (mean over ",  B_OPTIMISM_B, " boot): ", round(opt_est, 4)),
  paste0("Optimism-corrected train AUC approx: ", round(app_train_auc - opt_est, 4))
), file.path(outdir_base, "bootstrap_optimism_summary.txt"))
safe_rds(list(opt_diffs = opt_diffs, app_train_auc = app_train_auc, opt_est = opt_est),
         file.path(outdir_base, "bootstrap_optimism_objects.rds"))

### 12) INTERCEPT-ONLY RECALIBRATION (SENSITIVITY) ----------------------------
recal_mod <- glm(y7_test ~ offset(logit_safe(p7_test)), family = binomial())
alpha_hat <- coef(recal_mod)[1]
p7_recal  <- plogis(logit_safe(p7_test) + alpha_hat)
safe_rds(list(alpha_hat = alpha_hat, recal_mod = recal_mod),
         file.path(outdir_base, "recalibration_intercept_only.rds"))

res7_recal <- eval_binary_probs(y7_test, p7_recal, n_bins = 5, seed = 777)
safe_rds(res7_recal, file.path(outdir_base, "PRM_GLM7_test_recalibrated.rds"))

cat(sprintf("\nRecalibrated AUC   = %.4f (unchanged by design)\n", res7_recal$auc))
cat(sprintf("Recalibrated Brier = %.5f  (original: %.5f)\n",
            res7_recal$brier, res7$brier))
cat(sprintf("Recalibrated cal intercept = %.3f  slope = %.3f\n",
            res7_recal$cal_slope_intercept[1], res7_recal$cal_slope_intercept[2]))

### 13) SEPARATION DIAGNOSTICS ------------------------------------------------
sep_outdir <- file.path(outdir_base, "separation_diagnostics")
dir.create(sep_outdir, recursive = TRUE, showWarnings = FALSE)

## 13a) Fitted value extremes on TRAIN ----------------------------------------
fitted_train <- fitted(glm7)
n_train      <- nrow(train7)
n_extreme_lo <- sum(fitted_train < 0.001)
n_extreme_hi <- sum(fitted_train > 0.999)
n_extreme    <- n_extreme_lo + n_extreme_hi
pct_extreme  <- round(100 * n_extreme / n_train, 1)

cat(sprintf("\nTRAIN (n=%d):\n", n_train))
cat(sprintf("  Fitted probs < 0.001 : %d\n", n_extreme_lo))
cat(sprintf("  Fitted probs > 0.999 : %d\n", n_extreme_hi))
cat(sprintf("  Total extreme        : %d / %d (%.1f%%)\n\n",
            n_extreme, n_train, pct_extreme))

sep_flag <- pct_extreme > 20
if (sep_flag) {
  cat("  ⚠  WARNING: >20% of training fitted values near 0/1 — quasi-complete separation.\n\n")
} else if (n_extreme > 0) {
  cat("  ℹ  Some extreme fitted values (<20%) — moderate concern.\n\n")
} else {
  cat("  ✓  No extreme fitted values — no separation detected.\n\n")
}

# Save separation diagnostics CSV (same format as calibration_all_models.R)
sep_csv <- data.frame(
  total_n         = n_train,
  pred_below_0001 = n_extreme_lo,
  pred_above_0999 = n_extreme_hi,
  n_extreme       = n_extreme,
  pct_extreme     = pct_extreme
)
write.csv(sep_csv,
          file.path(outdir_base, "separation_PRM_train.csv"),
          row.names = FALSE)

writeLines(c(
  sprintf("Train n = %d", n_train),
  sprintf("Fitted probs < 0.001: %d", n_extreme_lo),
  sprintf("Fitted probs > 0.999: %d", n_extreme_hi),
  sprintf("Total extreme: %d (%.1f%%)", n_extreme, pct_extreme),
  sprintf("Separation flag (>20%%): %s", ifelse(sep_flag, "YES — quasi-complete separation", "NO"))
), file.path(sep_outdir, "01_fitted_value_extremes.txt"))

## 13b) EPV on training data --------------------------------------------------
n_events_train <- sum(train7$group01 == 1)
n_pred         <- length(features7_global)
epv_train      <- n_events_train / n_pred
epv_flag       <- epv_train < 10

cat(sprintf("EPV CHECK:\n  Events=%d | Predictors=%d | EPV=%.2f\n",
            n_events_train, n_pred, epv_train))
if (epv_train < 5)  cat("  ⚠  EPV < 5: HIGH risk of overfitting and separation.\n\n") else
  if (epv_flag)     cat("  ⚠  EPV < 10: Moderate risk — report in supplement.\n\n") else
    cat("  ✓  EPV >= 10 — adequate.\n\n")

writeLines(c(
  sprintf("Events in train: %d", n_events_train),
  sprintf("Predictors: %d", n_pred),
  sprintf("EPV: %.2f", epv_train),
  sprintf("EPV flag: %s", ifelse(epv_flag, "BELOW guideline (<10)", "OK"))
), file.path(sep_outdir, "04_EPV_train.txt"))

## 13c) Coefficient magnitudes: GLM vs Ridge ----------------------------------
coef_glm <- coef(glm7)[-1]

X_train  <- as.matrix(train7[, features7_global])
y_train  <- train7$group01
set.seed(SEED_MAIN)
cv_ridge <- glmnet::cv.glmnet(X_train, y_train, alpha = 0, family = "binomial",
                              nfolds = min(10, nrow(train7)))
coef_ridge <- as.numeric(coef(cv_ridge, s = "lambda.1se"))[-1]
names(coef_ridge) <- features7_global

ratio_glm_ridge <- ifelse(abs(coef_ridge) < 1e-8, NA,
                          abs(coef_glm) / abs(coef_ridge))

coef_table <- data.frame(
  feature    = features7_global,
  coef_GLM   = round(coef_glm,   4),
  coef_Ridge = round(coef_ridge, 4),
  ratio_abs  = round(ratio_glm_ridge, 2),
  stringsAsFactors = FALSE
) %>% arrange(desc(ratio_abs))

cat("COEFFICIENT COMPARISON (GLM vs Ridge, sorted by |GLM/Ridge| ratio):\n")
print(coef_table, row.names = FALSE)
cat("\n")

max_ratio <- max(ratio_glm_ridge, na.rm = TRUE)
max_feat  <- features7_global[which.max(ratio_glm_ridge)]
coef_flag <- max_ratio > 3

cat(sprintf("  Largest ratio: %.2f-fold inflation for %s\n\n", max_ratio, max_feat))
if (coef_flag) cat("  ⚠  Ratios >3× — consistent with separation-driven inflation.\n\n") else
  cat("  ✓  Coefficients comparable — no strong inflation detected.\n\n")

write.csv(coef_table,
          file.path(sep_outdir, "02_coefficient_GLM_vs_Ridge.csv"), row.names = FALSE)

## 13d) Firth logistic regression as separation-robust reference AUC ----------
cat("FIRTH PENALISED LOGISTIC REGRESSION:\n")

form7_firth <- as.formula(paste("group01 ~", paste(features7_global, collapse = " + ")))
firth_mod   <- tryCatch(logistf::logistf(form7_firth, data = train7, pl = FALSE),
                        error = function(e) { message("logistf failed: ", e$message); NULL })

firth_flag <- NA
if (!is.null(firth_mod)) {
  X_test_mat <- model.matrix(form7_firth, data = test7)[, -1, drop = FALSE]
  lp_firth   <- as.numeric(cbind(1, X_test_mat) %*% coef(firth_mod))
  p_firth    <- plogis(lp_firth)
  
  roc_firth <- pROC::roc(y7_test, p_firth, quiet = TRUE)
  auc_firth <- as.numeric(pROC::auc(roc_firth))
  ci_firth  <- as.numeric(pROC::ci.auc(roc_firth, method = "delong"))
  
  roc_glm7  <- pROC::roc(y7_test, p7_test, quiet = TRUE)
  auc_glm7  <- as.numeric(pROC::auc(roc_glm7))
  ci_glm7   <- as.numeric(pROC::ci.auc(roc_glm7, method = "delong"))
  
  cat(sprintf("  GLM-7  AUC = %.3f (95%% CI: %.3f–%.3f)\n", auc_glm7, ci_glm7[1], ci_glm7[3]))
  cat(sprintf("  Firth  AUC = %.3f (95%% CI: %.3f–%.3f)\n", auc_firth, ci_firth[1], ci_firth[3]))
  cat(sprintf("  ΔAUC (GLM − Firth) = %.3f\n\n", auc_glm7 - auc_firth))
  
  firth_flag <- abs(auc_glm7 - auc_firth) > 0.02
  if (!firth_flag) cat("  ✓  AUC nearly identical — separation does not inflate test AUC.\n\n") else
    cat("  ⚠  ΔAUC >0.02 — separation may be inflating GLM test AUC.\n\n")
  
  cal_firth <- tryCatch(coef(glm(y7_test ~ logit_safe(p_firth), family = binomial())),
                        error = function(e) c(NA, NA))
  cal_glm_s <- tryCatch(coef(glm(y7_test ~ logit_safe(p7_test), family = binomial())),
                        error = function(e) c(NA, NA))
  
  writeLines(c(
    sprintf("GLM-7  AUC = %.3f (95%% CI: %.3f–%.3f)", auc_glm7, ci_glm7[1], ci_glm7[3]),
    sprintf("Firth  AUC = %.3f (95%% CI: %.3f–%.3f)", auc_firth, ci_firth[1], ci_firth[3]),
    sprintf("ΔAUC (GLM − Firth) = %.4f", auc_glm7 - auc_firth),
    sprintf("GLM   cal on test: intercept=%.3f, slope=%.3f", cal_glm_s[1], cal_glm_s[2]),
    sprintf("Firth cal on test: intercept=%.3f, slope=%.3f", cal_firth[1], cal_firth[2])
  ), file.path(sep_outdir, "03_Firth_vs_GLM_AUC_calibration.txt"))
} else {
  cat("  Firth model could not be fitted — skipping.\n\n")
}

## 13e) Calibration slope on test set -----------------------------------------
cat("CALIBRATION SLOPE ON 20% TEST SET:\n")

cal_test <- tryCatch(coef(glm(y7_test ~ logit_safe(p7_test), family = binomial())),
                     error = function(e) { cat("  Cal slope model failed\n"); c(NA, NA) })
hl_test  <- tryCatch(ResourceSelection::hoslem.test(x = y7_test, y = p7_test, g = 5),
                     error = function(e) NULL)

cat(sprintf("  Intercept : %.3f  (ideal = 0)\n", cal_test[1]))
cat(sprintf("  Slope     : %.3f  (ideal = 1)\n", cal_test[2]))
if (!is.null(hl_test)) cat(sprintf("  H-L p     : %.4f\n", hl_test$p.value))
cat("\n")

cal_flag <- !is.na(cal_test[2]) && cal_test[2] < 0.8
if (!is.na(cal_test[2])) {
  if (cal_test[2] < 0.6)    cat("  ⚠  Slope < 0.6: strong overconfidence — separation likely.\n\n") else
    if (cal_test[2] < 0.8)    cat("  ⚠  Slope 0.6–0.8: moderate overconfidence.\n\n") else
      cat("  ✓  Slope >= 0.8: acceptable calibration.\n\n")
}

writeLines(c(
  sprintf("Calibration intercept (test): %.4f", cal_test[1]),
  sprintf("Calibration slope (test):     %.4f", cal_test[2]),
  if (!is.null(hl_test)) sprintf("H-L p (test, g=5): %.4f", hl_test$p.value) else "H-L: not computed"
), file.path(sep_outdir, "05_calibration_slope_test.txt"))

## 13f) Overall verdict --------------------------------------------------------
n_flags <- sum(c(sep_flag, coef_flag, isTRUE(firth_flag), epv_flag, cal_flag), na.rm = TRUE)

verdict <- if (n_flags == 0) {
  "NO separation concern — AUC on the 20% test set appears robust."
} else if (n_flags <= 2) {
  paste0("MODERATE concern (", n_flags, "/5 flags). ",
         "AUC is likely real but may be modestly inflated. ",
         "Report Firth AUC alongside GLM AUC and note EPV limitation.")
} else {
  paste0("HIGH concern (", n_flags, "/5 flags). ",
         "Quasi-complete separation likely explains high AUC. ",
         "Primary inference should rest on external validation. ",
         "Report Firth AUC, calibration slope, and EPV in the rebuttal.")
}
cat("\n=== OVERALL VERDICT ===\n", verdict, "\n\n")
writeLines(c(
  sprintf("Flags triggered: %d / 5", n_flags),
  sprintf("sep_flag   : %s", sep_flag),
  sprintf("coef_flag  : %s", coef_flag),
  sprintf("firth_flag : %s", isTRUE(firth_flag)),
  sprintf("epv_flag   : %s", epv_flag),
  sprintf("cal_flag   : %s", cal_flag),
  "",
  verdict
), file.path(sep_outdir, "06_overall_verdict.txt"))

### 14) SAVE MAIN OBJECTS ------------------------------------------------------
saveRDS(list(glm7          = glm7,
             res7          = res7,
             res7_recal    = res7_recal,
             aucs_loo      = aucs_loo,
             perm_auc_test = perm_auc_test,
             perm_auc_full = perm_auc_full,
             opt_diffs     = opt_diffs),
        file = file.path(outdir_base, "all_saved_objects_full_checks.rds"))

### 15) PERMUTATION TEST ON TRAIN/TEST AUC GAP (overfitting check) ------------
auc_train_obs <- as.numeric(pROC::auc(pROC::roc(train7$group01,
                                                predict(glm7, newdata = train7, type = "response"),
                                                quiet = TRUE)))
auc_test_obs  <- as.numeric(pROC::auc(pROC::roc(y7_test, p7_test, quiet = TRUE)))
gap_obs       <- auc_train_obs - auc_test_obs

cat(sprintf("Train AUC : %.3f\n", auc_train_obs))
cat(sprintf("Test  AUC : %.3f\n", auc_test_obs))
cat(sprintf("Gap       : %.3f\n", gap_obs))

set.seed(SEED_PERM)
N_PERM_GAP <- 1000
full7       <- rbind(train7, test7)
n_full      <- nrow(full7)
n_tr        <- nrow(train7)

perm_gaps <- replicate(N_PERM_GAP, {
  idx        <- sample(seq_len(n_full), size = n_tr, replace = FALSE)
  perm_train <- full7[idx,  , drop = FALSE]
  perm_test  <- full7[-idx, , drop = FALSE]
  
  mod <- tryCatch(glm(form7, data = perm_train, family = binomial()),
                  error = function(e) NULL)
  if (is.null(mod)) return(NA_real_)
  
  auc_tr <- tryCatch(as.numeric(pROC::auc(pROC::roc(perm_train$group01,
                                                    predict(mod, newdata = perm_train, type = "response"),
                                                    quiet = TRUE))),
                     error = function(e) NA_real_)
  auc_te <- tryCatch(as.numeric(pROC::auc(pROC::roc(perm_test$group01,
                                                    predict(mod, newdata = perm_test,  type = "response"),
                                                    quiet = TRUE))),
                     error = function(e) NA_real_)
  auc_tr - auc_te
})
perm_gaps <- perm_gaps[!is.na(perm_gaps)]
p_gap <- (sum(perm_gaps >= gap_obs) + 1) / (length(perm_gaps) + 1)

cat(sprintf("Permutation p (gap): %.4f\n", p_gap))
cat(sprintf("  -> %s\n\n",
            ifelse(p_gap < 0.05,
                   "Gap is significantly larger than chance — evidence of overfitting.",
                   "Gap is within chance variation — no strong overfitting signal.")))

gap_df <- data.frame(gap = perm_gaps)
p_gap_plot <- ggplot(gap_df, aes(x = gap)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(xintercept = gap_obs, color = "#c70000", linewidth = 1.2) +
  annotate("text", x = gap_obs + 0.005, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("Observed gap = %.3f\np = %.4f", gap_obs, p_gap),
           color = "#c70000", size = 3.5) +
  labs(title = "Permutation test: Train/Test AUC gap (Global GLM-7)",
       subtitle = "Red line = observed gap. Distribution = gaps under random split.",
       x = "AUC gap (train - test)", y = "Count") +
  theme_bw(base_size = 12)

svglite(file.path(outdir_base, "overfitting_gap_permutation.svg"), width = 4, height = 4)
print(p_gap_plot)
dev.off()

writeLines(c(
  sprintf("Train AUC:        %.4f", auc_train_obs),
  sprintf("Test  AUC:        %.4f", auc_test_obs),
  sprintf("Observed gap:     %.4f", gap_obs),
  sprintf("Permutation n:    %d",   length(perm_gaps)),
  sprintf("Permutation p:    %.4f", p_gap),
  sprintf("Verdict: %s", ifelse(p_gap < 0.05, "Overfitting signal detected", "No strong overfitting"))
), file.path(outdir_base, "overfitting_gap_permutation.txt"))

### ============================================================
### 16) FULL MODEL ARCHIVE — GLM-7 for manuscript submission
### ============================================================
# This section saves everything a reader would need to reproduce
# predictions, inspect the model, or apply it to new data.
# Mirrors the archive structure of recalibration_global_glm21.R.

archive_dir <- file.path(outdir_base, "model_archive")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== ARCHIVING MODEL FOR MANUSCRIPT ===\n")

# ── 16a) Frozen model object --------------------------------------------------
safe_rds(glm7, file.path(archive_dir, "model_PRM_GLM7_frozen.rds"))

# ── 16b) Coefficient table (human-readable, with ORs) ------------------------
coef_summary <- as.data.frame(coef(summary(glm7)))
coef_summary$predictor <- rownames(coef_summary)
coef_summary <- coef_summary[, c("predictor", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
colnames(coef_summary) <- c("predictor", "beta", "se", "z", "p_value")
coef_summary$OR      <- round(exp(coef_summary$beta), 4)
coef_summary$OR_lo95 <- round(exp(coef_summary$beta - 1.96 * coef_summary$se), 4)
coef_summary$OR_hi95 <- round(exp(coef_summary$beta + 1.96 * coef_summary$se), 4)
coef_summary$beta    <- round(coef_summary$beta,  4)
coef_summary$se      <- round(coef_summary$se,    4)
coef_summary$z       <- round(coef_summary$z,     3)
coef_summary$p_value <- signif(coef_summary$p_value, 3)
write.csv(coef_summary, file.path(archive_dir, "model_Global_PRM_coefficients_OR.csv"),
          row.names = FALSE)

# ── 16c) Feature list ---------------------------------------------------------
write.csv(data.frame(feature  = features7_global,
                     position = seq_along(features7_global)),
          file.path(archive_dir, "model_PRM_GLM7_feature_list.csv"),
          row.names = FALSE)

# ── 16d) Performance metrics (all key stats in one table) --------------------
n_pred_7       <- length(features7_global)
n_events_7     <- sum(train7$group01)
epv_7          <- n_events_7 / n_pred_7

perf_table <- data.frame(
  metric  = c("n_train", "n_test", "events_train", "events_test",
              "n_predictors", "EPV",
              "AUC_test", "AUC_test_DeLong_lo", "AUC_test_DeLong_hi",
              "AUC_test_boot_lo", "AUC_test_boot_hi",
              "AUC_train_apparent", "optimism_estimate", "AUC_train_optimism_corrected",
              "AUC_full_cohort",
              "Brier_test", "Brier_test_lo", "Brier_test_hi",
              "Brier_baseline",
              "Cal_intercept_test", "Cal_slope_test",
              "HL_p_test",
              "Cal_mean_pred", "Cal_mean_obs", "Cal_in_the_large",
              "perm_p_test", "perm_p_full",
              "sep_pct_extreme_train", "sep_flag",
              "coef_flag", "epv_flag", "cal_flag", "firth_flag",
              "n_flags_total",
              "seed_main", "seed_perm", "seed_optim",
              "timestamp"),
  value   = c(nrow(train7), nrow(test7),
              n_events_7, sum(test7$group01),
              n_pred_7, round(epv_7, 2),
              round(res7$auc, 4),
              round(res7$auc_ci_delong[1], 4), round(res7$auc_ci_delong[3], 4),
              round(res7$auc_ci_boot[1], 4),   round(res7$auc_ci_boot[2], 4),
              round(auc_train_obs, 4), round(opt_est, 4),
              round(auc_train_obs - opt_est, 4),
              round(auc_full, 4),
              round(res7$brier, 5),
              round(res7$brier_ci[1], 5), round(res7$brier_ci[2], 5),
              round(mean(y7_test) * (1 - mean(y7_test)), 5),
              round(res7$cal_slope_intercept[1], 4),
              round(res7$cal_slope_intercept[2], 4),
              if (!is.null(res7$hosmer_lemeshow))
                formatC(res7$hosmer_lemeshow$p.value, digits = 4, format = "f") else NA,
              round(mean(p7_test), 4), round(mean(y7_test), 4),
              round(mean(p7_test) - mean(y7_test), 4),
              p_perm_text, p_perm_full_text,
              pct_extreme, as.character(sep_flag),
              as.character(coef_flag), as.character(epv_flag),
              as.character(cal_flag),  as.character(isTRUE(firth_flag)),
              n_flags,
              SEED_MAIN, SEED_PERM, SEED_OPTIM,
              timestamp),
  stringsAsFactors = FALSE
)
write.csv(perf_table,
          file.path(archive_dir, "model_PRM_GLM7_performance_archive.csv"),
          row.names = FALSE)

# ── 16e) Frozen test-set predictions (for reproduction) ----------------------
pred_archive <- data.frame(
  obs_index              = seq_along(y7_test),
  truth                  = y7_test,
  pred_prob              = round(p7_test, 6),
  logit_pred             = round(logit_safe(p7_test), 6),
  pred_prob_recalibrated = round(p7_recal, 6)
)
write.csv(pred_archive,
          file.path(archive_dir, "model_PRM_GLM7_test_predictions_frozen.csv"),
          row.names = FALSE)

# ── 16f) ROC curve points (for figure reproduction) --------------------------
roc_df <- data.frame(
  specificity = rev(res7$roc_obj$specificities),
  sensitivity = rev(res7$roc_obj$sensitivities),
  threshold   = rev(res7$roc_obj$thresholds)
)
write.csv(roc_df,
          file.path(archive_dir, "model_PRM_GLM7_ROC_test.csv"),
          row.names = FALSE)

# ── 16g) Calibration bins (for figure reproduction) --------------------------
write.csv(res7$calibration,
          file.path(archive_dir, "model_PRM_GLM7_calibration_bins_test.csv"),
          row.names = FALSE)

# ── 16h) Separation diagnostics summary (matches PRM format) -----------------
write.csv(data.frame(
  model           = "PRM GLM-7",
  training_set    = "Lasso80train (80% split)",
  total_n_train   = n_train,
  pred_below_0001 = n_extreme_lo,
  pred_above_0999 = n_extreme_hi,
  n_extreme       = n_extreme,
  pct_extreme     = pct_extreme,
  sep_flag        = sep_flag,
  n_events_train  = n_events_7,
  n_predictors    = n_pred_7,
  EPV             = round(epv_7, 2)
), file.path(archive_dir, "model_PRM_GLM7_separation_EPV_summary.csv"),
row.names = FALSE)

# ── 16i) Complete RDS bundle (everything in one file) ------------------------
safe_rds(
  list(
    # Model
    model             = glm7,
    features          = features7_global,
    formula           = form7,
    # Performance
    res_test          = res7,
    res_test_recal    = res7_recal,
    perf_table        = perf_table,
    coef_summary      = coef_summary,
    # Predictions
    y_test            = y7_test,
    p_test            = p7_test,
    p_test_recal      = p7_recal,
    # Robustness
    aucs_loo          = aucs_loo,
    perm_auc_test     = perm_auc_test,
    perm_auc_full     = perm_auc_full,
    opt_diffs         = opt_diffs,
    # Separation
    sep_csv           = sep_csv,
    coef_table_ridge  = coef_table,
    cv_ridge          = cv_ridge,
    # Recalibration
    alpha_hat         = alpha_hat,
    # Config
    seeds             = list(main = SEED_MAIN, perm = SEED_PERM, optim = SEED_OPTIM),
    timestamp         = timestamp,
    session_info      = sessionInfo()
  ),
  file.path(archive_dir, "model_PRM_GLM7_FULL_ARCHIVE.rds")
)

