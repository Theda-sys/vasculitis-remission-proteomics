# =============================================================================
# Relapse exploratory pipeline
# Primary test: univariable Cox with protein_score as continuous predictor
# Visualization: KM by median split and tertiles (illustration only)
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(broom)

# --- 0. Paths ----------------------------------------------------------------

data_csv      <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
model_dir     <- "../proteomics_ml/data/models_and_artifacts/final"
out_dir  <- "~/Documents/MDC/ml_proteomics/proteomics_ml/data/COX"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load data ------------------------------------------------------------

df <- read.csv(data_csv, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# --- 2. Protein score from frozen model --------------------------------------
# Always predict from the frozen Train60 model — never refit here.
# This guarantees the score is identical to the one reported in the main analysis.

proteins <- c("AHSG_001", "CLEC3B_020", "COMP_023", "F9_027",
              "LRG1_042", "MCAM_047", "MRC1_073..")

missing_prots <- setdiff(proteins, colnames(df))
if (length(missing_prots) > 0) stop("Missing protein columns: ", paste(missing_prots, collapse = ", "))

model_7panel    <- readRDS(file.path(model_dir, "model_7panel_Train60_glm.rds"))
df$protein_score <- predict(model_7panel, newdata = df, type = "response")
cat("protein_score computed for", sum(!is.na(df$protein_score)), "rows.\n")

# --- 3. Restrict to remission patients with 24-month endpoint ----------------

flare_24mo_col    <- "FOLLOWUP_flare_since_sampling_less24Mo_YES1"
flare_days_col    <- "FOLLOWUP_flare_since_sampling_days"
lost_days_col     <- "FOLLOWUP.lost_to_followup_days"

stopifnot(flare_24mo_col %in% colnames(df))

rem_df <- df %>%
  filter(disease_numeric == 1) %>%
  mutate(
    flare_24mo    = as.integer(.data[[flare_24mo_col]]),
    flare_days    = suppressWarnings(as.numeric(.data[[flare_days_col]])),
    lost_days     = if (lost_days_col %in% colnames(.)) suppressWarnings(as.numeric(.data[[lost_days_col]])) else NA_real_,
    time_to_event = case_when(
      !is.na(flare_days)                    ~ flare_days,
      is.na(flare_days) & !is.na(lost_days) ~ lost_days
    ),
    event = as.integer(!is.na(flare_days))
  )

# Analysis dataset: remission patients with protein_score AND 24mo endpoint
rem_df2 <- rem_df %>% filter(!is.na(protein_score), !is.na(flare_24mo))

cat("\nRemission rows total:                  ", nrow(rem_df),
    "\nWith protein_score + 24mo endpoint:    ", nrow(rem_df2),
    "\nFlares within 24 months:               ", sum(rem_df2$flare_24mo),
    "\nNo flare / censored within 24 months:  ", sum(rem_df2$flare_24mo == 0), "\n")

if (nrow(rem_df2) < 10 || sum(rem_df2$flare_24mo) < 3) {
  stop("Insufficient data: n = ", nrow(rem_df2), ", events = ", sum(rem_df2$flare_24mo))
}

# --- 4. Primary analysis: logistic regression + AUC -------------------------

fit_logistic <- glm(flare_24mo ~ protein_score, data = rem_df2, family = binomial)

logistic_tidy <- broom::tidy(fit_logistic, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\nLogistic regression (flare within 24mo ~ protein_score):\n")
print(logistic_tidy)
write_csv(logistic_tidy, file.path(out_dir, "relapse_24mo_logistic.csv"))

# AUC
roc_obj <- pROC::roc(rem_df2$flare_24mo, rem_df2$protein_score,
                     quiet = TRUE, ci = TRUE)
auc_ci  <- as.numeric(pROC::ci(roc_obj))

cat(sprintf("\nAUC: %.3f (95%% CI %.3f–%.3f)\n", auc_ci[2], auc_ci[1], auc_ci[3]))

auc_out <- tibble(
  AUC    = round(auc_ci[2], 3),
  CI_low = round(auc_ci[1], 3),
  CI_high= round(auc_ci[3], 3),
  n      = nrow(rem_df2),
  n_flare= sum(rem_df2$flare_24mo)
)
write_csv(auc_out, file.path(out_dir, "relapse_24mo_AUC.csv"))

# ROC plot
roc_plot <- pROC::ggroc(roc_obj, colour = "#2C7BB6", size = 1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.35, y = 0.1,
           label = sprintf("AUC = %.3f\n(95%% CI %.3f–%.3f)", auc_ci[2], auc_ci[1], auc_ci[3]),
           size = 4.5) +
  labs(title = "ROC: 7-protein score predicting flare within 24 months",
       x = "Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "ROC_protein_score_24mo.svg"), roc_plot, width = 5, height = 5)

# --- 5. Illustration: KM by median split (visualization only) ----------------
# NOTE: the binary logistic + AUC above is the primary analysis.
# KM is shown for illustration only — do not use for inference.

rem_km <- rem_df2 %>%
  filter(!is.na(time_to_event)) %>%
  mutate(score_group = ifelse(protein_score >= median(protein_score), "High", "Low"))

fit_km       <- survfit(Surv(time_to_event, event) ~ score_group, data = rem_km)
p_logrank    <- 1 - pchisq(survdiff(Surv(time_to_event, event) ~ score_group, data = rem_km)$chisq, df = 1)

cat(sprintf("\nKM median split: n = %d | log-rank p = %.3f\n", nrow(rem_km), p_logrank))

km_plot <- ggsurvplot(
  fit_km, data = rem_km,
  risk.table = TRUE, pval = TRUE, conf.int = FALSE,
  palette = c("#2C7BB6", "#D7191C"),
  legend.title = "Protein score",
  legend.labs  = c("High", "Low"),
  xlab = "Days since sampling",
  ylab = "Probability remaining flare-free",
  ggtheme = theme_bw(base_size = 14),
  risk.table.height = 0.22,
  tables.theme = theme_bw(base_size = 10)
)
ggsave(file.path(out_dir, "KM_protein_score_24mo_mediansplit.svg"), km_plot$plot, width = 8, height = 6)

# --- 6. Summary and suggested text -------------------------------------------

write_csv(
  tibble(
    n_remission_with_endpoint = nrow(rem_df2),
    n_flare_24mo              = sum(rem_df2$flare_24mo),
    n_no_flare_24mo           = sum(rem_df2$flare_24mo == 0),
    AUC                       = round(auc_ci[2], 3),
    AUC_CI_low                = round(auc_ci[1], 3),
    AUC_CI_high               = round(auc_ci[3], 3),
    OR_protein_score          = logistic_tidy$estimate[logistic_tidy$term == "protein_score"],
    OR_CI_low                 = logistic_tidy$conf.low[logistic_tidy$term == "protein_score"],
    OR_CI_high                = logistic_tidy$conf.high[logistic_tidy$term == "protein_score"],
    p_value                   = logistic_tidy$p.value[logistic_tidy$term == "protein_score"],
    logrank_p_km              = round(p_logrank, 3)
  ),
  file.path(out_dir, "relapse_24mo_summary.csv")
)

cat(sprintf(
  "\n--- Suggested text ---
Among %d patients sampled in remission with available 24-month follow-up,
%d (%d%%) experienced a flare within 24 months. The 7-protein predicted probability
(protein_score) showed an AUC of %.3f (95%% CI %.3f–%.3f) for discriminating
patients who flared within 24 months. In a logistic regression model,
OR = %.3f (95%% CI %.3f–%.3f), p = %s.
KM curves by median split are shown for illustration (log-rank p = %.3f).

CAVEAT: the protein_score was trained to classify active vs remission at sampling
and was not optimised for relapse prediction. These results are therefore exploratory.\n",
  nrow(rem_df2),
  sum(rem_df2$flare_24mo),
  round(100 * mean(rem_df2$flare_24mo)),
  auc_ci[2], auc_ci[1], auc_ci[3],
  logistic_tidy$estimate[logistic_tidy$term == "protein_score"],
  logistic_tidy$conf.low[logistic_tidy$term  == "protein_score"],
  logistic_tidy$conf.high[logistic_tidy$term == "protein_score"],
  format.pval(logistic_tidy$p.value[logistic_tidy$term == "protein_score"], digits = 3),
  p_logrank
))

cat("All outputs written to:", out_dir, "\n")

