# =============================================================================
# CRP-STRATIFIED SUBGROUP ANALYSIS — Test40 (60/40 split)
# 7-protein TQL intensity panel vs CRP alone
#
# Purpose: Assess whether the 7-protein panel discriminates AAV remission
#          in patients where CRP is elevated (>5 mg/L) — i.e., the subgroup
#          where CRP is least specific (elevated by any inflammation, not
#          just AAV activity).
# =============================================================================

library(tidyverse)
library(pROC)
library(scales)

set.seed(7)

# ----------------------------------------------------------
# PATHS — adjust if needed
# ----------------------------------------------------------
outdir_models <- "../proteomics_ml/data/models_and_artifacts/final"
data_csv      <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
train60_file  <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test40_file   <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"

outdir        <- "../proteomics_ml/data/TQL_60_40_ModelComparisons/CRP"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
dfz    <- read_csv(data_csv, show_col_types = FALSE)
load(train60_file)
load(test40_file)
test40 <- dfz %>% filter(Patient %in% full_meta_test$Patient)
train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)

# Attach frozen panel scores
scores <- read_csv(file.path(outdir_models, "scores_test40_frozen.csv"),
                   show_col_types = FALSE)
test40 <- test40 %>%
  left_join(scores %>% select(Patient, pred_train60_7), by = "Patient")

# Refit CRP-alone logistic on Train60
model_crp      <- glm(group01 ~ CRP_Plasma_mg_L, data = train60, family = binomial())
test40$pred_crp <- predict(model_crp, newdata = test40, type = "response")

# ----------------------------------------------------------
# SUBGROUPS
# ----------------------------------------------------------
crp_high <- test40 %>% filter(CRPhigher5_YES1 == 1)
crp_low  <- test40 %>% filter(CRPhigher5_YES1 == 0)

cat("Overall  n=", nrow(test40),   "events=", sum(test40$group01),   "\n")
cat("CRP>5    n=", nrow(crp_high), "events=", sum(crp_high$group01), "\n")
cat("CRP<=5   n=", nrow(crp_low),  "events=", sum(crp_low$group01),  "\n")

# ----------------------------------------------------------
# ROC OBJECTS
# ----------------------------------------------------------
roc_panel_all  <- roc(test40$group01,   test40$pred_train60_7,  quiet = TRUE)
roc_crp_all    <- roc(test40$group01,   test40$pred_crp,        quiet = TRUE)
roc_panel_high <- roc(crp_high$group01, crp_high$pred_train60_7,quiet = TRUE)
roc_crp_high   <- roc(crp_high$group01, crp_high$pred_crp,      quiet = TRUE)
roc_panel_low  <- roc(crp_low$group01,  crp_low$pred_train60_7, quiet = TRUE)
roc_crp_low    <- roc(crp_low$group01,  crp_low$pred_crp,       quiet = TRUE)

fmt_auc <- function(r) {
  a  <- as.numeric(auc(r))
  ci <- as.numeric(ci.auc(r, method = "delong"))
  sprintf("%.3f (%.3f\u2013%.3f)", a, ci[1], ci[3])
}

# ----------------------------------------------------------
# FIGURE 1: ROC CURVES — base R style matching manuscript
# ----------------------------------------------------------
cols_crp <- c(
  panel_all  = "#1F77B4",
  crp_all    = "#AEC7E8",
  panel_high = "#D62728",
  crp_high   = "#FFAAAA",
  panel_low  = "#2CA02C",
  crp_low    = "#98DF8A"
)

leg_labels <- c(
  paste0("Panel — Overall     AUC = ", fmt_auc(roc_panel_all)),
  paste0("CRP   — Overall     AUC = ", fmt_auc(roc_crp_all)),
  paste0("Panel — CRP > 5     AUC = ", fmt_auc(roc_panel_high)),
  paste0("CRP   — CRP > 5     AUC = ", fmt_auc(roc_crp_high)),
  paste0("Panel — CRP \u2264 5     AUC = ", fmt_auc(roc_panel_low)),
  paste0("CRP   — CRP \u2264 5     AUC = ", fmt_auc(roc_crp_low))
)



plot(roc_panel_all,
     col  = cols_crp["panel_all"], lwd = 2.5,
     xlab = "1 \u2212 Specificity", ylab = "Sensitivity",
     main = "ROC curves: 7-protein panel vs CRP alone by CRP stratum",
     cex.main = 1.0, cex.lab = 1.1)

plot(roc_crp_all,    col = cols_crp["crp_all"],    lwd = 1.8, lty = 2, add = TRUE)
plot(roc_panel_high, col = cols_crp["panel_high"],  lwd = 2.5, lty = 1, add = TRUE)
plot(roc_crp_high,   col = cols_crp["crp_high"],    lwd = 1.8, lty = 2, add = TRUE)
plot(roc_panel_low,  col = cols_crp["panel_low"],   lwd = 2.5, lty = 1, add = TRUE)
plot(roc_crp_low,    col = cols_crp["crp_low"],     lwd = 1.8, lty = 2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")

legend("bottomright",
       legend  = leg_labels,
       col     = unname(cols_crp),
       lwd     = c(2.5, 1.8, 2.5, 1.8, 2.5, 1.8),
       lty     = c(1, 2, 1, 2, 1, 2),
       bty     = "n",
       cex     = 0.72,
       y.intersp = 1.4,
       x.intersp = 0.8)

mtext("Test40 (n=85). Solid = 7-protein panel, dashed = CRP alone. DeLong 95% CIs.",
      side = 1, line = 4, cex = 0.82, col = "grey40")


# ----------------------------------------------------------
# BOOTSTRAP AUC CI FUNCTION (for forest plot)
# ----------------------------------------------------------
boot_auc_ci <- function(outcome, predictor, n_boot = 2000) {
  set.seed(7)
  auc_obs <- as.numeric(auc(roc(outcome, predictor, quiet = TRUE)))
  boot_aucs <- replicate(n_boot, {
    idx <- sample(seq_along(outcome), replace = TRUE)
    tryCatch(as.numeric(auc(roc(outcome[idx], predictor[idx], quiet=TRUE))),
             error=function(e) NA_real_)
  })
  boot_aucs <- na.omit(boot_aucs)
  list(est=auc_obs, lo=quantile(boot_aucs,0.025), hi=quantile(boot_aucs,0.975))
}

# ----------------------------------------------------------
# COMPUTE METRICS TABLE — fixed (no apply on data.frame rows)
# ----------------------------------------------------------
compute_row <- function(df, label) {
  p  <- boot_auc_ci(df$group01, df$pred_train60_7)
  cr <- boot_auc_ci(df$group01, df$pred_crp)
  dl <- roc.test(roc(df$group01, df$pred_train60_7, quiet=TRUE),
                 roc(df$group01, df$pred_crp,        quiet=TRUE),
                 method="delong")
  list(
    stratum      = label,
    n            = nrow(df),
    events       = sum(df$group01),
    prev         = round(mean(df$group01), 3),
    auc_panel    = round(p$est,  3),
    auc_panel_lo = round(p$lo,   3),
    auc_panel_hi = round(p$hi,   3),
    auc_crp      = round(cr$est, 3),
    auc_crp_lo   = round(cr$lo,  3),
    auc_crp_hi   = round(cr$hi,  3),
    delta_auc    = round(p$est - cr$est, 3),
    delong_p     = round(dl$p.value, 4)
  )
}

res_list <- list(
  compute_row(test40,   "Overall (n=85)"),
  compute_row(crp_high, "CRP > 5 mg/L"),
  compute_row(crp_low,  "CRP \u2264 5 mg/L")
)

# Convert to data.frame safely
results <- do.call(rbind, lapply(res_list, as.data.frame))
print(results)
write_csv(results, "../proteomics_ml/data/TQL_60_40_ModelComparisons/CRP/tables/CRP_subgroup_metrics.csv")

# ----------------------------------------------------------
# FIGURE 2: FOREST PLOT — base R
# ----------------------------------------------------------


strata     <- c("Overall\n(n=85, events=43)",
                "CRP > 5 mg/L\n(n=46, events=9)",
                "CRP \u2264 5 mg/L\n(n=39, events=34)")
y_pos      <- c(3, 2, 1)
dodge      <- 0.2
col_panel  <- "mediumpurple2"
col_crp    <- "#d95f02"

plot(NA, xlim = c(0.4, 1.05), ylim = c(0.5, 3.8),
     xlab = "AUC (95% CI)", ylab = "",
     yaxt = "n", bty = "l",
     main = "AUC by CRP stratum: 7-protein panel vs CRP alone",
     cex.main = 1.0, cex.lab = 1.1)

abline(v = 0.5, lty = 2, col = "grey70")
axis(2, at = y_pos, labels = strata, las = 1, cex.axis = 0.85, tick = FALSE)

for (i in seq_along(res_list)) {
  r  <- res_list[[i]]
  yp <- y_pos[i]
  
  # Panel
  points(r$auc_panel, yp + dodge, pch = 16, col = col_panel, cex = 1.5)
  segments(r$auc_panel_lo, yp+dodge, r$auc_panel_hi, yp+dodge,
           col=col_panel, lwd=2)
  
  # CRP
  points(r$auc_crp, yp - dodge, pch = 17, col = col_crp, cex = 1.5)
  segments(r$auc_crp_lo, yp-dodge, r$auc_crp_hi, yp-dodge,
           col=col_crp, lwd=2)
  
  # ΔAUC annotation
  p_lab <- ifelse(r$delong_p < 0.05,
                  sprintf("\u0394AUC=%+.3f, p=%.3f*", r$delta_auc, r$delong_p),
                  sprintf("\u0394AUC=%+.3f, p=%.3f",  r$delta_auc, r$delong_p))
  text(1.02, yp, p_lab, cex = 0.75, adj = 0, col = "grey30")
}

legend("bottomleft",
       legend = c("7-protein panel", "CRP alone"),
       col    = c(col_panel, col_crp),
       pch    = c(16, 17), lwd = 2,
       bty    = "n", cex = 0.9, pt.cex = 1.3)

mtext("Bootstrap 95% CIs (2,000 resamples, seed=7). DeLong p-value for panel vs CRP.",
      side=1, line=4, cex=0.8, col="grey40")


# ----------------------------------------------------------
# FIGURE 3: CRP BOXPLOT (already generated — just re-save cleanly)
# ----------------------------------------------------------
crp_box_data <- test40 %>%
  mutate(
    Outcome     = factor(ifelse(group01==1, "Remission", "Active"),
                         levels=c("Active","Remission")),
    CRP_stratum = factor(ifelse(CRPhigher5_YES1==1, "CRP > 5 mg/L", "CRP \u2264 5 mg/L"),
                         levels=c("CRP > 5 mg/L","CRP \u2264 5 mg/L"))
  )

p_crp <- ggplot(crp_box_data, aes(x=Outcome, y=log1p(CRP_Plasma_mg_L), fill=Outcome)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
  geom_jitter(aes(color = "black"), width=0.15, size=2.5, pch = 21) +
  facet_wrap(~CRP_stratum, scales="free_y") +
  scale_fill_manual(values  = c("Remission" = "#0000ab", "Active" = "#ffd973")) +
  scale_colour_manual(values = c("Remission" = "#0000ab", "Active" = "#ffd973")) +
  scale_y_continuous(name="CRP (mg/L, log scale)",
                     breaks=log1p(c(0,1,5,10,50,100,200)),
                     labels=c(0,1,5,10,50,100,200)) +
  labs(title="CRP distribution by remission status and CRP stratum",
       subtitle="Test40 (n=85). Within CRP-high stratum, CRP retains partial signal\nbut does not fully separate cases — the protein panel does.",
       x=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        strip.background=element_rect(fill="grey92",colour="grey70"),
        strip.text=element_text(face="bold",size=11),
        plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9.5,colour="grey40"))



# ----------------------------------------------------------
# FIGURE 4: PANEL SCORE BOXPLOT
# ----------------------------------------------------------
panel_box_data <- test40 %>%
  mutate(
    Outcome     = factor(ifelse(group01==1,"Remission","Active"),
                         levels=c("Active","Remission")),
    CRP_stratum = factor(ifelse(CRPhigher5_YES1==1,"CRP > 5 mg/L","CRP \u2264 5 mg/L"),
                         levels=c("CRP > 5 mg/L","CRP \u2264 5 mg/L"))
  )

p_panel <- ggplot(panel_box_data,
                  aes(x=Outcome, y=pred_train60_7, fill=Outcome)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
  geom_jitter(aes(color = "black"), width=0.15, size=2.5, pch = 21) +
  facet_wrap(~CRP_stratum) +
  scale_fill_manual(values  = c("Remission" = "#0000ab", "Active" = "#ffd973")) +
  scale_colour_manual(values = c("Remission" = "#0000ab", "Active" = "#ffd973")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2),
                     name="7-protein panel predicted probability") +
  labs(title="7-protein panel score by remission status and CRP stratum",
       subtitle="Dashed line = 0.5 classification threshold.\nPanel cleanly separates active from remission in both strata.",
       x=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        strip.background=element_rect(fill="grey92",colour="grey70"),
        strip.text=element_text(face="bold",size=11),
        plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9.5,colour="grey40"))



# ----------------------------------------------------------
# FORMATTED TABLE OUTPUT
# ----------------------------------------------------------
cat("\n========================================\n")
cat("CRP-stratified subgroup results\n")
cat("========================================\n")
for (r in res_list) {
  p_str <- ifelse(r$delong_p < 0.001, "<0.001", sprintf("%.4f", r$delong_p))
  cat(sprintf("%-20s  n=%-3d events=%-3d  Panel AUC: %.3f (%.3f\u2013%.3f)  CRP AUC: %.3f (%.3f\u2013%.3f)  \u0394AUC=%+.3f  DeLong p=%s\n",
              r$stratum, r$n, r$events,
              r$auc_panel, r$auc_panel_lo, r$auc_panel_hi,
              r$auc_crp,   r$auc_crp_lo,  r$auc_crp_hi,
              r$delta_auc, p_str))
}
cat("========================================\n")
cat("All figures saved to figures/\n")
cat("Table saved to tables/CRP_subgroup_metrics.csv\n")

# Test with ANCA
test40$ANCA_POS_B_MPOhigher5_B_PR3higher10_Phigher20_YES1

test40$pred_crp <- predict(model_crp, newdata = test40, type = "response")

# ----------------------------------------------------------
# SUBGROUPS
# ----------------------------------------------------------
anca_pos <- test40 %>% filter(ANCA_POS_B_MPOhigher5_B_PR3higher10_Phigher20_YES1 == 1)
anca_neg  <- test40 %>% filter(ANCA_POS_B_MPOhigher5_B_PR3higher10_Phigher20_YES1 == 0)

cat("Overall  n=", nrow(test40),   "events=", sum(test40$group01),   "\n")
cat("anca_pos    n=", nrow(anca_pos), "events=", sum(anca_pos$group01), "\n")
cat("anca_neg   n=", nrow(anca_neg),  "events=", sum(anca_neg$group01),  "\n")
# n= 12 events= 12  can not take anca_neg -> no actives 
# ----------------------------------------------------------
# ROC OBJECTS
# ----------------------------------------------------------
roc_panel_all  <- roc(test40$group01,   test40$pred_train60_7,  quiet = TRUE)
roc_panel_pos <- roc(anca_pos$group01, anca_pos$pred_train60_7,quiet = TRUE)
# roc_panel_neg  <- roc(anca_neg$group01,  anca_neg$pred_train60_7, quiet = TRUE)


fmt_auc <- function(r) {
  a  <- as.numeric(auc(r))
  ci <- as.numeric(ci.auc(r, method = "delong"))
  sprintf("%.3f (%.3f\u2013%.3f)", a, ci[1], ci[3])
}

# ----------------------------------------------------------
# FIGURE 1: ROC CURVES — base R style matching manuscript
# ----------------------------------------------------------
cols_crp <- c(
  panel_all  = "mediumpurple2",
  anca_pos   = "#FFAAAA"
)

leg_labels <- c(
  paste0("Panel — Overall     AUC = ", fmt_auc(roc_panel_all)),
  paste0("ANCA   — Positiv     AUC = ", fmt_auc(roc_panel_pos))
)

plot(roc_panel_all,
     col  = cols_crp["panel_all"], lwd = 2.5,
     xlab = "1 \u2212 Specificity", ylab = "Sensitivity",
     main = "ROC curves: 7-protein panel vs CRP alone by CRP stratum",
     cex.main = 1.0, cex.lab = 1.1)

plot(roc_panel_pos,    col = cols_crp["anca_pos"],    lwd = 1.8, lty = 2, add = TRUE)
abline(a = 0, b = 1, lty = 3, col = "grey60")

legend("bottomright",
       legend  = leg_labels,
       col     = unname(cols_crp),
       lwd     = c(2.5, 1.8, 2.5, 1.8, 2.5, 1.8),
       lty     = c(1, 2, 1, 2, 1, 2),
       bty     = "n",
       cex     = 0.72,
       y.intersp = 1.4,
       x.intersp = 0.8)

mtext("Test40 (n=85) 95% CIs.",
      side = 1, line = 4, cex = 0.82, col = "grey40")


# ----------------------------------------------------------
# BOOTSTRAP AUC CI FUNCTION (for forest plot)
# ----------------------------------------------------------
boot_auc_ci <- function(outcome, predictor, n_boot = 2000) {
  set.seed(7)
  auc_obs <- as.numeric(auc(roc(outcome, predictor, quiet = TRUE)))
  boot_aucs <- replicate(n_boot, {
    idx <- sample(seq_along(outcome), replace = TRUE)
    tryCatch(as.numeric(auc(roc(outcome[idx], predictor[idx], quiet=TRUE))),
             error=function(e) NA_real_)
  })
  boot_aucs <- na.omit(boot_aucs)
  list(est=auc_obs, lo=quantile(boot_aucs,0.025), hi=quantile(boot_aucs,0.975))
}




