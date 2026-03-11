# Benötigte Libraries
library(vsn)
library(preprocessCore)  # Für Quantile Normalisierung
library(ggplot2)
library(dplyr)
library(reshape2)

# Rohdaten einlesen (Spalten 2–8 = Peptide)
peptide_data <- as.matrix(concentration_calibartion_curve[, 2:8])

# 1. Log2-Transformation
# Log2-Transformation (+1, um Nullen zu vermeiden)
log2_data <- log2(peptide_data + 1)

# Visualisierung (Boxplot und Mean-SD)
par(mfrow = c(1, 2))
boxplot(log2_data, main = "Log2-Transformiert")
meanSdPlot(log2_data)

# 2. Quantile-Normalisierung (optional, falls gewünscht)
quantile_norm_data <- normalize.quantiles(peptide_data)
colnames(quantile_norm_data) <- colnames(peptide_data)
rownames(quantile_norm_data) <- rownames(peptide_data)

# Boxplot
boxplot(quantile_norm_data, main = "Quantile Normalisiert")

# 3. Loess-Normalisierung (lokale Regression)
loess_normalize <- function(data) {
  ref <- apply(data, 1, median)
  norm_data <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    fit <- loess(data[, i] ~ ref)
    norm_data[, i] <- data[, i] - predict(fit) + mean(data[, i], na.rm = TRUE)
  }
  colnames(norm_data) <- colnames(data)
  rownames(norm_data) <- rownames(data)
  return(norm_data)
}

# 4. Lineare Regression-Normalisierung
linear_reg_norm <- function(data) {
  ref <- apply(data, 1, median)
  norm_data <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    fit <- lm(data[, i] ~ ref)
    norm_data[, i] <- data[, i] - fitted(fit) + mean(data[, i], na.rm = TRUE)
  }
  colnames(norm_data) <- colnames(data)
  rownames(norm_data) <- rownames(data)
  return(norm_data)
}

loess_norm_data <- loess_normalize(peptide_data)
linreg_norm_data <- linear_reg_norm(peptide_data)

library(pROC)

# Korrigiere die Zielvariable (am besten ganz oben im Script, so ist die Codierung einheitlich):
concentration_calibartion_curve$group_bin <- ifelse(concentration_calibartion_curve$group == "active", 0, 1)

# Dann in der Schleife so referenzieren:
dat <- data.frame(
  Conc = norm_mat[, peptid],
  group = concentration_calibartion_curve$group_bin
)

dat$group <- as.numeric(as.character(dat$group))
table(dat$group[train_idx])                 # Nur 0/1 erlaubt!
stopifnot(all(dat$group[train_idx] %in% 0:1))

for (meth in names(norm_methods)) {
  norm_mat <- norm_methods[[meth]]
  for (peptid in colnames(norm_mat)) {
    dat <- data.frame(
      Conc = norm_mat[, peptid],
      group = concentration_calibartion_curve$group_bin  # Sicher als 0/1 Integer (nicht Text!)
    )
    train_idx <- concentration_calibartion_curve$cleaned %in% row.names(data_train)
    test_idx <- concentration_calibartion_curve$cleaned %in% row.names(data_test)
    # Werte-Prüfung
    stopifnot(all(dat$group[train_idx] %in% 0:1))
    model <- glm(group ~ Conc, data = dat[train_idx, ], family = binomial)
    pred <- predict(model, newdata = dat[test_idx, ], type = "response")
    roc_obj <- roc(dat$group[test_idx], pred)
    auc_results <- rbind(auc_results, data.frame(
      Peptid = peptid,
      Methode = meth,
      AUC = as.numeric(auc(roc_obj))
    ))
  }
}

auc_results[order(auc_results$AUC, decreasing = TRUE), ]


norm_methods <- list(
  Raw = peptide_data,
  Log2 = log2_data,
  Quantile = quantile_norm_data,
  Loess = loess_norm_data,
  LinReg = linreg_norm_data,
  VSN = vsn_raw
)

# Benutze für alle Methoden dieselben Metadaten
meta_cols <- c("Row.names", "group", "cleaned")

library(pROC)
library(caret)

panel_metrics <- data.frame(
  Methode = character(),
  AUC = numeric(),
  F1 = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  stringsAsFactors = FALSE
)

# Die Reihenfolge deiner Panel-Features:
panel_vars <- c("AHSG.1", "CLEC3B.20", "COMP.23.", "F9.27", "LRG1.42", "MCAM.47",  "MRC1.73")

for(meth in names(norm_methods)) {
  nm <- norm_methods[[meth]]
  # baue Dataframes wie bei vsn_raw
  panel_df <- data.frame(
    concentration_calibartion_curve[, meta_cols],
    nm[, panel_vars]
  )
  # CLEAN: Ensure group is coded 0/1 for GLM
  panel_df$group_bin <- ifelse(panel_df$group == "active", 0, 1)
  
  train_idx <- panel_df$cleaned %in% row.names(data_train)
  test_idx <- panel_df$cleaned %in% row.names(data_test)
  
  # GLM Panelmodell
  formula_panel <- as.formula(paste("group_bin ~", paste(panel_vars, collapse = " + ")))
  model <- glm(formula_panel, data = panel_df[train_idx, ], family = binomial)
  
  pred <- predict(model, newdata = panel_df[test_idx, ], type = "response")
  true <- panel_df$group_bin[test_idx]
  
  # Schwellenwert: 0.5
  predicted <- ifelse(pred > 0.5, 1, 0)
  predicted <- factor(predicted, levels = c(0,1))
  actual <- factor(true, levels = c(0,1))
  
  conf_matrix <- caret::confusionMatrix(predicted, actual, positive = "1")
  
  ppv <- conf_matrix$byClass["Pos Pred Value"]
  npv <- conf_matrix$byClass["Neg Pred Value"]
  sens <- conf_matrix$byClass["Sensitivity"]
  spec <- conf_matrix$byClass["Specificity"]
  f1 <- 2 * (ppv * sens) / (ppv + sens)
  roc_obj <- pROC::roc(true, pred)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  panel_metrics <- rbind(panel_metrics, data.frame(
    Methode = meth,
    AUC = auc_val,
    F1 = f1,
    Sensitivity = sens,
    Specificity = spec,
    PPV = ppv,
    NPV = npv
  ))
}

print(panel_metrics)

library(ggplot2)

ggplot(panel_metrics, aes(x = Methode, y = AUC, fill = Methode)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  labs(
    y = "AUC",
    x = "Normalization Method"
  ) +
  ylim(0, 1.0) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2", name = "Normalization") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(40, 5, 5, 5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1)
  )


# Prepare long format data
library(reshape2)
ppv_npv_long <- melt(panel_metrics, id.vars = "Methode", measure.vars = c("PPV", "NPV"),
                     variable.name = "Metric", value.name = "Value")

ggplot(ppv_npv_long, aes(x = Methode, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  labs(
    y = "Value",
    x = "Normalization Method"
  ) +
  ylim(0, 0.95) +
  theme_bw() +
  scale_fill_manual(values = c("PPV" = "#1f78b4", "NPV" = "#33a02c"),
                    name = "Metric",
                    labels = c("PPV", "NPV")) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(40, 5, 5, 5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1)
  )


# Bootstrapping
library(boot)
library(pROC)
library(caret)
library(dplyr)

norm_methods <- list(
  Raw = peptide_data,
  Log2 = log2_data,
  Quantile = quantile_norm_data,
  Loess = loess_norm_data,
  LinReg = linreg_norm_data,
  VSN = vsn_raw
)

panel_vars <- c("AHSG.1", "CLEC3B.20", "COMP.23.", "F9.27", "LRG1.42", "MCAM.47", "MRC1.73")
panel_bootstrap_results <- list()

for(meth in names(norm_methods)) {
  nm <- norm_methods[[meth]]
  panel_df <- data.frame(
    concentration_calibartion_curve[, c("Row.names", "group", "cleaned")],
    nm[, panel_vars]
  )
  panel_df$group_bin <- ifelse(panel_df$group == "active", 0, 1)
  train_idx <- panel_df$cleaned %in% row.names(data_train)
  test_idx <- panel_df$cleaned %in% row.names(data_test)
  formula_panel <- as.formula(paste("group_bin ~", paste(panel_vars, collapse = " + ")))
  model <- glm(formula_panel, data = panel_df[train_idx, ], family = binomial)
  pred_test <- predict(model, newdata = panel_df[test_idx, ], type = "response")
  
  calib_train <- data.frame(set = "Train", obs = panel_df$group_bin[train_idx],
                            pred = predict(model, newdata = panel_df[train_idx, ], type = "response"))
  calib_test <- data.frame(set = "Test", obs = panel_df$group_bin[test_idx], pred = pred_test)
  calib_all <- bind_rows(calib_train, calib_test)
  
  calculate_boot_ci <- function(data, R = 1000) {
    if (nrow(data) == 0) return(data.frame())
    boot_func <- function(d, indices) { mean(d$obs[indices]) }
    boot_results <- boot(data, statistic = boot_func, R = R)
    ci <- boot.ci(boot_results, type = "perc", conf = 0.95)
    data.frame(ci_low = ci$percent[4], ci_high = ci$percent[5])
  }
  calib_binned <- calib_all %>%
    group_by(set) %>%
    mutate(bin = ntile(pred, 5)) %>%
    group_by(set, bin) %>%
    group_modify(~ {
      bind_cols(
        summarise(.x, mean_pred = mean(pred), obs_rate = mean(obs), n_bin = n()),
        calculate_boot_ci(.x)
      )
    }) %>%
    ungroup()
  
  brier_train <- mean((calib_train$pred - calib_train$obs)^2)
  brier_test  <- mean((calib_test$pred - calib_test$obs)^2)
  roc_train <- roc(response = calib_train$obs, predictor = calib_train$pred)
  roc_test <- roc(response = calib_test$obs, predictor = calib_test$pred)
  auc_train <- as.numeric(auc(roc_train))
  auc_test  <- as.numeric(auc(roc_test))
  
  brier_boot <- boot(calib_test, function(data, i) {
    mean((data$pred[i] - data$obs[i])^2)
  }, R = 1000)
  brier_ci <- quantile(brier_boot$t, c(0.025, 0.975))
  
  panel_bootstrap_results[[meth]] <- list(
    calib_binned = calib_binned,
    brier_train = brier_train,
    brier_test  = brier_test,
    auc_train   = auc_train,
    auc_test    = auc_test,
    brier_ci    = brier_ci
  )
}

# Ergebnisübersicht
result_summary <- do.call(rbind, lapply(names(panel_bootstrap_results), function(m) {
  res <- panel_bootstrap_results[[m]]
  data.frame(
    Methode = m,
    Brier_Train = res$brier_train,
    Brier_Test  = res$brier_test,
    Brier_CI_Low = res$brier_ci[1],
    Brier_CI_High= res$brier_ci[2],
    AUC_Train    = res$auc_train,
    AUC_Test     = res$auc_test
  )
}))
print(result_summary)


library(ggplot2)
library(reshape2)

# 1. Brier score plot (Train/Test, English labels)
brier_long <- melt(result_summary, id.vars = "Methode", measure.vars = c("Brier_Train", "Brier_Test"),
                   variable.name = "Dataset", value.name = "Brier_Score")

ggplot(brier_long, aes(x = Methode, y = Brier_Score, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    y = "Brier Score (lower = better)",
    x = "Normalization Method"
  ) +
  ylim(0, 0.15) +
  theme_bw() +
  scale_fill_manual(
    values = c("Brier_Train" = "#99ccff", "Brier_Test" = "#8835ff"),
    name = "Dataset",
    labels = c("Train", "Test")
  )+
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(40, 5, 5, 5),  # Viel mehr Abstand oben (40 statt 10)
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1)
  )





