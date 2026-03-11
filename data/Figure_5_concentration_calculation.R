# Plasma concentrationm calculated for the complete dataset usig the concentration calibartion curve
# Load necessary packages
library(ROCR)
library(pROC)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(xgboost)
library(SHAPforxgboost)
library(shapviz)
library(ggridges)
library(cowplot)
library(ggbeeswarm)#
library(ResourceSelection)
# Functions

# Nagelkerke R² and DOC
calculate_nagelkerke_r2 <- function(model) {
  null_formula <- formula(paste(names(model$model)[1], "~ 1"))
  null_model <- glm(null_formula, data = model$data, family = "binomial")
  ll_model <- logLik(model)
  ll_null <- logLik(null_model)
  n <- nobs(model)
  r2_cox_snell <- 1 - exp((2/n) * (as.numeric(ll_null) - as.numeric(ll_model)))
  max_r2 <- 1 - exp((2/n) * as.numeric(ll_null))
  r2_nagelkerke <- r2_cox_snell / max_r2
  return(r2_nagelkerke)
}

calculate_doc_robust <- function(actual, predicted_probs) {
  # make sure the input is numeric and NAs are accounted for
  if (!is.numeric(actual)) actual <- as.numeric(as.character(actual))
  if (!is.numeric(predicted_probs)) predicted_probs <- as.numeric(predicted_probs)
  
  complete_cases <- complete.cases(actual, predicted_probs)
  actual <- actual[complete_cases]
  predicted_probs <- predicted_probs[complete_cases]
  
  if (length(actual) == 0) {
    warning("None complete cases for DOC-Calculation.")
    return(NA)
  }
  
  unique_actual <- unique(actual)
  if (length(unique_actual) < 2) {
    warning("Only one class in 'actual' after NA-accouning – DOC undefined")
    return(NA)
  }
  # Assumption: Classes are 0 and 1 after conversion/adjustment
  pos_class <- max(unique_actual)
  neg_class <- min(unique_actual)
  
  pos_probs <- predicted_probs[actual == pos_class]
  neg_probs <- predicted_probs[actual == neg_class]
  
  if (length(pos_probs) == 0 || length(neg_probs) == 0) {
    # sould be covered by the above provided check, but make sure again
    warning("No pairs from both classes available.")
    return(NA)
  }
  
  # Use expand.grid for pairs (more efficient for large data than outer)
  pairs <- expand.grid(pos_probs = pos_probs, neg_probs = neg_probs)
  
  # Count concordant, discordant pairs and ties 
  concordant <- sum(pairs$pos_probs > pairs$neg_probs)
  discordant <- sum(pairs$pos_probs < pairs$neg_probs)
  ties <- sum(pairs$pos_probs == pairs$neg_probs)
  
  # DOC = (concordant + 0.5 * ties) / total number of pairs
  doc <- (concordant + 0.5 * ties) / nrow(pairs)
  
  return(doc)
}


# For the Calibation Based Method
concentration_calibartion_curve <- read.csv2("../Uwes_wishes/Nature_Code/data/concentration_calibartion_curve.csv", row.names = 1)
concentration_calibartion_curve <- data.frame(t(concentration_calibartion_curve))

concentration_calibartion_curve$Experimental_Group <- row.names(concentration_calibartion_curve)
concentration_calibartion_curve <-  concentration_calibartion_curve %>%
  mutate(group = case_when(
    grepl("active", Experimental_Group) ~ "active",
    grepl("Active", Experimental_Group) ~ "active",
    grepl("Remission_", Experimental_Group) ~ "remission",
    grepl("remisson", Experimental_Group) ~ "remission",
    grepl("Remission", Experimental_Group) ~ "remission",
    grepl("HC", Experimental_Group) ~ "healthy",
    TRUE ~ NA_character_ # for words that do not match any of the conditions
  ))



load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_trainn.r")
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_test.r")


# Example: row.names of data_test and concentration_calibartion_curve
# Vectors with the row names
# Create dataframe with original and adjusted row names
# Dataframe with original row names
name_mapping <- data.frame(
  original = row.names(concentration_calibartion_curve),
  cleaned = NA_character_,
  stringsAsFactors = FALSE
)

# Maske für Prag-Row-Names
is_prag <- grepl("^Prag", name_mapping$original)
name_mapping$cleaned[is_prag] <- name_mapping$original[is_prag]

# clean up Berlin
name_mapping$cleaned[!is_prag] <- name_mapping$original[!is_prag] |>
  sub("^Berlin_", "", x = _) |>         # "Berlin_" entfernen
  gsub("\\.", "_", x = _)   |>          # Punkt zu Unterstrich
  sub("(_P\\d+_.*)$", "", x = _)        # alles ab _P... löschen

# look at the results 
head(name_mapping, 10)

rename_mapping <- c(
  "PR3_active_A20_40" = "MPO_active_A20_40",
  "PR3_active_A21_14" = "MPO_active_A21_14" ,
  "PR3_active_A21_6" = "MPO_active_A21_6" ,
  "PR3_active_A20_21" =  "PR3_active_A20_21"
)
row.names(data_test) <- ifelse(row.names(data_test) %in% names(rename_mapping),
                          rename_mapping[row.names(data_test)],
                          row.names(data_test))

row.names(name_mapping) <- name_mapping$original
concentration_calibartion_curve <- merge(concentration_calibartion_curve,name_mapping, by = 0 )

concentration_calibartion_curve_test <- concentration_calibartion_curve[concentration_calibartion_curve$cleaned %in% row.names(data_test), ]

not_in_data_test <- concentration_calibartion_curve_test$cleaned[!(concentration_calibartion_curve_test$cleaned %in% row.names(data_test))]
not_in_conc <- row.names(data_test)[!(row.names(data_test) %in% concentration_calibartion_curve_test$cleaned)]



rename_mapping <- c(
  "PR3_active_A20_40" = "MPO_active_A20_40",
  "PR3_active_A21_14" = "MPO_active_A21_14" ,
  "PR3_active_A21_6" = "MPO_active_A21_6" ,
  "PR3_active_A20_21_1" =  "PR3_active_A20_21"
)
row.names(data_train) <- ifelse(row.names(data_train) %in% names(rename_mapping),
                               rename_mapping[row.names(data_train)],
                               row.names(data_train))

concentration_calibartion_curve_train <- concentration_calibartion_curve[concentration_calibartion_curve$cleaned %in% row.names(data_train), ]

not_in_data_train <- concentration_calibartion_curve_train$cleaned[!(concentration_calibartion_curve_train$cleaned %in% row.names(data_train))]
not_in_conc <- row.names(data_train)[!(row.names(data_train) %in% concentration_calibartion_curve_train$cleaned)]


concentration_calibartion_curve_train$factorgroup_numeric <- factor(concentration_calibartion_curve_train$group , 
levels = c("active", "remission"), 
labels = c(0, 1))


concentration_calibartion_curve_test$group <- factor(concentration_calibartion_curve_test$group , 
                                                      levels = c("active", "remission"), 
                                                      labels = c(0, 1))



formula_model4   <- group ~ AHSG.1 + CLEC3B.20 + COMP.23. + F9.27 + LRG1.42 + MCAM.47 + MRC1.73


# Train on 60% train data
model_7panel <- glm(formula_model4, data = concentration_calibartion_curve_train, family = binomial())

# Prediction 
pred_7panel <- predict(model_7panel, newdata = concentration_calibartion_curve_test, type = "response")
# 7Panel
actual_7 <- factor(as.numeric(as.character(concentration_calibartion_curve_test$group)), levels = c(0, 1))
predicted_7 <- ifelse(pred_7panel > 0.5, 1, 0)
predicted_7 <- factor(predicted_7, levels = c(0, 1))

conf_matrix_7 <- caret::confusionMatrix(predicted_7, actual_7, positive = "1")

ppv_7 <- conf_matrix_7$byClass["Pos Pred Value"]
npv_7 <- conf_matrix_7$byClass["Neg Pred Value"]
sens_7 <- conf_matrix_7$byClass["Sensitivity"]
spec_7 <- conf_matrix_7$byClass["Specificity"]
f1_7 <- 2 * (ppv_7 * sens_7) / (ppv_7 + sens_7)

roc_7 <- pROC::roc(response = as.numeric(actual_7), predictor = pred_7panel)
auc_7 <- pROC::auc(roc_7)
r2_7 <- calculate_nagelkerke_r2(model_7panel)
doc_7 <- calculate_doc_robust(as.numeric(actual_7), pred_7panel)

results_df_7 <- data.frame(
  Metric = c("AUC", "PPV", "Sensitivity", "Specificity", "NPV", "F1 Score", "Nagelkerke R²", "DOC"),
  Value = c(round(auc_7, 3), round(ppv_7, 3), round(sens_7, 3), round(spec_7, 3),
            round(npv_7, 3), round(f1_7, 3), round(r2_7, 3), round(doc_7, 3))
)

#openxlsx::write.xlsx(results_df_7, "../Revision/Prep_dataframes/panelcalibrationCurveCOncentration_60_40.xlsx")


plot(roc_7, col = "blue", lwd = 2, main = "ROC – 7er Panel concentration calibartion curve", print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")

### look into more 
formula_model4   <- group ~ AHSG.1 + CLEC3B.20 + COMP.23. + F9.27 + LRG1.42 + MCAM.47 + MRC1.73
model_7panel_ratio <- glm(formula_model4, data = concentration_calibartion_curve_train, family = binomial)

# coefficient plot — it shows the weight (impact) of each protein on the log-odds of remission.
# Extract and format coefficients
coefs <- as.data.frame(coef(summary(model_7panel_ratio)))
coefs$Feature <- rownames(coefs)
colnames(coefs) <- c("Estimate", "StdError", "Zvalue", "Pvalue", "Feature")

# Remove intercept
coefs <- coefs[coefs$Feature != "(Intercept)", ]

# Order by effect size
coefs$Feature <- factor(coefs$Feature, levels = coefs$Feature[order(coefs$Estimate)])

# coefs
# Estimate   StdError      Zvalue      Pvalue   Feature
# AHSG.1     0.004002926 0.07906430  0.05062874 0.959621359    AHSG.1
# CLEC3B.20  0.307848808 0.16219039  1.89807053 0.057686791 CLEC3B.20
# COMP.23.   0.376838297 0.46548342  0.80956330 0.418191208  COMP.23.
# F9.27     -9.043819768 4.13654113 -2.18632415 0.028791898     F9.27
# LRG1.42   -0.047656653 0.02680739 -1.77774338 0.075446011   LRG1.42
# MCAM.47    8.561102306 3.18163400  2.69078791 0.007128350   MCAM.47
# MRC1.73   -0.853197924 0.29473972 -2.89475044 0.003794601   MRC1.73

# Desired order (reversed)
desired_order <- rev(c("MRC1.73", "MCAM.47", "LRG1.42", "F9.27", 
                       "COMP.23.", "CLEC3B.20", "AHSG.1"))

# Apply the reversed order
coefs$Feature <- factor(coefs$Feature, levels = desired_order)
# Plot
library(ggplot2)

# mark significance
coefs <- coefs %>%
  mutate(
    Significance = case_when(
      Pvalue < 0.001 ~ "***",
      Pvalue < 0.01 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Direction = ifelse(Estimate > 0, "Positive", "Negative"),
    SigColor = case_when(
      Pvalue < 0.05 & Estimate > 0 ~ "#023e8a",  # Grün
      Pvalue < 0.05 & Estimate < 0 ~ "#ed2923",  # Rot
      TRUE ~ "#023e8a"
    )
  )

# Plot
coefs$Feature <- gsub("_", " ", coefs$Feature)

coef_plot <- ggplot(coefs, aes(x = Feature, y = Estimate)) +
  geom_col(aes(fill = SigColor), color = "black", width = 0.7,  alpha = 0.7) +
  geom_errorbar(aes(
    ymin = Estimate - 1.96 * StdError,
    ymax = Estimate + 1.96 * StdError
  ), width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = Significance), hjust = ifelse(coefs$Estimate > 0, -0.3, 1.3), size = 6) +
  coord_flip() +
  scale_fill_identity() +
  labs(
    y = "Log-Odds Estimate (Remission ↑)",
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(hjust = 0.5)
  )

# test for normal distribution - if yes = zscore otherwise normalisation tests
# Analyse data distribution (assumption: concentration_calibation_based is the dataframe)

# Shapiro-Wilk test for all peptides (original and log data)
shapiro_results <- data.frame(
  Peptid = character(),
  P_Wert_Original = numeric(),
  P_Wert_Log = numeric(),
  stringsAsFactors = FALSE
)

ar(mfrow = c(2, 2)) # 2x2 plot grids per peptide

for(peptid in colnames(concentration_calibartion_curve)[2:8]) {
  # original data analysis 
  original_data <- concentration_calibartion_curve[[peptid]]
  log_data <- log10(original_data + 1)
  
  # Plots
  hist(original_data, main = paste("Distribution", peptid), 
       xlab = "Concentration", col = "lightblue", breaks = 5)
  qqnorm(original_data, main = paste("Q-Q Plot", peptid))
  qqline(original_data, col = "red")
  boxplot(original_data, main = paste("Boxplot", peptid), horizontal = TRUE, col = "steelblue")
  plot(density(log_data, na.rm = TRUE), main = paste("Log-Transformiert", peptid), col = "darkgreen")
  
  # Shapiro-Wilk-Tests
  sw_original <- shapiro.test(original_data)$p.value
  sw_log <- shapiro.test(log_data)$p.value
  
  # save results 
  shapiro_results <- rbind(shapiro_results, data.frame(
    Peptid = peptid,
    P_Wert_Original = round(sw_original, 4),
    P_Wert_Log = round(sw_log, 4)
  ))
  

  cat(paste0(
    "Peptid: ", peptid, "\n",
    "  Shapiro-Wilk (Original): p = ", round(sw_original, 4), 
    ifelse(sw_original < 0.05, " (none normal)", " (normality possible)"), "\n",
    "  Shapiro-Wilk (Log):     p = ", round(sw_log, 4),
    ifelse(sw_log < 0.05, " (none normal)", " (normality possible)"), "\n\n"
  ))
}

# Finale Ergebnistabelle
print(shapiro_results)
# Peptid: F9.27
# Shapiro-Wilk (Original): p = 0 (none normal)
# Shapiro-Wilk (Log):     p = 0 (none normal)
# 
# Peptid: LRG1.42
# Shapiro-Wilk (Original): p = 7e-04 (none normal)
# Shapiro-Wilk (Log):     p = 0.068 (normality possible)
# 
# Peptid: MCAM.47
# Shapiro-Wilk (Original): p = 0 (none normal)
# Shapiro-Wilk (Log):     p = 0.0015 (none normal)
# 
# Peptid: COMP.23.
# Shapiro-Wilk (Original): p = 7e-04 (none normal)
# Shapiro-Wilk (Log):     p = 0.1837 (normality possible)
# 
# Peptid: AHSG.1
# Shapiro-Wilk (Original): p = 0 (none normal)
# Shapiro-Wilk (Log):     p = 0 (none normal)
# 
# Peptid: MRC1.73
# Shapiro-Wilk (Original): p = 0 (none normal)
# Shapiro-Wilk (Log):     p = 0 (none normal)
# 
# Peptid: CLEC3B.20
# Shapiro-Wilk (Original): p = 0.6346 (normality possible)
# Shapiro-Wilk (Log):     p = 0.0033 (none normal)

par(mfrow = c(2, 2))  

print(shapiro_results)
# Peptid P_Wert_Original P_Wert_Log
# 1     F9.27          0.0000     0.0000
# 2   LRG1.42          0.0007     0.0680
# 3   MCAM.47          0.0000     0.0015
# 4  COMP.23.          0.0007     0.1837
# 5    AHSG.1          0.0000     0.0000
# 6   MRC1.73          0.0000     0.0000
# 7 CLEC3B.20          0.6346     0.0033
# > 

# Identify peptide columns (columns 2-8 based on head()) 
peptide_columns <- colnames(concentration_calibartion_curve)[2:8]

for(peptid in peptide_columns) {
  # Histogram with density curve
  hist(concentration_calibartion_curve[[peptid]],
       main = paste("Distibution", peptid),
       xlab = "Concentration",
       col = "lightblue",
       breaks = 5) # Fewer bins for small samples 
  
  # Q-Q plot for normal distribution check
  qqnorm(concentration_calibartion_curve[[peptid]], main = paste("Q-Q Plot", peptid))
  qqline(concentration_calibartion_curve[[peptid]], col = "red")
  
  # Boxplot with outlier detection
  boxplot(concentration_calibartion_curve[[peptid]],
          main = paste("Boxplot", peptid),
          horizontal = TRUE,
          col = "steelblue")
  
  # Log-transformed version (if required)
  log_values <- log10(concentration_calibartion_curve[[peptid]] + 1)  # +1 to avoid 0 
  plot(density(log_values, na.rm = TRUE),
       main = paste("Log-Transformed", peptid),
       col = "darkgreen")
}


# See 00_test_different_normalisation_calibration.R there different normalisations where tested 
# All four methods (Raw, Log2, Quantile, VSN) perform excellently in your panel, 
# but VSN offers maximum methodological transparency, robustness, and is tailor-made for proteomics normalization. 
library(vsn)

# VSN normalisation on raw data (columns 2-8)
vsn_raw <- justvsn(as.matrix(concentration_calibartion_curve[, 2:8]))

# add metdata
final_data_raw <- cbind(
  concentration_calibartion_curve[, c("Row.names", "group", "cleaned")],
  vsn_raw
)

# check quality
meanSdPlot(vsn_raw)  # sould be more or less flatt-ish
boxplot(vsn_raw, las = 2, main = "VSN-normalised raw data")


# Daten from VSN prep for ggplot
plot_data <- data.frame(
  mean_rank = 1:length(apply(vsn_raw, 1, mean)),
  sd = apply(vsn_raw, 1, sd)
)


# improved Mean-SD plot
ggplot(plot_data, aes(x = mean_rank, y = sd)) +
  geom_point(alpha = 0.7, size = 4, pch = 21, fill = "grey") +
  geom_smooth(color = "red", method = "loess", se = FALSE) +
  labs(
    x = "Mean Rank", 
    y = "Standard Deviation") +
  theme_bw()+
  theme (axis.title.x = element_text (size = 13), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position="right") 

# prepare data for ggplot -> long format 
vsn_long <- reshape2::melt(vsn_raw, 
                           id.vars = NULL, 
                           variable.name = "Peptid", 
                           value.name = "normalised_value")

vsn_long$Var2 <- gsub("_", "", vsn_long$Var2)

# 
ggplot(vsn_long, aes(x = Var2, y = Normalisierter_Wert)) +
  geom_point(fill = "blue",
             pch = 21,                 # filled circle with outline
             alpha = 0.5, 
             size = 2,# soft fill            # bold outline
             position = position_jitter()
  )+
  geom_boxplot(alpha = 0.2, fill = "blue", outliers = F) +
  
  labs(
    y = "VSN normalised value") +
  theme_bw()+
  theme (axis.title.y = element_text (size = 13), 
         axis.text.x = element_text (size = 13, angle = 45), 
         axis.text.y = element_text (size = 13), 
         axis.title.x = element_blank (),
         legend.position="right") ->vsnNorm

# The plots confirm that the VSN normalisation of the data set was very effective. 
# The data is now homogeneous, comparable and optimally prepared for further analyses or machine learning.
# VSN turns ‘unstable’ and difficult to compare proteomics data into ‘calm’ data,
# uniformly comparable values - an important prerequisite for reliable analyses


# write.csv(final_data_raw, file = "../Uwes_wishes/Nature_Code/data/variance_final_data_raw.csv", quote = F)

vsn_data_based_test <- final_data_raw[final_data_raw$cleaned %in% row.names(data_test), ]
vsn_data_based_train <- final_data_raw[final_data_raw$cleaned %in% row.names(data_train), ]

# 1) Make a clean numeric outcome
vsn_data_based_train$group01 <- ifelse(vsn_data_based_train$group == "remission", 1L, 0L)
vsn_data_based_test$group01  <- ifelse(vsn_data_based_test$group == "remission", 1L, 0L)

table(vsn_data_based_train$group01, useNA = "ifany")  # sanity check

# 2) Use group01 in the formula
formula_model4 <- group01 ~ AHSG.1 + CLEC3B.20 + COMP.23. + F9.27 + LRG1.42 + MCAM.47 + MRC1.73

model_7panel_vsn <- glm(formula_model4,
                        data = vsn_data_based_train,
                        family = binomial())


# coefficient plot — it shows the weight (impact) of each protein on the log-odds of remission.

# Extract and format coefficients
coefs <- as.data.frame(coef(summary(model_7panel_vsn)))
coefs$Feature <- rownames(coefs)
colnames(coefs) <- c("Estimate", "StdError", "Zvalue", "Pvalue", "Feature")

# Remove intercept
coefs <- coefs[coefs$Feature != "(Intercept)", ]

# Order by effect size
coefs$Feature <- factor(coefs$Feature, levels = coefs$Feature[order(coefs$Estimate)])

# Desired order (reversed)
desired_order <- rev(c("MRC1.73", "MCAM.47", "LRG1.42", "F9.27", 
                       "COMP.23.", "CLEC3B.20", "AHSG.1"))

# Apply the reversed order
coefs$Feature <- factor(coefs$Feature, levels = desired_order)
# coefs
# Estimate  StdError    Zvalue      Pvalue   Feature
# AHSG.1    -0.9013899 1.1564301 -0.779459 0.435709364    AHSG.1
# CLEC3B.20  3.8245496 1.3474988  2.838258 0.004536048 CLEC3B.20
# COMP.23.   2.9994271 1.2319268  2.434745 0.014902303  COMP.23.
# F9.27     -1.9770400 1.2629601 -1.565402 0.117488743     F9.27
# LRG1.42    1.2955697 1.0498390  1.234065 0.217178636   LRG1.42
# MCAM.47    0.9758509 0.9523317  1.024696 0.305506461   MCAM.47
# MRC1.73   -4.6909732 1.8628482 -2.518172 0.011796553   MRC1.73


coefs <- coefs %>%
  mutate(
    Significance = case_when(
      Pvalue < 0.001 ~ "***",
      Pvalue < 0.01 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Direction = ifelse(Estimate > 0, "Positive", "Negative"),
    SigColor = case_when(
      Pvalue < 0.05 & Estimate > 0 ~ "#023e8a",  # Grün
      Pvalue < 0.05 & Estimate < 0 ~ "#ed2923",  # Rot
      TRUE ~ "#023e8a"
    )
  )

# Plot
coefs$Feature <- gsub("_", " ", coefs$Feature)

coef_plot <- ggplot(coefs, aes(x = Feature, y = Estimate)) +
  geom_col(aes(fill = SigColor), color = "black", width = 0.7,  alpha = 0.7) +
  geom_errorbar(aes(
    ymin = Estimate - 1.96 * StdError,
    ymax = Estimate + 1.96 * StdError
  ), width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = Significance), hjust = ifelse(coefs$Estimate > 0, -0.3, 1.3), size = 6) +
  coord_flip() +
  scale_fill_identity() +
  labs(
    y = "Log-Odds Estimate (Remission ↑)",
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(hjust = 0.5)
  )


# Prediction 
pred_7panel_vsn <- predict(formula_model4, newdata = vsn_data_based_test, type = "response")

# 7Panel
actual_7_vsn <- factor(as.numeric(as.character(vsn_data_based_test$group)), levels = c(0, 1))
predicted_7_vsn <- ifelse(pred_7panel_vsn > 0.5, 1, 0)
predicted_7_vsn <- factor(predicted_7_vsn, levels = c(0, 1))

conf_matrix_7_vsn <- caret::confusionMatrix(predicted_7_vsn, actual_7_vsn, positive = "1")

ppv_7 <- conf_matrix_7_vsn$byClass["Pos Pred Value"]
npv_7 <- conf_matrix_7_vsn$byClass["Neg Pred Value"]
sens_7 <- conf_matrix_7_vsn$byClass["Sensitivity"]
spec_7 <- conf_matrix_7_vsn$byClass["Specificity"]
f1_7 <- 2 * (ppv_7 * sens_7) / (ppv_7 + sens_7)

roc_7 <- pROC::roc(response = as.numeric(actual_7_vsn), predictor = pred_7panel_vsn)
auc_7 <- pROC::auc(roc_7)
r2_7 <- calculate_nagelkerke_r2(model_7panel_vsn)
doc_7 <- calculate_doc_robust(as.numeric(actual_7_vsn), actual_7_vsn)

results_df_7 <- data.frame(
  Metric = c("AUC", "PPV", "Sensitivity", "Specificity", "NPV", "F1 Score", "Nagelkerke R²", "DOC"),
  Value = c(round(auc_7, 3), round(ppv_7, 3), round(sens_7, 3), round(spec_7, 3),
            round(npv_7, 3), round(f1_7, 3), round(r2_7, 3), round(doc_7, 3))
)


# openxlsx::write.xlsx(results_df_7, "Revision/Prep_dataframes/7panelVSNNormalisationcalibartion60_40.xlsx")


plot(roc_7, col =  "#c70000", lwd = 2, print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")


# costamized Plot for Uwe
all_samples_df <- data.frame(
  Sample = vsn_data_based_test$Row.names,
  Predicted = predicted_7_vsn,
  Actual = vsn_data_based_test$group,
  Probability = pred_7panel_vsn
)

all_class <- all_samples_df %>%
  mutate(
    Actual = as.numeric(as.character(Actual)),
    Predicted = as.numeric(as.character(Predicted)),
    Jittered_Actual = Actual + runif(n(), -0.4, 0.4),
    Jittered_Predicted = Predicted + runif(n(), -0.4, 0.4)
  )

uwe_plot <- ggplot(all_class, aes(x = Jittered_Actual, y = Jittered_Predicted, fill = as.factor(Predicted))) +
  geom_point(size = 4, alpha = 0.8, pch = 21, color = "black") +
  scale_fill_manual(values = c("#ffd973", "#0000ab"))+
  scale_x_continuous(breaks = c(0, 1), labels = c("Active", "Remission"), limits = c(-0.5, 1.5)) +
  scale_y_continuous(breaks = c(0, 1), labels = c("Active", "Remission"), limits = c(-0.5, 1.5)) +
  labs(
    x = "Actual Class",
    y = "Predicted Class",
    fill = "Prediction"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.grid = element_blank()
  )

# ggsave(uwe_plot, file = "Revision/Figures_RevisionCalibration/Corrected_Calibartion60_40confusionMatrix_glm_UweStyle.svg", width = 5, height = 4)

ˇ# Lets make the shap-ish costamized plot

get_glm_beeswarm_data_v2 <- function(model, data, feature_list) {
  pred_probs <- predict(model, newdata = data, type = "response")
  
  data_list <- lapply(feature_list, function(feature) {
    data.frame(
      PatientID = rownames(data),
      Protein = feature,
      Value = data[[feature]],
      Prediction = pred_probs,
      Group = factor(data$group, labels = c("active", "remission"))
    )
  })
  
  do.call(rbind, data_list)
}


# Define your protein list
proteins <- c("MRC1.73", "MCAM.47", "LRG1.42", "F9.27", 
              "COMP.23.", "CLEC3B.20", "AHSG.1")

# Run it
model_7panel_vsn <- glm(formula_model4, data = vsn_data_based_train, family = binomial)
beeswarm_data <- get_glm_beeswarm_data_v2(model_7panel_vsn, vsn_data_based_test, proteins)



median_lines_grouped <- beeswarm_data %>%
  group_by(Protein, Group) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = "drop")

# write.xlsx(median_lines_grouped, "Revision/intermediate/Corrected_60_40_median_lines_vsncalibartion.xlsx")



beeswarm_data$PredictedGroup <- ifelse(beeswarm_data$Prediction > 0.5, "Predicted Remission", "Predicted Active")
beeswarm_data$PredictedGroup <- factor(beeswarm_data$PredictedGroup, levels = c("Predicted Active", "Predicted Remission"))

median_lines_predicted <- beeswarm_data %>%
  group_by(Protein, PredictedGroup) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = "drop")

# write.xlsx(median_lines_predicted, "Revision/intermediate/Corrected_60_40_median_lines_vsncalibartion_predicted40.xlsx")


# beeswarm_data$Protein <- gsub("_", " ", beeswarm_data$Protein)
# beeswarm_data$Protein <- gsub("\\..", "", beeswarm_data$Protein)
ggplot(beeswarm_data, aes(x = Value, y = Protein, fill = PredictedGroup)) +
  # ggdist::stat_halfeye(
  #   adjust = 0.5,       # smoothness of the density
  #   justification = -0.3, # move density slightly to the left
  #   .width = 0,
  #   point_colour = NA,
  #   alpha = 0.6
  # ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = PredictedGroup),
    size = 2,
    alpha = 0.7,
    width = 0.2,
    groupOnX = FALSE
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7)  +
  # acet_wrap(~ PredictedGroup) +
  scale_fill_manual(values = c("Predicted Active" = "#cc5500", "Predicted Remission" = "#a2ab2d")) +
  scale_color_manual(values = c("Predicted Active" = "#cc5500", "Predicted Remission" = "#a2ab2d")) +
  labs(
    title = "60% train - 40% test",
    x = "TQL ratio",
    y = NULL,
    fill = "True Group",
    color = "True Group"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15)
  )



# Zuerst: rownames als Spalte in den Datensatz bringen
vsn_data_based_test$PatientID <- vsn_data_based_test$Row.names
example_patient_ids <- c("Berlin_MPO_active_A19.24_P02_F01", "Prag_PR3_Remission_3784_P02_E04")

example_long <- vsn_data_based_test %>%
  filter(PatientID %in% example_patient_ids) %>%
  select(PatientID, group, all_of(proteins)) %>%
  tidyr::pivot_longer(cols = all_of(proteins), names_to = "Protein", values_to = "Value") %>%
  mutate(
    Protein = factor(Protein, levels = rev(proteins)),
    Group = factor(group, levels = c(0, 1), labels = c("Active", "Remission")),
    Prediction = predict(model_7panel_vsn, 
                         newdata = vsn_data_based_test[vsn_data_based_test$PatientID %in% example_patient_ids, ], type = "response")[rep(1:2, each = length(proteins))]
  )

example_long$PredictedGroup <- ifelse(example_long$Prediction > 0.5, 
                                      "Predicted Remission", "Predicted Active")


# beeswarm_data$Protein <- gsub("\\.", "", beeswarm_data$Protein)
# example_long$Protein <- gsub("\\.", "", example_long$Protein)
# example_long$Protein <- gsub("[0-9]+.*$", "", example_long$Protein)
# beeswarm_data$Protein <- gsub("[0-9]+.*$", "", beeswarm_data$Protein)

ggplot(beeswarm_data, aes(x = Value, y = Protein, fill = PredictedGroup)) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -0.3,
    .width = 0,
    point_colour = NA,
    alpha = 0.6
  ) +
  # geom_boxplot(
  #   width = 0.15,
  #   outlier.shape = NA,
  #   alpha = 0.5,
  #   position = position_nudge(y = -0.2)
  # ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = PredictedGroup),
    size = 1.5,
    alpha = 0.7,
    width = 0.2,
    groupOnX = FALSE
  ) +
  # patient-line
  geom_path(
    data = example_long,
    aes(x = Value, y = Protein, group = PatientID,  color = PredictedGroup), linewidth = 0.5
  ) +
  geom_point(
    data = example_long,
    aes(x = Value, y = Protein, fill = Group),
    shape = 21, color = "black", size = 4
  ) +
  scale_fill_manual(values = c(
    "Predicted Active" = "#ffd973",
    "Predicted Remission" = "#0000ab",
    "Active" = "#ffd973",
    "Remission" = "#0000ab"
  )) +
  scale_color_manual(values = c(
    "Predicted Active" = "#ffd973",
    "Predicted Remission" = "#0000ab"
  )) +
  labs(
    x = "VSN-normalized peptide concentration",
    y = NULL,
    fill = "Group",
    color = "Group"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15)
  )


# rownames to columns
concentration_calibartion_curve$PatientID <- concentration_calibartion_curve$Row.names
example_patient_ids <- c("Berlin_MPO_active_A19.24_P02_F01", "Berlin_PR3_Remission_1826_P02_B03")

example_long <- concentration_calibartion_curve %>%
  filter(PatientID %in% example_patient_ids) %>%
  select(PatientID, group, all_of(proteins)) %>%
  tidyr::pivot_longer(cols = all_of(proteins), names_to = "Protein", values_to = "Value")

# example_long$Protein <- gsub("_", " ", example_long$Protein)
# example_long$Protein <- gsub("\\.", "", example_long$Protein)


# write.csv2(example_long, file = "../Revision/Prep_dataframes/Example_samples_ForMarie.csv")


# Predict for all patients (e.g., test set)
all_probs <- predict(model_7panel_vsn, newdata = vsn_data_based_test, type = "response")

prediction_df <- data.frame(
  PatientID = vsn_data_based_test$Row.names,
  Prediction = all_probs,
  TrueGroup = factor(vsn_data_based_test$group, labels = c("Active", "Remission"))
) %>%
  arrange(Prediction) %>%
  mutate(Index = 1:n())  # for x-axis

highlight_ids  <- c("Berlin_MPO_active_A19.24_P02_F01", "Berlin_PR3_Remission_1826_P02_B03")

highlight_data <- prediction_df[prediction_df$PatientID %in% highlight_ids, ]

ggplot(prediction_df, aes(x = Index, y = Prediction)) +
  geom_line(color = "gray80") +
  geom_point(color = "gray60", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  
  # Highlight patients
  geom_point(data = highlight_data, aes(x = Index, y = Prediction, fill = TrueGroup),
             color = "black", size = 3, shape = 21, stroke = 1.5) +
  
  geom_text(data = highlight_data, aes(x = Index, y = Prediction + 0.05,
                                       label = paste0("Pat ", highlight_ids, "\n", round(Prediction, 2))),
            color = "black", size = 3.5, vjust = 0) +
  
  scale_fill_manual(values = c("Active" = "#ffd973", "Remission" = "#0000ab")) +
  labs(
    title = "Model Predictions for Example Patients",
    y = "Predicted Probability of Remission",
    x = "Patients (sorted by prediction)",
    fill = "True Class"
  ) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
  ) -> second


# Main Plot
main_plot <- ggplot(prediction_df, aes(x = Index, y = Prediction)) +
  geom_line(color = "gray80") +
  geom_point(color = "gray60", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  
  # Highlight patients
  geom_point(data = highlight_data, aes(x = Index, y = Prediction, fill = TrueGroup),
             color = "black", size = 3, shape = 21, stroke = 1.5) +
  
  geom_text(data = highlight_data, aes(x = Index, y = Prediction + 0.05,
                                       label = paste0("Pat ", highlight_ids, "\n", round(Prediction, 2))),
            color = "black", size = 3.5, vjust = 0) +
  
  scale_fill_manual(values = c("Active" = "#ffd973", "Remission" = "#0000ab")) +
  labs(
    title = "Model Predictions for Example Patients",
    y = "Predicted Probability of Remission",
    x =  "Patients (sorted by prediction)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 20)  # Mehr Abstand unten (20pt)
  )

# Separater Rug-Plot test
rug_plot <- ggplot(prediction_df, aes(x = Index, color = TrueGroup)) +
  # longer Rug-Linien using higher length-values
  geom_rug(sides = "top", length = unit(2, "npc"), size = 0.2) +
  scale_color_manual(values = c("Active" = "#ffd973", "Remission" = "#0000ab")) +
  labs(x = "Patients (sorted by prediction)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10), 
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  )


rug_plot <- ggplot(prediction_df, aes(x = Index, color = TrueGroup)) +
  geom_rug(sides = "top", length = unit(2, "npc"), size = 1) +
  scale_color_manual(values = c("Active" = "#ffd973", "Remission" = "#0000ab")) +
  labs(x = "Patients (sorted by prediction)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(40, 5, 5, 5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  )


# Combine plots
combined_plot <- plot_grid(
  main_plot, rug_plot, 
  ncol = 1,
  align = "v",
  rel_heights = c(6, 1),  
  axis = "lr",
  # increase distance between the plots
  greedy = FALSE
)




