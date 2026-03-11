# ============================================
# 0) Packages
# ============================================
library(glmnet)     # LASSO
library(pROC)       # ROC + AUC CI
library(boot)       # bootstrap helper
library(caret)      # confusionMatrix, createDataPartition if needed
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)  
# ============================================
# 1) Load data
# ============================================

reduction_csv <- "../Uwes_wishes/Nature_Code/data/135_confoundR_cleaned.csv"   # adjust if located elsewhere
proteome_file  <- "../Uwes_wishes/Nature_Code/data/605_data_unique_genenames.tsv"  # proteome data (wide)
meta_file      <- "../Uwes_wishes/Nature_Code/data/all_metadata.tsv"   # metadata
outdir         <- "../feature_reduction"
dir.create(outdir, showWarnings = FALSE)

# 2) Load reduction list (you already printed it in console)
reduction <- read.csv2(reduction_csv, stringsAsFactors = FALSE)
# drop empty last row if present
reduction <- reduction[!is.na(reduction$Proteins) & nzchar(reduction$Proteins), , drop = FALSE]
candidate_pre <- as.character(reduction$Proteins)
cat("Pre-LASSO candidate count:", length(candidate_pre), "\n")

# 3) Load proteome and metadata (canonicalise names)
df.rf <- read.table(proteome_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
df.rf <- data.frame(t(df.rf), check.names = FALSE)
# canonicalise rownames in proteome like your earlier code
rownames(df.rf) <- gsub("_", ".", rownames(df.rf))

meta <- read.table(meta_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
meta <- meta %>% filter(!group == "HC")
rownames(meta) <- meta$Sample_ID_Proteomics
meta$Sample_ID_Proteomics <- NULL

# Safe canonicalisation: keep only letters and digits, convert to lower-case
canon <- function(x) {
  # coerce to character and remove any non-alphanumeric characters
  x2 <- as.character(x)
  x2 <- tolower(gsub("[^A-Za-z0-9]", "", x2))
  return(x2)
}

# Quick test
canon(c("IGHV1.24","IGHV1_24","IGHV1-24","IGHV1 24"))
# expected: "ighv124" for all

# assume `reduction` is already loaded (your pro_139_sig -> reduction)
candidate_pre <- as.character(reduction$Proteins)
# load proteome column names (adjust path if needed)
df.rf <- data.frame(t(df.rf), check.names = FALSE)
rownames(df.rf) <- gsub("_", ".", rownames(df.rf))
prot_cols <- colnames(df.rf)

cand_map <- data.frame(candidate = candidate_pre, matched = NA_character_, stringsAsFactors = FALSE)
for(i in seq_along(candidate_pre)) {
  idx <- which(canon(candidate_pre[i]) == canon(prot_cols))
  if(length(idx) >= 1) cand_map$matched[i] <- prot_cols[idx[1]]
}
write.csv(cand_map, "candidate_to_proteome_mapping_fixed.csv", row.names = FALSE)
print(head(cand_map, 20))
if(any(is.na(cand_map$matched))) {
  warning("Some candidates not matched. Inspect candidate_to_proteome_mapping_fixed.csv to correct names.")
} else {
  message("All candidates matched. Proceed with the rest of the script.")
}

# assume cand_map exists in workspace and df.rf and meta loaded earlier
matched_cols <- cand_map$matched
df_reduced <- df.rf[, matched_cols, drop = FALSE]
df_reduced <- data.frame(df_reduced, check.names = FALSE)

row.names(meta) <- gsub("_", ".", row.names(meta))
df_reduced <- df_reduced[row.names(df_reduced) %in% row.names(meta),]
meta_need <- meta[, c("group", "disease.nr")]

df_reduced <- merge(df_reduced, meta_need, by = 0)
row.names(df_reduced) <- df_reduced$Row.names
df_reduced$Row.names <- NULL
# quick checks
cat("Reduced data dim (samples x features+meta):", dim(df_reduced), "\n")
cat("Sample groups distribution:\n"); print(table(df_reduced$group))
# Sample groups distribution:
#   
#   MPO_A MPO_R PR3_A PR3_R 
# 22    16    36    39 

# convert group to numeric factor 0/1 as in pipeline
df_reduced$group <- factor(df_reduced$disease.nr, levels = c("active","remission"), labels = c(0,1))
# keep only matched predictors + group
model_df <- df_reduced[, c(matched_cols, "group")]

# create train/test if not already saved
train_rds <- "../Uwes_wishes/Nature_Code/data/Lasso80train.data.r"  # your .r files
test_rds  <- "../Uwes_wishes/Nature_Code/data/Lasso80test.data.r"
if(file.exists(train_rds) && file.exists(test_rds)) {
  load(train_rds); load(test_rds)
  cat("Loaded existing train/test splits.\n")
} else {
  set.seed(123)
  library(caret)
  idx <- createDataPartition(model_df$group, p = 0.8, list = FALSE)
  train.data <- model_df[idx, , drop = FALSE]
  test.data  <- model_df[-idx, , drop = FALSE]
  save(train.data, file = train_rds)
  save(test.data, file = test_rds)
  cat("Saved train/test splits to:\n", train_rds, "\n", test_rds, "\n")
}
cat("Train/test sizes:", nrow(train.data), nrow(test.data), "\n")
# Train/test sizes: 91 22
na <- row.names(train.data)
nat <- row.names(test.data)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# ============================================
# 3) Lasso on global
# ============================================
load(file = "../Uwes_wishes/Nature_Code/data/Lasso80train.data.r")
load(file = "../Uwes_wishes/Nature_Code/data/Lasso80test.data.r")
# in_panel21 <- "../Uwes_wishes/Nature_Code/data/important_features_sorted_lasso_21.csv"
outdir <- "../feature_reduction/lasso"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot("group" %in% colnames(train.data))

# ============================================
# 4) Outcome y (0/1) with checks
# ============================================
# Common pitfall: as.numeric(factor) gives 1/2 not 0/1. Avoid that.
# This handles numeric 0/1, character "0"/"1", and factor levels.

y_raw <- train.data$group

if (is.factor(y_raw)) {
  y_raw <- as.character(y_raw)
}
if (is.character(y_raw)) {
  y_train <- suppressWarnings(as.integer(y_raw))
} else {
  y_train <- as.integer(y_raw)
}

stopifnot(length(y_train) == nrow(train.data))
stopifnot(all(y_train %in% c(0L, 1L)))
# Optional: ensure both classes exist (otherwise binomial fit can fail or be degenerate)
stopifnot(length(unique(y_train)) == 2)

# ============================================
# 5) Predictors + design matrix
# ============================================
predictors_all <- setdiff(colnames(train.data), "group")
stopifnot(length(predictors_all) >= 1)

x_train <- model.matrix(~ . - 1, data = train.data[, predictors_all, drop = FALSE])
x_test  <- model.matrix(~ . - 1, data = test.data[,  predictors_all, drop = FALSE])

# 2) outcome as 0/1
y_train <- as.integer(as.character(train.data$group))
y_test  <- as.integer(as.character(test.data$group))
stopifnot(all(y_train %in% c(0L,1L)))
stopifnot(all(y_test  %in% c(0L,1L)))

pf_all1 <- rep(1, ncol(x_train))
names(pf_all1) <- colnames(x_train)

# 5) CV + final fit (standard LASSO)
K <- 10
set.seed(42)
foldid <- sample(rep(1:K, length.out = nrow(x_train)))

cv_lasso_all <- glmnet::cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  foldid = foldid,
  nfolds = K,
  type.measure = "deviance",
  standardize = TRUE
  # penalty.factor omitted -> default is all 1
  # or add: , penalty.factor = pf_all1
)

cat("lambda.min:", cv_lasso_all$lambda.min,
    "lambda.1se:", cv_lasso_all$lambda.1se, "\n")

lasso_all_final <- glmnet::glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  lambda = cv_lasso_all$lambda.min,
  standardize = TRUE
  # penalty.factor omitted -> default is all 1
  # or add: , penalty.factor = pf_all1
)

# 6) extract coefficients & selected features
coefs_mat <- as.matrix(coef(lasso_all_final))           # includes intercept
coefs_no_int <- coefs_mat[setdiff(rownames(coefs_mat), "(Intercept)"), , drop = FALSE]
coefs_no_int <- data.frame(coefs_no_int)
selected <- coefs_no_int %>% 
  filter(s0 != 0)
# write_csv(selected, file = "../Uwes_wishes/data/important_features_sorted_lasso_21.csv")

# ============================================
# 1) Lasso on PRM Peptide consolidation rule
# ============================================
# load the PRM data

prm_train <- read.csv("../Uwes_wishes/data/PRM_semi_pure_lassofeatures80iger_split.csv",
                      sep = ";", dec = ",")
rownames(prm_train) <- prm_train$Experimental_Group
prm_train$Experimental_Group <- NULL

prm_test <- read.csv("../Uwes_wishes/data/PRM_semi_pure_lassofeatures20iger_split.csv",
                     sep = ";", dec = ",")
rownames(prm_test) <- prm_test$Experimental_Group
prm_test$Experimental_Group <- NULL

# --- 2. Prepare matrices -----------------------------------------------------

X_train <- as.matrix(prm_train %>% select(-activ_vs_remission))
y_train <- as.numeric(prm_train$activ_vs_remission)

X_test  <- as.matrix(prm_test  %>% select(-activ_vs_remission))
y_test  <- as.numeric(prm_test$activ_vs_remission)

storage.mode(X_train) <- "double"
storage.mode(X_test)  <- "double"

# --- 3. Lasso (binomial) with cross-validation -------------------------------

set.seed(123)
cv_fit <- cv.glmnet(x = X_train, y = y_train, family = "binomial", alpha = 1)

best_lambda <- cv_fit$lambda.min
message("Best lambda (min): ", round(best_lambda, 6))

# --- 4. Final model and feature importance -----------------------------------

lasso_model <- glmnet(
  x      = X_train,
  y      = y_train,
  family = "binomial",
  alpha  = 1,
  lambda = best_lambda
)

coef_mat           <- as.matrix(coef(lasso_model, s = best_lambda))
important_features <- coef_mat[coef_mat[, 1] != 0, , drop = FALSE]

message("Selected features (incl. intercept): ", nrow(important_features))
print(important_features)

# Save importance table (sorted by absolute coefficient)
importance_df <- data.frame(
  feature     = rownames(important_features),
  coefficient = important_features[, 1]
) %>%
  arrange(desc(abs(coefficient)))


# write.csv(importance_df,
#           file      = "lasso_PRM_important_features.csv",
#           row.names = FALSE)


library(Boruta)
set.seed(1234)  # reproducible RF/Boruta
# --- load split ---
prm_train <- read.csv("../Uwes_wishes/data/PRM_semi_pure_lassofeatures80iger_split.csv", sep=";", dec=",")
rownames(prm_train) <- prm_train$Experimental_Group
prm_train$Experimental_Group <- NULL

prm_test <- read.csv("../Uwes_wishes/data/PRM_semi_pure_lassofeatures20iger_split.csv", sep=";", dec=",")
rownames(prm_test) <- prm_test$Experimental_Group
prm_test$Experimental_Group <- NULL
# 1) Make sure group is factor and create explicit train/test split if not already done
# (If you already have stable train/test files, confirm that they were created without peeking test.)


train_df$group <- factor(train_df$group, levels = c("0","1"))
test_df$group  <- factor(test_df$group,  levels = c("0","1"))

h <- row.names(train_df)
hi <- row.names(test_df)

# 2) Remove ID/rownames from columns if present
if("Experimental_Group" %in% colnames(train_df)) train_df$Experimental_Group <- NULL
if("Experimental_Group" %in% colnames(test_df))  test_df$Experimental_Group  <- NULL
selected_features <- c("MCAM_047" ,   "AHSG_004" ,    "H6PD_031", "LRG1_042" ,"TNXB_101" , "CDH2_014",
                       "MRC1_073", "KRT78_039", "CLEC3B_020",  "C9_076", "COMP_023", "ITIH4_103",
                       "F9_026" , "B2M_011", "IGHV1.24_036","PKM_089" ,"KRT6B_037", "group" )

# Optained  from the training cohort (e.g. using CV within training, stability, etc.)

# 3) Restrict predictors to the selected_features but ensure 'group' not included in predictors
selected_features_clean <- setdiff(selected_features, "group")
selected_features_clean <- intersect(selected_features_clean, colnames(train_df))
train_small <- train_df[, c(selected_features_clean, "group")]
test_small  <- test_df[,  c(selected_features_clean, "group")]

# 4) Preprocessing: fit on train, apply to both
library(caret)
pp <- preProcess(train_small[, selected_features_clean], method = c("center","scale"), verbose = FALSE)
train_pp <- predict(pp, train_small[, selected_features_clean])
test_pp  <- predict(pp, test_small[, selected_features_clean])
train_pp$group <- train_small$group
test_pp$group  <- test_small$group

# 5) Run Boruta only on the TRAINING data
set.seed(1234)  # Boruta reproducible
boruta_out <- Boruta(group ~ ., data = train_pp, doTrace = 2, maxRuns = 100)
print(boruta_out)
final_vars <- getSelectedAttributes(boruta_out, withTentative = T)

# print(boruta_output$finalDecision)
# AHSG_004      B2M_010       C9_076   CDH2_014     CLEC3B_020     COMP_023       F9_026     H6PD_031 IGHV1.24_036 
# Confirmed    Tentative    Tentative     Rejected    Confirmed    Confirmed    Confirmed    Rejected     Rejected 
# ITIH4_103    KRT78_038    LRG1_042     MCAM_048     MRC1_073     TNXB_101  PKM_089     KRT6B_037 
# Rejected    Tentative    Confirmed    Confirmed    Confirmed    Rejected    Rejected    Rejected  


