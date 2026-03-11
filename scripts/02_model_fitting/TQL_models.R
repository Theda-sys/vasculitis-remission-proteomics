# ============================================================
# SAVE ALL THREE MODELS — Berlin, Train60, Train60_Concentration
# Run this script ONCE to freeze all models
# ============================================================
library(pROC); library(dplyr); library(readr); library(jsonlite)
set.seed(7)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
data_csv <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
outdir   <- "../proteomics_ml/data/models_and_artifacts/final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
dfz <- read_csv(data_csv)

if (!"cohort" %in% colnames(dfz)) stop("column 'cohort' not found in dfz")

# ----------------------------------------------------------
# DEFINE PREDICTORS
# ----------------------------------------------------------
# 7-protein panel
prots_7 <- c("AHSG_001","CLEC3B_020","COMP_023",
             "F9_027","LRG1_042","MCAM_047","MRC1_073..")


# Check all predictors exist
check_cols <- function(df, vars) {
  miss <- setdiff(vars, colnames(df))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
}
check_cols(dfz, prots_7)

# ----------------------------------------------------------
# FORMULAS
# ----------------------------------------------------------
formula_7    <- as.formula(paste("group01 ~", paste(prots_7,    collapse = " + ")))

# ----------------------------------------------------------
# DEFINE SPLITS
# ----------------------------------------------------------
data_trainB <- dfz %>% filter(cohort == "Berlin")

train60_file <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test40_file  <- "../Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"

if (file.exists(train60_file) && file.exists(test40_file)) {
  load(train60_file)   # provides full_meta_train
  load(test40_file)    # provides full_meta_test
  train60 <- dfz %>% filter(Patient %in% full_meta_train$Patient)
  test40  <- dfz %>% filter(Patient %in% full_meta_test$Patient)
} else {
  stop("Train/test patient list files not found")
}

cat("Berlin n=",   nrow(data_trainB), "\n")
cat("Train60 n=",  nrow(train60),     "\n")
cat("Test40 n=",   nrow(test40),      "\n")

# ----------------------------------------------------------
# FIT ALL THREE MODELS
# ----------------------------------------------------------

# 1) Berlin 7-protein (frozen external validation model)
model_berlin_7      <- glm(formula_7,    data = data_trainB, family = binomial())
if (!isTRUE(model_berlin_7$converged))      warning("Berlin 7-panel did not converge.")

# 2) Train60 7-protein (60/40 split)
model_train60_7     <- glm(formula_7,    data = train60,     family = binomial())
if (!isTRUE(model_train60_7$converged))     warning("Train60 7-panel did not converge.")

# ----------------------------------------------------------
# SAVE MODEL RDS OBJECTS
# ----------------------------------------------------------
saveRDS(model_berlin_7,     file.path(outdir, "model_7panel_Berlin_glm.rds"))
saveRDS(model_train60_7,    file.path(outdir, "model_7panel_Train60_glm.rds"))

# ----------------------------------------------------------
# SAVE COEFFICIENT TABLES (CSV)
# ----------------------------------------------------------
save_coef_csv <- function(model, filename) {
  tab <- as.data.frame(coef(summary(model)))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  write_csv(tab, file.path(outdir, filename))
}

save_coef_csv(model_berlin_7,     "coef_7panel_Berlin.csv")
save_coef_csv(model_train60_7,    "coef_7panel_Train60.csv")

# ----------------------------------------------------------
# SAVE COEFFICIENTS AS JSON
# ----------------------------------------------------------
write_json(as.list(coef(model_berlin_7)),
           file.path(outdir, "coef_7panel_Berlin.json"),     pretty = TRUE, auto_unbox = TRUE)
write_json(as.list(coef(model_train60_7)),
           file.path(outdir, "coef_7panel_Train60.json"),    pretty = TRUE, auto_unbox = TRUE)
# ----------------------------------------------------------
# FREEZE PREDICTED SCORES AT POINT OF SAVING
# ----------------------------------------------------------
score <- function(df, model) as.numeric(predict(model, newdata = df, type = "response"))

# test40 — all three models
test40$pred_train60_7     <- score(test40, model_train60_7)


write_csv(
  test40 %>% select(Patient, cohort, group01, pred_train60_7),
  file.path(outdir, "scores_test40_frozen.csv")
)

# Prague (data_testP) — Berlin model only (external validation)
data_testP <- dfz %>% filter(cohort == "Prague")
data_testP$pred_berlin_7 <- score(data_testP, model_berlin_7)

write_csv(
  data_testP %>% select(Patient, cohort, group01, pred_berlin_7),
  file.path(outdir, "scores_Prague_frozen.csv")
)

# Berlin model for seperation
data_testB <- dfz %>% filter(cohort == "Berlin")
data_testB$pred_berlin_7 <- score(data_testB, model_berlin_7)
write_csv(
  data_testB %>% select(Patient, cohort, group01, pred_berlin_7),
  file.path(outdir, "scores_data_testB_frozen.csv")
) 

# 60 model for seperation
# test40 — all three models
train60$pred_train60_7     <- score(train60, model_train60_7)
write_csv(
  train60 %>% select(Patient, cohort, group01, pred_train60_7),
  file.path(outdir, "scores_train60_frozen.csv")
)

