library(pROC)
library(dplyr)
library(readr)
library(ResourceSelection)
library(boot)
library(openxlsx)
outdir   <- "../proteomics_ml/data/models_and_artifacts/final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# For the Calibation Based Method
concentration_calibartion_curve <- read.csv2("../Uwes_wishes/Nature_Code/data/concentration_calibartion_curve.csv", row.names = 1)
concentration_calibartion_curve <- data.frame(t(concentration_calibartion_curve))

# --- 0) Config: paths & options ----------------------------------------------
data_csv         <- "~/Documents/MDC/ml_proteomics/Uwes_wishes/Nature_Code/data/260216 V3_2026 von 250228_DataFreeze.csv"
train_meta_path  <- "~/Documents/MDC/ml_proteomics/Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_train.r"
test_meta_path   <- "~/Documents/MDC/ml_proteomics/Uwes_wishes/Nature_Code/data/2549_tqlfull_meta_test.r"

set.seed(7)  # report this seed in Methods

# 0) Read & basic split (if you already have data_train / data_test you can skip this block)
df <- read.csv(data_csv, stringsAsFactors = FALSE, check.names = FALSE)
if ("Row.names" %in% colnames(df)) df <- df[!grepl("^HC", df$Row.names), , drop = FALSE]
df$diseasestatus <- NA_character_
df$diseasestatus[grepl("active", df$Row.names, ignore.case = TRUE)] <- "active"
df$diseasestatus[grepl("remission", df$Row.names, ignore.case = TRUE)] <- "remission"
df$group <- ifelse(df$diseasestatus == "remission", 1L, 0L)  # group: 1 = remission, 0 = active
# berlin / prague split (only if objects not already present)
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqlmeta_full.r")
fullzANCA <- meta_full[, c("Patient", "ANCA_z_combined")]
colnames(fullzANCA) <- c("Patient", "zANCA")
row.names(fullzANCA) <- NULL
fullzANCA <- data.frame(fullzANCA)
fullzANCA <- fullzANCA %>%
  mutate(zANCA = as.numeric(zANCA))
# write_csv(fullzANCA, file = "../helperfullzANCA.csv")

dfz <- fullzANCA %>% 
  merge(df, by = "Patient")

final_data_raw <- read.csv(file = "../Uwes_wishes/Nature_Code/data/variance_final_data_raw.csv")

rn_dfz   <- dfz$Row.names
rn_vsn   <- final_data_raw$cleaned

length(rn_dfz)
length(rn_vsn)

# In both
length(intersect(rn_dfz, rn_vsn))

# Present only in dfz$Row.names
setdiff(rn_dfz, rn_vsn)

# Present only in final_data_raw$cleaned
setdiff(rn_vsn, rn_dfz)

# Any duplicates?
rn_dfz[duplicated(rn_dfz)]
rn_vsn[duplicated(rn_vsn)]


# vector of IDs to remove (the setdiff you just saw)
ids_to_remove <- c("MPO_active_A20_36",
                   "PR3_active_A20_21",
                   "Prag_MPO_Active_4140_P01_C09",
                   "Prag_MPO_Active_5665_P02_D03",
                   "Prag_MPO_Remission_4376_P01_D01",
                   "Prag_MPO_Remission_5537_P02_A06",
                   "Prag_MPO_Remission_5711_P02_B05")

# remove from final_data_raw by cleaned
final_data_raw_filtered <- final_data_raw[!(final_data_raw$cleaned %in% ids_to_remove), ]
final_data_raw_filtered$X <- NULL
final_data_raw_filtered$Row.names <- NULL
final_data_raw_filtered$group <- NULL

dfvsn <- merge(final_data_raw_filtered,
               dfz,
               by.x = "cleaned",
               by.y = "Row.names")


dfvsn <- dfvsn %>% 
  filter(is.na(finale_outliers_250228_UJ_YES1))

# canonical numeric group: 0 = active, 1 = remission
dfvsn$group01 <- ifelse(tolower(as.character(dfvsn$diseasestatus)) == "remission", 1L,
                        ifelse(tolower(as.character(dfvsn$diseasestatus)) == "active", 0L, NA_integer_))
if (any(is.na(dfvsn$group01))) warning("Some rows have NA in group01; inspect dfz$diseasestatus unique values.")

dfvsn <- dfvsn[, -c(102:117)]
# final <- dfvsn
# saveRDS(final, file = "../proteomics_ml/data/variance_dfvsn_final.rds")

load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_trainn.r")
load(file = "../Uwes_wishes/Nature_Code/data/2549_tqldata_test.r")


# Create VSN‑based train/test using the old variance_final_data_raw
vsn_data_based_test  <- dfvsn[dfvsn$cleaned %in% row.names(data_test), ]
vsn_data_based_train <- dfvsn[dfvsn$cleaned %in% row.names(data_train), ]


# 2) Fit 7‑panel on VSN data
formula_model4 <- group01 ~ AHSG.1 + CLEC3B.20 + COMP.23. + F9.27 +
  LRG1.42 + MCAM.47 + MRC1.73

model_7panel_vsn <- glm(formula_model4,
                        data   = vsn_data_based_train,
                        family = binomial())

saveRDS(model_7panel_vsn, file.path(outdir, "model_conc_Train60_glm.rds"))

save_coef_csv <- function(model, filename) {
  tab <- as.data.frame(coef(summary(model)))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  write_csv(tab, file.path(outdir, filename))
}

save_coef_csv(model_7panel_vsn, "coef_conc_Train60.csv")

jsonlite::write_json(as.list(coef(model_7panel_vsn)),
           file.path(outdir, "coef_conc_Train60.json"),      pretty = TRUE, auto_unbox = TRUE)


score <- function(df, model) as.numeric(predict(model, newdata = df, type = "response"))

vsn_data_based_test$pred_train60_conc <- score(vsn_data_based_test, model_7panel_vsn)

write_csv(
  vsn_data_based_test %>% select(Patient, cohort, group01, pred_train60_conc),
  file.path(outdir, "scores_test40_frozen_concentration.csv")
)


vsn_data_based_train$pred_train60_conc <- score(vsn_data_based_train, model_7panel_vsn)
write_csv(
  vsn_data_based_train %>% select(Patient, cohort, group01, pred_train60_conc),
  file.path(outdir, "scores_train60_concentration_frozen.csv")
)
