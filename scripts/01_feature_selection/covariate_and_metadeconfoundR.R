# ============================================================
# covaariate testing and metadeconfoundR
# ============================================================
# metadeconfoundR run was done on the global proteom in the 
# discovery cohort
# ============================================================

library(metadeconfoundR)  
library(dplyr)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(cowplot)  
library(ggpubr)    
library(scales)
library(svglite)

data.pca <- read_csv(file = "../26_RevisionNature/data605_datafreez.csv")
data.pca <- as.data.frame(data.pca)
row.names(data.pca) <- data.pca$Row.names.y


data.pca <- data.pca[!grepl("^HC", data.pca$Row.names), ]
data.pca$diseasestatus <- NA
data.pca$diseasestatus[grepl("active", data.pca$Row.names.y, ignore.case = TRUE)] <- "active"
data.pca$diseasestatus[grepl("remission", data.pca$Row.names.y, ignore.case = TRUE)] <- "remission"
#active = 0, remission = 1
data.pca$disease_numeric <- ifelse(data.pca$diseasestatus == "active", 0, 1)
data.pca$disease_numeric <- as.numeric(data.pca$disease_numeric)

data.pca <- data.pca %>% 
  filter(cohort == "Berlin")


feature <- data.pca[, c(1:605)]
meta <- data.pca[, c(606:725 )]

# choose the metadata columns you actually want
metaMat <- meta[, c(
  "disease_numeric",           # outcome: 0 = remission, 1 = active (must be first)
  "Age_years",
  "Sex_female_male",
  "Creatinine_mg_dL",
  "eGFR_CKD_EPI_2009", 
  "eGFR", # or your CKD-EPI 2021 column name
  "CRP_Plasma_mg_L",
  "Hemoglobin_g_dL",
  "Hematokrit_per",
  "Leukocytes",
  "Neutrophils",
  "Platelets",
  "ANCA_titer_PRAG",
  "ANCA_ELISA_BERLIN",
  "Kidney_involved_YES1",
  "Lung_involved_YES1",
  "ENT_involved_YES1",
  "Muscle_Joints_involved_YES1",
  "Skin_Mouth_Eye_involved_YES1",
  "CNS_involved_YES1"
)]

# make sure first column is 0/1 case status

# recode sex and cohort etc. to numeric 0/1 or factors as needed
metaMat$Sex_female_male <- ifelse(metaMat$Sex_female_male == "female", 1, 0)

# metaMat$cohort <- ifelse(metaMat$cohort == "Berlin", 0, 1)  # or use model.matrix later

# ensure all binary YES1 columns are 0/1 numeric
bin_cols <- grep("YES1$", names(metaMat), value = TRUE)
metaMat[bin_cols] <- lapply(metaMat[bin_cols], function(x) as.numeric(x))

# reduce the dataset
ad.index.keep <- which(colSums(feature)*100/(sum(feature)) > 0.01)
feature <- feature[, ad.index.keep]
feature <- feature[order(rownames(feature)), ]
metaMat <- metaMat[order(rownames(metaMat)), ]


library(metadeconfoundR)

set.seed(1)
mc_out <- MetaDeconfound(
  featureMat = feature,
  metaMat    = metaMat,
  nnodes     = 4,          # or number of cores you have
  QCutoff    = 0.1,
  DCutoff    = 0,          # or a small effect-size cutoff if you prefer
  logLevel   = "ERROR",
  returnLong = TRUE
)

# write_csv(mc_out, file = "../MoreSup/metdeconfoundR_ouput/metadeconfrounderR.csv")

mc_out <- read_csv(file = "../MoreSup/metdeconfoundR_ouput/metadeconfrounderR.csv")
## 1. Prepare data: cap Ds to [-1, 1] for colour scale
out_disease <- mc_out %>%
  filter(
    metaVariable == "disease_numeric",   # only disease
    status != "NS",                      # remove NS
    !str_detect(status, "^C:")           # remove all C:... statuses
  ) # 134 disease status associated 

mc_out_disease <- mc_out[mc_out$feature %in% out_disease$feature, ]
key_covariates <- c("CRP_Plasma_mg_L", "Hemoglobin_g_dL",
                    "Hematokrit_per","Platelets", "disease_numeric","eGFR")

mc_out_disease <- mc_out_disease[mc_out_disease$metaVariable %in% key_covariates, ]
nice_names <- c(
  "CRP_Plasma_mg_L"    = "C-reactive protein (mg/L)",
  "eGFR"               = "eGFR",
  "Hemoglobin_g_dL"    = "hemoglobin (g/dL)",
  "Hematokrit_per"     = "hematocrit (%)",
  "Platelets"          = "platelet count",
  "disease_numeric"    = "remission vs active"
)

mc_plot <- mc_out_disease %>%
  mutate(
    metaVariable = factor(nice_names[metaVariable], levels = nice_names),
    signif_lab = case_when(
      Qs < 0.001 ~ "***",
      Qs < 0.01  ~ "**",
      Qs < 0.1   ~ "*",
      TRUE       ~ ""
    ),
    has_C = str_detect(status, "C"),
    star_colour = ifelse(has_C, "grey40", "black")
  )

p <- ggplot(mc_plot, aes(x = metaVariable, y = feature, fill = Ds)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Effect size (D)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank())

p + geom_text(
  data = mc_plot %>% filter(signif_lab != ""),
  aes(label = signif_lab, colour = star_colour),
  size = 3
) +
  scale_colour_identity() +
  labs(x = "Metadata", y = "Feature") -> p

# ggsave(p, file = "../MoreSup/metdeconfoundR_ouput/heatmaps/MetadeconfRheatmap134.svg", width = 4.5, height = 28)


lasso21 <- read.table(file = "../Uwes_wishes/Nature_Code/data/21_lasso_para_withgene_names.csv", sep = ",", header =T)
lasso21_feat <- lasso21$Names_from_dataframe

mc_out_21 <- mc_out[mc_out$feature %in% lasso21_feat, ]

key_covariates <- c("CRP_Plasma_mg_L", "Hemoglobin_g_dL",
                    "Hematokrit_per","Platelets", "disease_numeric","eGFR")

mc_out_21 <- mc_out_21[mc_out_21$metaVariable %in% key_covariates, ]
nice_names <- c(
  "CRP_Plasma_mg_L"    = "C-reactive protein (mg/L)",
  "eGFR"               = "eGFR",
  "Hemoglobin_g_dL"    = "hemoglobin (g/dL)",
  "Hematokrit_per"     = "hematocrit (%)",
  "Platelets"          = "platelet count",
  "disease_numeric"    = "remission vs active"
)

mc_plot <- mc_out_21 %>%
  mutate(
    metaVariable = factor(nice_names[metaVariable], levels = nice_names),
    signif_lab = case_when(
      Qs < 0.001 ~ "***",
      Qs < 0.01  ~ "**",
      Qs < 0.1   ~ "*",
      TRUE       ~ ""
    ),
    has_C = str_detect(status, "C"),
    star_colour = ifelse(has_C, "grey40", "black")
  )

p <- ggplot(mc_plot, aes(x = metaVariable, y = feature, fill = Ds)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Effect size (D)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank())

p + geom_text(
  data = mc_plot %>% filter(signif_lab != ""),
  aes(label = signif_lab, colour = star_colour),
  size = 3
) +
  scale_colour_identity() +
  labs(x = "Metadata", y = "Feature") -> p

# ggsave(p, file = "../MoreSup/metdeconfoundR_ouput/heatmaps/MetadeconfRheatmap21.svg", width = 4.5, height = 7)

panel_proteins <- c("AHSG","CLEC3B","COMP",
                    "F9","LRG1","MCAM","MRC1")
mc_out_7 <- mc_out[mc_out$feature %in% panel_proteins, ]

key_covariates <- c( "CRP_Plasma_mg_L", 
                     "Hemoglobin_g_dL","Hematokrit_per","Platelets", "disease_numeric","eGFR")

mc_out_7 <- mc_out_7[mc_out_7$metaVariable %in% key_covariates, ]
nice_names <- c(
  "CRP_Plasma_mg_L"    = "C-reactive protein (mg/L)",
  "eGFR"               = "eGFR",
  "Hemoglobin_g_dL"    = "hemoglobin (g/dL)",
  "Hematokrit_per"     = "hematocrit (%)",
  "Platelets"          = "platelet count",
  "disease_numeric"    = "remission vs active"
)

# write_csv(mc_out_7, file = "../MoreSup/MetadeconfRoutput7.csv")
mc_plot <- mc_out_7 %>%
  mutate(
    metaVariable = factor(nice_names[metaVariable], levels = nice_names),
    signif_lab = case_when(
      Qs < 0.001 ~ "***",
      Qs < 0.01  ~ "**",
      Qs < 0.1   ~ "*",
      TRUE       ~ ""
    ),
    has_C = str_detect(status, "C"),
    star_colour = ifelse(has_C, "grey40", "black")
  )

p <- ggplot(mc_plot, aes(x = metaVariable, y = feature, fill = Ds)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Effect size (D)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank())

p + geom_text(
  data = mc_plot %>% filter(signif_lab != ""),
  aes(label = signif_lab, colour = star_colour),
  size = 3
) +
  scale_colour_identity() +
  labs(x = "Metadata", y = "Feature") -> p

# ggsave(p, file = "../MoreSup/metdeconfoundR_ouput/heatmaps/MetadeconfRheatmap7.svg", width = 3.5, height = 3)


# Spearman correlations between candidate proteins and CRP, eGFR, 
# hematocrit, haemoglobin, platelet count, and ANCA to evaluate potential redundancy (done post-hoc)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
artifact_dir <- "../proteomics_ml/data/models_and_artifacts/final/"
data_csv     <- "../proteomics_ml/data/TQLData_combinedCohorts.csv"
outdir       <- "../proteomics_ml/data/covariate_assignment"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prots_7 <- c("AHSG_001", "CLEC3B_020", "COMP_023",
             "F9_027", "LRG1_042", "MCAM_047", "MRC1_073..")

# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
dfzt        <- read_csv(data_csv)
data_trainB <- dfzt %>% filter(cohort == "Berlin")
data_testP  <- dfzt %>% filter(cohort == "Prague")

message("Berlin n=", nrow(data_trainB),
        "; Prague n=", nrow(data_testP),
        " (events in Prague = ", sum(data_testP$group01 == 1, na.rm = TRUE), ")")

key_covariates <- c("CRP_Plasma_mg_L", "Hemoglobin_g_dL",
                    "Hematokrit_per","Platelets",
                    "zANCA", "eGFR_CKD_EPI_2009")

proteins_7 <- c("AHSG_001","CLEC3B_020","COMP_023",
                "F9_027","LRG1_042","MCAM_047","MRC1_073..")

# Sanity: make sure variables exist and numeric in Berlin
stopifnot(all(proteins_7 %in% colnames(data_trainB)))
stopifnot(all(key_covariates %in% colnames(data_trainB)))

is_num <- sapply(data_trainB[, c(proteins_7, key_covariates)], is.numeric)
if (!all(is_num)) {
  warning("Some proteins or covariates are not numeric in Berlin; non-numeric pairs will be skipped.")
}

# Spearman correlations: 7 proteins x key covariates (Berlin)
all_corr <- expand.grid(
  protein   = proteins_7,
  covariate = key_covariates,
  stringsAsFactors = FALSE
)

all_corr$rho <- NA_real_
all_corr$p   <- NA_real_

for (i in seq_len(nrow(all_corr))) {
  p_name <- all_corr$protein[i]
  c_name <- all_corr$covariate[i]
  
  if (!is.numeric(data_trainB[[p_name]]) ||
      !is.numeric(data_trainB[[c_name]])) {
    next
  }
  
  sub <- data_trainB[, c(p_name, c_name)]
  sub <- sub[complete.cases(sub), ]
  if (nrow(sub) < 5) next
  
  ct <- suppressWarnings(cor.test(sub[[1]], sub[[2]],
                                  method = "spearman"))
  all_corr$rho[i] <- unname(ct$estimate)
  all_corr$p[i]   <- ct$p.value
}

library(dplyr)
all_corr <- all_corr %>%
  group_by(covariate) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

# write.csv(all_corr,
#           file.path(outdir, "Spearman_7panel_vs_keycovariates_Berlin.csv"),
#           row.names = FALSE)

num_cols <- sapply(all_corr, is.numeric)
all_corr[, num_cols] <- lapply(all_corr[, num_cols], round, 3)


# write.csv(all_corr,
#           file.path(outdir, "Spearman_7panel_vs_keycovariates_Berlin3Digits.csv"),
#           row.names = FALSE)

