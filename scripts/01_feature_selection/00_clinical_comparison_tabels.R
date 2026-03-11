
# Statistic for clinical data 
# t-Test (or Wilcoxon) for numeric: Age_years
# CRP_Plasma_mg_L
# Creatinine_mg_dL
# Hemoglobin_g_dL
# Hematokrit_per
# Leukocytes
# Neutrophils
# Platelets
# BVAS

# Chi²-Test (or Fisher-Test) for categorical data: Sex_female_male
# Kidney_involved_YES1
# Lung_involved_YES1
# ENT_involved_YES1
# Muscle_Joints_involved_YES1
# Skin_Mouth_Eye_involved_YES1
# CNS_involved_YES1
# split Prague and Berlin 

library(dplyr)
library(tidyr)
library(gtsummary)
library(flextable)
library(officer)

final_meta <- read.table(file = "24828_final_tabels/250228_Data Freeze TQL(BER_PRG)_meta_clean_outliers_editUJ_new_ML.csv", 
                         sep = ";", 
                         header = TRUE, 
                         na.strings = c("NA", ""), 
                         dec = ",", row.names = 1)

final_meta$group <- factor(final_meta$ACTIVE_1Yes, levels = c(0,1), labels = c("remission", "active"))
final_meta$Skin_Mouth_Eye_involved_YES1 <- as.factor(final_meta$Skin_Mouth_Eye_involved_YES1)
bin_vars <- c("Kidney_involved_YES1", "Lung_involved_YES1",
              "ENT_involved_YES1", "Muscle_Joints_involved_YES1",
              "Skin_Mouth_Eye_involved_YES1", "CNS_involved_YES1")

for (var in bin_vars) {
  final_meta[[var]] <- as.character(final_meta[[var]])
  final_meta[[var]] <- dplyr::case_when(
    final_meta[[var]] %in% c("1", "yes") ~ 1,
    final_meta[[var]] %in% c("0", "no") ~ 0,
    TRUE ~ NA_real_
  )
}

final_meta$Sex_female_male <- factor(final_meta$Sex_female_male, levels = c("female", "male"))
final_meta$group <- factor(final_meta$ACTIVE_1Yes, levels = c(0,1), labels = c("remission", "active"))


# Berlin subset
final_meta_berlin <- final_meta %>% 
  filter(cohort == "Berlin") %>%
  select(-c(62:76))  # TQL remove


# Ensure AAV_group and group are set
final_meta_berlin$AAV_group <- factor(final_meta_berlin$AAV_group, levels = c("MPO", "PR3"))
final_meta_berlin$group <- factor(final_meta_berlin$ACTIVE_1Yes, levels = c(0,1), labels = c("remission", "active"))


# Binary variables cleanup (if not already done)
bin_vars <- c("Kidney_involved_YES1", "Lung_involved_YES1",
              "ENT_involved_YES1", "Muscle_Joints_involved_YES1",
              "Skin_Mouth_Eye_involved_YES1", "CNS_involved_YES1")

for (var in bin_vars) {
  final_meta_berlin[[var]] <- as.character(final_meta_berlin[[var]])
  final_meta_berlin[[var]] <- dplyr::case_when(
    final_meta_berlin[[var]] %in% c("1", "yes") ~ 1,
    final_meta_berlin[[var]] %in% c("0", "no") ~ 0,
    TRUE ~ NA_real_
  )
}

# Make sure Sex is properly formatted
final_meta_berlin$Sex_female_male <- factor(final_meta_berlin$Sex_female_male, levels = c("female", "male"))

# Full list of variables
clinical_vars_berlin_group <- c("group", "Age_years", "CRP_Plasma_mg_L", "Creatinine_mg_dL", "Hemoglobin_g_dL",
                                "Hematokrit_per", "Leukocytes", "Neutrophils", "Platelets", "BVAS",
                                "Sex_female_male", "Kidney_involved_YES1", "Lung_involved_YES1",
                                "ENT_involved_YES1", "Muscle_Joints_involved_YES1",
                                "Skin_Mouth_Eye_involved_YES1", "CNS_involved_YES1")

# Drop group if comparing by AAV_group
clinical_vars_berlin_aav <- clinical_vars_berlin_group[clinical_vars_berlin_group != "group"]

create_safe_comparison_table <- function(data, compare_var, filename, variable_list) {
  # Remove compare_var from variables (if present)
  variable_list <- setdiff(variable_list, compare_var)
  
  # Ensure grouping variable is factor
  data[[compare_var]] <- as.factor(data[[compare_var]])
  
  # Remove rows with NA in grouping variable
  data <- data %>% filter(!is.na(.data[[compare_var]]))
  
  # Identify valid variables with at least some variability across groups
  valid_vars <- variable_list[
    sapply(variable_list, function(var) {
      temp <- data %>% select(all_of(var), !!sym(compare_var))
      tryCatch({
        n_lvls <- temp %>%
          group_by(!!sym(compare_var)) %>%
          summarise(n = n_distinct(.data[[var]], na.rm = TRUE)) %>%
          pull(n)
        all(n_lvls >= 1) && sum(n_lvls >= 2) >= 1
      }, error = function(e) FALSE)
    })
  ]
  
  # Warn user if variables were dropped
  dropped <- setdiff(variable_list, valid_vars)
  if (length(dropped) > 0) {
    message("⚠️ Dropped variables due to lack of variability in groups: ", paste(dropped, collapse = ", "))
  }
  
  # Generate dynamic label list
  label_map <- list(
    Age_years = "Age — year",
    CRP_Plasma_mg_L = "CRP [mg/l]",
    Creatinine_mg_dL = "Creatinine [mg/dl]",
    Hemoglobin_g_dL = "Hemoglobin [g/dl]",
    Hematokrit_per = "Hematokrit [%]",
    Leukocytes = "Leukocytes [/nl]",
    Neutrophils = "Neutrophils [/nl]",
    Platelets = "Platelets [/nl]",
    BVAS = "BVAS (0–63)",
    Sex_female_male = "Sex",
    Kidney_involved_YES1 = "Kidney",
    Lung_involved_YES1 = "Lung",
    ENT_involved_YES1 = "Ear/Nose/Throat",
    Muscle_Joints_involved_YES1 = "Muscle/Joints",
    Skin_Mouth_Eye_involved_YES1 = "Skin/Mouth/Eyes",
    CNS_involved_YES1 = "Central nervous system"
  )
  
  label_list <- label_map[names(label_map) %in% valid_vars]
  
  # Create and save table
  tbl <- data %>%
    select(all_of(valid_vars), !!sym(compare_var)) %>%
    tbl_summary(
      by = !!sym(compare_var),
      label = label_list,
      statistic = list(
        all_continuous() ~ "{mean} ± {sd}",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 2,
      missing = "no"
    ) %>%
    add_p(
      test = list(
        all_continuous() ~ "t.test",
        all_categorical() ~ "chisq.test"
      ),
      pvalue_fun = ~style_pvalue(.x, digits = 3)
    ) %>%
    add_overall() %>%
    bold_labels()
  
  save_as_docx(as_flex_table(tbl), path = filename)
}




data1 <- final_meta_berlin %>% filter(group == "active")
create_comparison_table(
  data = data1,
  compare_var = "AAV_group",
  filename = "Revision/Prep_dataframes/BER_Active_MPO_vs_PR3.docx",
  variable_list = clinical_vars_berlin_aav
)


# 1. Berlin: Remission MPO vs PR3
data_rem <- final_meta_berlin %>% filter(group == "remission")
create_safe_comparison_table(data_rem, "AAV_group", "Revision/Prep_dataframes/BER_Remission_MPO_vs_PR3.docx", clinical_vars_berlin_aav)

# 2. Berlin: PR3 Active vs Remission
data_pr3 <- final_meta_berlin %>% filter(AAV_group == "PR3")
create_safe_comparison_table(data_pr3, "group", "Revision/Prep_dataframes/BER_PR3_Active_vs_Remission.docx", clinical_vars_berlin_group)

# 3. Berlin: MPO Active vs Remission
data_mpo <- final_meta_berlin %>% filter(AAV_group == "MPO")
create_safe_comparison_table(data_mpo, "group", "Revision/Prep_dataframes/BER_MPO_Active_vs_Remission.docx", clinical_vars_berlin_group)



final_meta_prague <- final_meta %>%
  filter(cohort == "Prague") %>%
  select(-c(62:76))  # remove TQL etc.

final_meta_prague$group <- factor(final_meta_prague$ACTIVE_1Yes, levels = c(0, 1), labels = c("remission", "active"))
final_meta_prague$AAV_group <- factor(final_meta_prague$AAV_group, levels = c("MPO", "PR3"))

clinical_vars_prague_group <- c("group", "Age_years", "CRP_Plasma_mg_L", "Creatinine_mg_dL", "Hemoglobin_g_dL",
                                "Hematokrit_per", "Leukocytes", "Neutrophils", "Platelets", "BVAS",
                                "Sex_female_male", "Kidney_involved_YES1", "Lung_involved_YES1",
                                "ENT_involved_YES1", "Muscle_Joints_involved_YES1",
                                "Skin_Mouth_Eye_involved_YES1", "CNS_involved_YES1")

# remove 'group' if using AAV_group as comparison variable
clinical_vars_prague_aav <- clinical_vars_prague_group[clinical_vars_prague_group != "group"]

data_active_prague <- final_meta_prague %>% filter(group == "active")

create_safe_comparison_table(
  data = data_active_prague,
  compare_var = "AAV_group",
  filename = "Revision/Prep_dataframes/PRAGUE_Active_MPO_vs_PR3.docx",
  variable_list = clinical_vars_prague_aav
)

data_remission_prague <- final_meta_prague %>% filter(group == "remission")

create_safe_comparison_table(
  data = data_remission_prague,
  compare_var = "AAV_group",
  filename = "Revision/Prep_dataframes/PRAGUE_Remission_MPO_vs_PR3.docx",
  variable_list = clinical_vars_prague_aav
)


data_prague_pr3 <- final_meta_prague %>% filter(AAV_group == "PR3")

create_safe_comparison_table(
  data = data_prague_pr3,
  compare_var = "group",
  filename = "Revision/Prep_dataframes/PRAGUE_PR3_Active_vs_Remission.docx",
  variable_list = clinical_vars_prague_group
)

data_prague_mpo <- final_meta_prague %>% filter(AAV_group == "MPO")

create_safe_comparison_table(
  data = data_prague_pr3,
  compare_var = "group",
  filename = "Revision/Prep_dataframes/PRAGUE_MPO_Active_vs_Remission.docx",
  variable_list = clinical_vars_prague_group
)


## Missing ANCA titer and EryUrine 

# PRAGUE
final_meta_prague$Urin_ERY_POS_Bhigher0.Phigher20_YES1 <- ifelse(
  final_meta_prague$Urin_ERY_POS_Bhigher0.Phigher20_YES1 == "?", NA,
  final_meta_prague$Urin_ERY_POS_Bhigher0.Phigher20_YES1
)
final_meta_prague$Urin_ERY_POS_Bhigher0.Phigher20_YES1 <- as.numeric(final_meta_prague$Urin_ERY_POS_Bhigher0.Phigher20_YES1)
final_meta_prague$Urin_Erythrocytes_PRAG <- as.numeric(final_meta_prague$Urin_Erythrocytes_PRAG)

# BERLIN
final_meta_berlin$Urin_ERY_POS_Bhigher0.Phigher20_YES1 <- ifelse(
  final_meta_berlin$Urin_ERY_POS_Bhigher0.Phigher20_YES1 == "?", NA,
  final_meta_berlin$Urin_ERY_POS_Bhigher0.Phigher20_YES1
)
final_meta_berlin$Urin_ERY_POS_Bhigher0.Phigher20_YES1 <- as.numeric(final_meta_berlin$Urin_ERY_POS_Bhigher0.Phigher20_YES1)
final_meta_berlin$Urin_ERY_BERLIN <- as.numeric(final_meta_berlin$Urin_ERY_BERLIN)


generate_anca_urine_numeric_summary <- function(data, compare_var, filename,
                                                anca_var,
                                                urine1_var,
                                                urine2_var,
                                                label_anca,
                                                label_urine1,
                                                label_urine2) {
  compare_var <- rlang::ensym(compare_var)
  
  label_list <- list()
  label_list[[anca_var]] <- label_anca
  label_list[[urine1_var]] <- label_urine1
  label_list[[urine2_var]] <- label_urine2
  
  var_list <- c(anca_var, urine1_var, urine2_var)
  
  # Filter rows with compare_var present
  data <- data %>% filter(!is.na(.data[[rlang::as_string(compare_var)]]))
  
  # Variablen mit ausreichender Varianz
  valid_vars <- var_list[sapply(var_list, function(var) {
    temp <- data %>% select(all_of(var), !!compare_var)
    tryCatch({
      n_lvls <- temp %>%
        group_by(!!compare_var) %>%
        summarise(n = n_distinct(.data[[var]], na.rm = TRUE)) %>%
        pull(n)
      all(n_lvls >= 1) && sum(n_lvls >= 2) >= 1
    }, error = function(e) FALSE)
  })]
  
  # Tabelle
  tbl <- data %>%
    select(all_of(valid_vars), !!compare_var) %>%
    tbl_summary(
      by = !!compare_var,
      label = label_list[names(label_list) %in% valid_vars],
      statistic = list(all_continuous() ~ "{mean} ± {sd}"),
      digits = all_continuous() ~ 2,
      missing = "no"
    ) %>%
    add_p(
      test = all_continuous() ~ "t.test",
      pvalue_fun = ~style_pvalue(.x, digits = 3)
    ) %>%
    add_overall() %>%
    bold_labels()
  
  save_as_docx(as_flex_table(tbl), path = filename)
}

generate_anca_urine_numeric_summary(
  data = final_meta_berlin,
  compare_var = group,
  filename = "Revision/Prep_dataframes/BER_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_ELISA_BERLIN",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_ERY_BERLIN",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_berlin %>% filter(AAV_group == "PR3"),
  compare_var = group,
  filename = "Revision/Prep_dataframes/BER_PR3_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_ELISA_BERLIN",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_ERY_BERLIN",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_berlin %>% filter(AAV_group == "MPO"),
  compare_var = group,
  filename = "Revision/Prep_dataframes/BER_MPO_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_ELISA_BERLIN",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_ERY_BERLIN",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_berlin %>% filter(group == "active"),
  compare_var = AAV_group,
  filename = "Revision/Prep_dataframes/BER_active_MPO_vs_PR3_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_ELISA_BERLIN",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_ERY_BERLIN",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_berlin %>% filter(group == "remission"),
  compare_var = AAV_group,
  filename = "Revision/Prep_dataframes/BER_remission_MPO_vs_PR3_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_ELISA_BERLIN",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_ERY_BERLIN",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)


generate_anca_urine_numeric_summary(
  data = final_meta_prague,
  compare_var = group,
  filename = "Revision/Prep_dataframes/PRAGUE_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_titer_PRAG",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_Erythrocytes_PRAG",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_prague %>% filter(AAV_group == "PR3"),
  compare_var = group,
  filename = "Revision/Prep_dataframes/PRAGUE_PR3_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_titer_PRAG",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_Erythrocytes_PRAG",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_prague %>% filter(AAV_group == "MPO"),
  compare_var = group,
  filename = "Revision/Prep_dataframes/PRAGUE_MPO_active_vs_remission_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_titer_PRAG",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_Erythrocytes_PRAG",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_prague %>% filter(group == "active"),
  compare_var = AAV_group,
  filename = "Revision/Prep_dataframes/PRAGUE_active_MPO_vs_PR3_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_titer_PRAG",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_Erythrocytes_PRAG",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)

generate_anca_urine_numeric_summary(
  data = final_meta_prague %>% filter(group == "remission"),
  compare_var = AAV_group,
  filename = "Revision/Prep_dataframes/PRAGUE_remission_MPO_vs_PR3_ANCA_2Urine_NUMERIC.docx",
  anca_var = "ANCA_titer_PRAG",
  urine1_var = "Urin_ERY_POS_Bhigher0.Phigher20_YES1",
  urine2_var = "Urin_Erythrocytes_PRAG",
  label_anca = "ANCA Titer",
  label_urine1 = "Hematuria (YES1)",
  label_urine2 = "Erythrocytes (numeric)"
)
