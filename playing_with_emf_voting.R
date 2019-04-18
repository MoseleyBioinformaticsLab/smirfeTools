# helping me figure out labeled EMF voting
library(smirfeTools)
library(furrr)
plan(multiprocess)
set_internal_map(furrr::future_map)
#assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")
assigned_data = readRDS("lung_matched_tissue_raw_smirfe_assignments_2018-09-07_imf_adduct_fusion1_cancer.rds")

extracted_data = extract_assigned_data(assigned_data)
