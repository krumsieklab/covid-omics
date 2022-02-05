# we need to supplement these data seperately 
# script to generate preprocessed data for GGMs 
# corrected for status and intubated samples are added for maximal power in inferring GGMs 
library(maplet)
library(tidyverse)
library(magrittr)
library(glue)
setwd(glue::glue("{root}/1_preprocessing/"))

# where files located in 
if(!exists("paths")){
  paths = list(
    clinical_metadata = "gim_clinical_metadata.xlsx", 
    metabo = "gim_metabolomics.xlsx",
    proteo = "gim_proteomics.xlsx"
  )
}

# should SE object be saved
if(!exists("save_SE")) save_SE = F
# should preprocessed data be exported as xlsx
if(!exists("save_XLSX")) save_XLSX = F

# data loader and pre processing functions 
source("utils_GIM_dataloader.R")
source("utils_GIM_preprocessing.R")

# confounder formula
formula_conf <- "sex + age + Status"
# variables for correction
outcomes <- list(data.frame(outcome="sex", outcomeType="binary", outcomeMode="character"),
                 data.frame(outcome="age", outcomeType="numeric", outcomeMode="numeric"),
                 data.frame(outcome="Status", outcomeType="binary", outcomeMode="character")
) %>% do.call(rbind,.)

# correcting for confounders 
datasets <- c("proteo", "metabo")
Ds <- list()
for (dataset in datasets){
  
  file_path = paths[[dataset]]
  metadata_path = paths[["clinical_metadata"]]
  
  # load and preprocess
  D <-
    # load data
    GIM_data_loader(dataset=dataset,
                    file_path = file_path,
                    metadata_path = metadata_path) %>%
    # preprocess
    GIM_preprocessing() %>%
    # convert outcomes into numeric / factors 
    GIM_outcome_type_conversion (outcomes)  %>%
    # correct for confouders and groups
    mt_pre_confounding_correction(as.formula(glue("~ {formula_conf}"))) %>%
    # for coding convenience (empty statement to terminate %>% pipe)
    {.}
  Ds[[dataset]] <- D
}
rm(D)
names(Ds) <- c('Proteins', 'Metabolites')

save(file = "data_preprocessed_for_GGM.Rdata", Ds)


# # check if pre processed data is identical to previous version
# Ds_new = Ds
# load("~/Documents/R/_github_/GIM2/1_preprocessing/data_GGM_GIM_final.Rdata")
# 
# # proteomics
# sum( assay( Ds$Proteins ) - assay( Ds_new$Proteins ) )
# # metabolomics
# sum( assay( Ds$Metabolites ) - assay( Ds_new$Metabolites ) )
#
