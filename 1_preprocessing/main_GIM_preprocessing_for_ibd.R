# preprocess the data for the analysis of intubation-at-blood-draw(ibd)

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

# preprocess and export/save the data
library(readxl)

# get preprocessed datasets
Ds <-c("metabo", "proteo") %>% {structure(., names=.)} %>% 
  lapply(function(dataset){
    
    file_path = paths[[dataset]]
    metadata_path = paths[["clinical_metadata"]]
    
    #### loading & preprocessing ----
    # load data
    D <- 
      GIM_data_loader(
        dataset=dataset,
        file_path = file_path,
        metadata_path = metadata_path,
        exclude.intubated=F) %>%
      # preprocess
      GIM_preprocessing()
    
    D$int_at_bld_drw = c("No","Yes")[D$already_intubated + 1]
    D
  })

save(file = "data_preprocessed_for_ibd.Rdata", Ds)
 

# # check if pre processed data is identical to previous version
# Ds_new = Ds
# load("~/Documents/R/_github_/GIM2/1_preprocessing/data_intubation_at_blood_draw_GIM_final.Rdata")
# 
# # proteomics
# sum( assay( Ds$dataset1_proteo ) - assay( Ds_new$dataset1_proteo ) )
# table( Ds_new$dataset1_proteo$int_at_bld_drw, Ds$dataset1_proteo$int_at_bld_drw, useNA = "ifany")
# # metabolomics
# sum( assay( Ds$dataset123 ) - assay( Ds_new$dataset123 ) )
# table( Ds_new$dataset123$int_at_bld_drw, Ds$dataset123$int_at_bld_drw, useNA = "ifany")
