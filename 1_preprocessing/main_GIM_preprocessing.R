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

# export the data ------------------------------------------------------------
# preprocess and export/save the data

library(readxl)
library(openxlsx)
f_export <- function(dataset, outcomes, exclude.intubated, 
                     file_path = paths[[dataset]], 
                     metadata_path = paths[["clinical_metadata"]],
                     datestamp, dir, save_SE = FALSE, save_XLSX = FALSE){
  corradd <- ""
  intubadd <- if (exclude.intubated)""else"_withintub"
  file_name =  paste0(dir, sprintf("/data_%s_%s%s%s.xlsx", datestamp, dataset, corradd, intubadd) )
  
  #### loading & preprocessing ----
  D <- # load data
    GIM_data_loader(dataset=dataset, 
                    file_path=file_path,
                    metadata_path= metadata_path,
                    exclude.intubated=exclude.intubated) %>%
    # preprocess
    GIM_preprocessing() %>%
    # convert outcomes into numeric / factors 
    GIM_outcome_type_conversion( outcomes )
  
  if(save_SE){
    save(file = gsub("xlsx", "Rdata", file_name), D)
  }
  
  rowData(D)$kegg_db <- unlist(lapply(rowData(D)$kegg_db, function(x) toString(unlist(x))))
  
  if(save_XLSX){
    # write the SE into excel
    wb = createWorkbook()
    sheet = addWorksheet(wb, "data")
    writeData(wb, sheet=sheet, D %>% assay, rowNames = T, colNames = T)
    # this row is to use the original metabolite names
    writeData(wb, sheet=sheet, rowData(D)$BIOCHEMICAL, startCol = 1, startRow = 2)
    # this is needed for MT
    writeData(wb, sheet=sheet, "BIOCHEMICAL", startCol = 1, startRow = 1)
    sheet = addWorksheet(wb, "metanno")
    writeData(wb, sheet=sheet, D %>% rowData %>% as.data.frame, colNames = T, rowNames=F)
    sheet = addWorksheet(wb, "sampleanno")
    writeData(wb, sheet=sheet, D %>% colData %>% as.data.frame, colNames = T, rowNames=F)
    saveWorkbook(wb, file_name, overwrite = TRUE)
  }
  
  return(D)
}

# read outcomes and type conversion 
outcomes <- readxl::read_excel(path= path_outcome_list, sheet="phenotypes2test")

#--- params
# exclude intubated samples
exclude.intubated <- T
# data date stamp
datestamp <- paste0("preprocessed","")  #,"_",format(Sys.time(),"%Y%m%d")
#---
dir = glue::glue("{root}/1_preprocessing")

# run the analysis and gather the results
nods <- c("metabo", "proteo") %>% {structure(., names=.)} %>%
  lapply( f_export, outcomes = outcomes, datestamp = datestamp, 
          exclude.intubated = exclude.intubated, dir = dir, 
          save_SE=save_SE, save_XLSX=save_XLSX)



# # check if pre processed data is identical to previous version
# #metabolomics
# load("data_GIM_final_dataset123_nointub.Rdata")
# Dnew = D
# load("~/Documents/R/_github_/GIM2/1_preprocessing/data_GIM_final_dataset123_nointub.Rdata")
# sum( assay(D) - assay(Dnew) )
# 
# #proteomics
# load("data_GIM_final_dataset1_proteo_nointub.Rdata")
# Dnew = D
# load("~/Documents/R/_github_/GIM2/1_preprocessing/data_GIM_final_dataset1_proteo_nointub.Rdata")
# sum( assay(D) - assay(Dnew) )
#

