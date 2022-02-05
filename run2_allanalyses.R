#' 
#' runs association, GGM analysis, and creates composite outcome
#' 

rm(list = setdiff(ls(),"root")) 
# root directory 
setwd(root)

# outcome list
path_outcome_list = glue::glue("{root}/DATA/GIM_outcome_list.xlsx")

# preprocessed metabolomics data 
path_preprocessed_metabo = 
  glue::glue("{root}/1_preprocessing/data_preprocessed_metabo.Rdata")

# preprocessed proteomics data 
path_preprocessed_proteo =
  glue::glue("{root}/1_preprocessing/data_preprocessed_proteo.Rdata")

# preprocessed data for GGM analysis
path_preprocessed_GGM_data = 
  glue::glue("{root}/1_preprocessing/data_preprocessed_for_GGM.Rdata")

# preprocessed data for intubation-at-blood-draw analysis
path_preprocessed_ibd_data = 
  glue::glue("{root}/1_preprocessing/data_preprocessed_for_ibd.Rdata")

# working directory
wd = getwd()

# back up paths
objlist = ""; objlist = ls()

# run the case-control analysis -------------------------------------------
# check if the image of the results was already created
if( !dir.exists(glue::glue("{root}/2a_casecontrol/results") )){
  # run the case-control analysis 
  source(glue::glue("{root}/2a_casecontrol/main_GIM_casecontrol.R"))
}else{
  cat("case-control analysis has already been executed, results are in: \n\t")
  cat(list.files( glue::glue("{root}/2a_casecontrol/results")), sep = "\n\t")
  cat("\n")
}
setwd(wd)

# run clinical variables association analysis -----------------------------
# check if the image of the results was already created
if( !dir.exists(glue::glue("{root}/2b_clinical/results") )){
  # run the clinicals analysis 
  source(glue::glue("{root}/2b_clinical/main_GIM_clinical.R"))
}else{
  cat("clinical analysis has already been executed, results are in: \n\t")
  cat(list.files( glue::glue("{root}/2b_clinical/results")), sep = "\n\t")
  cat("\n")
}
setwd(wd)

# run GGM analysis --------------------------------------------------------
# check if the image of the results was already created
if( !dir.exists(glue::glue("{root}/2c_networks/results"))){
  # run the GGM analysis 
  source(glue::glue("{root}/2c_networks/main_GGM_analysis.R"))
}else{
  cat("GGM analysis has already been executed, results are in: \n\t")
  cat(list.files( glue::glue("{root}/2c_networks/results")), sep = "\n\t")
  cat("\n")
}
setwd(wd)

# intubated-at-blood-draw-samples not excluded  ---------------------------
# check if the image of the results was already created
if( !file.exists(glue::glue("{root}/2d_supplements/results/GIM_proteointubation_at_blood_draw.rds")) ){
  # run the analysis with intubated samples not excluded
  source(glue::glue("{root}/2d_supplements/main_GIM_intubation_at_blood_draw.R"))
}else{
  cat("Intubation-at-blood-draw analysis has already been executed, results are in: \n\t")
  cat(list.files( glue::glue("{root}/2d_supplements/results")), sep = "\n\t")
  cat("\n")
}
setwd(wd)

# supplements -------------------------------------------------------------
# check effect of Time also look at PCAs ----------------------------------
if( !file.exists(glue::glue("{root}/2d_supplements/results/supplement_PCA_Time.fig"))){
  source(glue::glue("{root}/2d_supplements/main_PCA_and_Time_analysis.R"))
}else{
  cat("PCA and Time analysis were already done, and saved in \n\t 'supplement_PCA_Time.fig'")
  cat("\n")
}
setwd(wd)




