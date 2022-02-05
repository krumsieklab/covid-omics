#' 
#' runs predictive modeling
#' 

# warning to the user, this might take very long
invisible(readline(prompt="\n!!!\nWarning: The run3_prediction.R script takes 18-24h to run on a 6-core Apple M1 CPU,\nand requires several gigabytes of memory!\nPress [enter] to continue."))


rm(list = setdiff(ls(),"root")) 
# root directory 
setwd(root)

# preprocessed metabolomics data 
path_preprocessed_metabo = 
  glue::glue("{root}/1_preprocessing/data_preprocessed_metabo.Rdata")

# minimal data ready for predictive modeling
path_data = glue::glue("{root}/3a_composite/data_for_predictive_modeling.Rdata")

# number of clusters for parallel running 
number_of_clusters = detectCores()-1

# working directory
wd = getwd()

# back up paths
objlist = ""; objlist = ls()

# create composite outcome  -----------------------------------------------
# check if the data is already ready
if( length(grep("predictive_modeling.Rdata", 
                list.files(path = glue::glue("{root}/3a_composite")))) == 0){
  # create composite outcome and save the data for predictive modeling 
  source(glue::glue("{root}/3a_composite/main_create_composite_outcome.R"))
}else{
  print("composite outcome created and data saved for predictive modeling")
}
setwd(wd)



# run predictive modeling with composite outcome --------------------------

# check if the model fits were already obtained 
# it can be time consuming to run the model everytime, so it makes sense to 
#   save the model fits to use later 
if( paste0("", grep("lasso_loo", list.files(glue::glue("{root}/3b_prediction")), value = T)) != 
    "lasso_loo_fits.Rdata"
  ){
  # run the predictive modeling with leave-one-out 
  source(glue::glue("{root}/3b_prediction/main_prediction_modeling.R"))
}else{
  print("pre-calculated model fits in the folder")
}
setwd(wd)

# check if the  lasso path  model fits were already obtained 
# it can be time consuming to run the model everytime, so it makes sense to 
#   save the model fits to use later 
if( paste0("",grep("lasso_path_loo", list.files(glue::glue("{root}/3b_prediction")), value = T)) != 
    "lasso_path_loo_fits.Rdata"
){
  # run the lasso path with leave-one-out 
  source(glue::glue("{root}/3b_prediction/main_prediction_modeling_lasso_path.R"))
}else{
  print("pre-calculated lasso path fits in the folder")
}
setwd(wd)

# check if the final model fits with whole data is obtained
# it can be very time consuming to run the model everytime, so it makes sense to 
#   save the model fit to use later 
if( paste0("",grep("signatures", list.files(glue::glue("{root}/3b_prediction")), value = T)) != 
    "reor4me_signatures.Rdata"
){
  # run the master fit with whole data 
  source(glue::glue("{root}/3b_prediction/main_prediction_master_fit.R"))
}else{
  print("model signatures in the folder")
}
setwd(wd)



