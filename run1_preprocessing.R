#' 
#' runs preprocessing scripts based on the data paths, 
#' and saves the preprocessed data as Rdata file or xlsx files
#' 

rm(list = setdiff(ls(),"root")) 
# root directory 
setwd(root)

# paths and arguments 
path_outcome_list = glue::glue("{root}/DATA/GIM_outcome_list.xlsx")

# paths where data resides in 
paths = list(
  clinical_metadata = glue::glue("{root}/DATA/gim_clinical_metadata.xlsx"), 
  metabo = glue::glue("{root}/DATA/gim_metabolomics.xlsx"),
  proteo = glue::glue("{root}/DATA/gim_proteomics.xlsx")
)

# to save preprocessed data as SE objects
save_SE = TRUE
# ro export proprocessed data as xlsx files
save_XLSX = FALSE

# working directory
wd = getwd()

# back up paths, settings etc
objlist = ""

objlist = ls()


# preprocess the data for general analysis --------------------------------
# check first if the data has already been generated
files = grep(pattern = "data_preprocessed", list.files(glue::glue("{root}/1_preprocessing")), value = T )
cat("Data is already preprocessed and saved, see: \n\t")
cat(files, sep = "\n\t")

if(length(files)<2){
  # preprocess the data for general analysis
  source(glue::glue("{root}/1_preprocessing/main_GIM_preprocessing.R"))
  # check if preprocessed data is produced
  grep(pattern = "data_preprocessed.*", 
       list.files(glue::glue("{root}/1_preprocessing")), value = T )
  # clean up for the next analysis 
  # rm(list = setdiff(ls(), objlist))
}
setwd(wd)

# preprocess the data for GGM ---------------------------------------------
# check first if the data has already been generated
# files = grep(pattern = "preprocessed_for_GGM", list.files(glue::glue("{root}/1_preprocessing")), value = T )
# cat(files, sep = "\n")
if(length(files)<1){
  # preprocess the data for GGM
  source(glue::glue("{root}/1_preprocessing/main_GIM_preprocessing_for_GGM.R"))
  # check if preprocessed data for GGM is produced
  grep(pattern = "preprocessed_for_GGM", 
       list.files(glue::glue("{root}/1_preprocessing")), value = T )
  # clean up for the next analysis 
  # rm(list = setdiff(ls(), objlist))
}
setwd(wd)

# preprocess the data for intubation-at-blood-draw ------------------------
# check first if the data has already been generated
# files = grep(pattern = "preprocessed_for_ibd", list.files(glue::glue("{root}/1_preprocessing")), value = T )
# cat(files, sep = "\n")
if(length(files)<1){
  # preprocess the data for intubation-at-blood-draw(ibd) analysis
  source(glue::glue("{root}/1_preprocessing/main_GIM_preprocessing_for_ibd.R"))
  # check if preprocessed data for ibd is produced
  grep(pattern = "preprocessed_for_ibd", 
       list.files(glue::glue("{root}/1_preprocessing")), value = T )
  # clean up for the next analysis 
  # rm(list = setdiff(ls(), objlist))
}
setwd(wd)
