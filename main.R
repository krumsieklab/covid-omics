rm(list = ls())

# set the root folder where all codes are located
root = dirname(rstudioapi::getActiveDocumentContext()$path) # "covid-omics"
setwd(root)

#'
#' renv package is needed to restore 
#' the package versions that are used in the code
#'  
renv::restore()

#' 
#' run pre-processing and get the preprocessed data
#' 1_preprocessing: preprocesses and saves the data 
#' 
source("run1_preprocessing.R")

#' 
#' run statistical association and GGM analysis
#' 2a_casecontrol: covid vs control analysis 
#' 2b_clinical  : association analysis with clinical variables
#' 2c_networks  : GGM network analysis
#' 2d_supplements: intubation at blood draw analysis
#'                PCAs and Time effect analysis
#' 
source("run2_allanalyses.R")

#' 
#' run predictive modeling and replication of the 
#'   model with independent datasets
#' 3a_composite  : creates composite outcomes
#' 3b_prediction : runs predictive modeling
#' 
source("run3_prediction.R")

#'
#'  4_replication: checks developed model on the independent dataset
#' 
source("run4_replication.R")

#' 
#' creates paper figures
#'
source("run5_figures.R")




