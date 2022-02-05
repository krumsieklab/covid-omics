# GIM: Perform clinical association analysis within the COVID+ samples

library(SummarizedExperiment)
library(maplet)
library(tidyverse)
library(magrittr)
library(glue)

# `root` directory and paths for processed data are needed to run this script
# load preprocessed datasets 
load(path_preprocessed_metabo)
D_met = D
load(path_preprocessed_proteo)
D_pro = D
rm(D)

# set working directory
root_dir = glue::glue("{root}/2b_clinical")
# set working directory to location of source code
setwd(root_dir)

# path_outcome_list = glue::glue("{root}/DATA/GIM_outcome_list.xlsx")
path_outcome_list = path_outcome_list 

# Init --------------------------------------------------------------------
# params
case_only <- T
exclude.intubated <- T

# confounder formula
formula_conf <- "Time" # "sex + age + Time"

# stat params
do.boxplots <- T
p.adj.cut <- 0.05
pboxcut <- 0.00001

# date of the analysis
datestamp <- format(Sys.time(),"%Y%m%d")

# outcomes 
outcomes <- readxl::read_excel(path = path_outcome_list , sheet="phenotypes2test")

# workhorse function
source("utils_GIM_clinical.R")

# run analysis ------------------------------------------------------------

nods <- list(metabo = list(D = D_met, dataset_name = "metabo"), 
             proteo = list(D = D_pro, dataset_name = "proteo")) %>%
  lapply( f_WORKHORSE, #COVID19_datapath = COVID19_datapath, 
          outcomes = outcomes,
          datestamp = datestamp,
          case_only = case_only,
          exclude.intubated = exclude.intubated, formula_conf = formula_conf, p.adj.cut = p.adj.cut,
          do.boxplots = do.boxplots, pboxcut = pboxcut,
          write_results = T,
          output_html = T)

# # for memory efficiency, don't save ggplot environment 
# nods = lapply(nods, function(x){
#   x$ggs_equalizer %<>% lapply(maplet:::fix_ggplot_env)
#   x
# })
# 
# root_dir = getwd()
# dir.create("images")
# # save the image of the results
# save.image(glue::glue("{root_dir}/images/GIM_clinicals_{datestamp}.RData"))


