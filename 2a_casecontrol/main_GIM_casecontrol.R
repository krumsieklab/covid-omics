# GIM, case vs. control analysis using LME models
library(maplet)
library(tidyverse)
library(magrittr)
library(glue)
library(SummarizedExperiment)

# `root` directory and paths for processed data are needed to run this script
# load preprocessed datasets 
load(path_preprocessed_metabo)
D_met = D
load(path_preprocessed_proteo)
D_pro = D
rm(D)

root_dir = glue::glue("{root}/2a_casecontrol")
# set working directory to location of source code
setwd(root_dir)

# Init --------------------------------------------------------------------
# params
# exclude already intubated samples
exclude.intubated <- T
# adjusted p cut-off
p.adj.cut <- 0.05
# confounder formula
formula_conf <- "Time" #"sex + age + Time"
# data stamp
datestamp <- format(Sys.time(),"%Y%m%d")

# workhorse function
source("utils_GIM_casecontrol.R")

# run analysis ------------------------------------------------------------
nods <- list(metabo = list(D = D_met, dataset_name = "metabo"), 
             proteo = list(D = D_pro, dataset_name = "proteo")) %>%
  lapply( f_WORKHORSE, #COVID19_datapath = COVID19_datapath, 
          datestamp = datestamp,
          exclude.intubated = exclude.intubated, formula_conf = formula_conf, p.adj.cut = p.adj.cut,
          write_results = T,
          output_html = T)

# # save the results image
# dir.create("images")
# # save the image of the results
# save.image(glue::glue("{root_dir}/images/GIM_case_control_{datestamp}.RData"))

