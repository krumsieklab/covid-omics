# GIM, effect of intubation at blood draw

library(maplet)
library(tidyverse)
library(magrittr)
library(glue)

# root = "/Users/mub4004/Documents/R/_github_/GIM2" # change later 

# load the preprocessed data 
load(path_preprocessed_ibd_data)

# set working directory to location of source code
setwd(glue::glue("{root}/2d_supplements"))

# adjusted p cut-off
p.adj.cut <- 0.05

# confounder formula
formula_conf <- "Time" 

# data stamp
datestamp <- format(Sys.time(),"%Y%m%d")

# run analysis ------------------------------------------------------------
# workhorse function
source("utils_GIM_intubation_at_blood_draw.R")

# run the analysis and gather the results
nods <- list(metabo = list(D = Ds$metabo, dataset_name = "metabo"), 
             proteo = list(D = Ds$proteo, dataset_name = "proteo")) %>%
  lapply( f_WORKHORSE_int_at_bld_drw, 
          datestamp = datestamp,
          formula_conf = formula_conf,
          p.adj.cut = p.adj.cut,
          write_results = T,
          output_html = T)

# 
# root_dir = getwd() #"~/Box/results/COVID19/Choi_COVID19/GIM/most_recent/"
# dir.create("images")
# 
# # save the image of the results 
# save.image(glue::glue("{root_dir}/images/GIM_intubation_at_blood_draw_{datestamp}.RData"))
# 
# volcano plots
# ggs <- lapply(nods, function(x)
#   metadata(x$D)$results %>% {.[[grep("volcano", x = names(.))]]$output[[1]]}  +
#     labs(color = "p.adj<0.05")
#   )
#
