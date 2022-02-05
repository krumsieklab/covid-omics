library(SummarizedExperiment)
library(magrittr)
library(survival)

# path_preprocessed_metabo = "~/GIM/1_preprocessing/data_preprocessed_metabo.Rdata"
# root = "~/GIM"
  
# load preprocessed metabolomics dataset 
load(path_preprocessed_metabo )
# set working directory to location of source code
setwd(glue::glue("{root}/3a_composite"))

D = D[,!is.na(D$death)]
D = D[,D$Status == "COVID"]


# prepare clinical variables ----------------------------------------------
# variables to be used to create composite outcome 
df = as.data.frame( colData(D) )
vars_of_interest = c(#"mv_at_discharge", 
  "death", "time_to_death",
  "rr_therapy", "pof_support", "intubation1", "renal_replacement", 
  "O2_device", "Fi_O2", "max_aki_stage", "hospitalization", "disposition")

df = df[,c("Patient_ID",vars_of_interest)] %>% unique()

# ordinal scale for each variable where higher is worse
df$death = as.numeric(df$death)
df$time_to_death = -df$time_to_death
df$rr_therapy = as.numeric(df$rr_therapy == "TRUE")
df$pof_support = as.numeric(df$pof_support == "TRUE")
df$intubation1 = as.numeric(df$intubation1)
df$renal_replacement = as.numeric(df$renal_replacement == "TRUE")
#'
#' RA	    room air
#' NCorTC	Nasal cannula or trach collar
#' NR	    non-rebreather mask
#' HFNC	  high-flow nasal cannula
#' NIMV	  non-invasive mechanical ventilation
#' MV	    mechanical ventilation
#'
df$O2_device = as.numeric(factor( 
  df$O2_device, levels = c("RA", "NCorTC", "NR", "HFNC", "NIMV", "MV") )
)
df$Fi_O2 = df$Fi_O2
df$max_aki_stage = df$max_aki_stage
df$hospitalization = df$hospitalization
#'
#' SNur	Skilled Nursing Facility
#' SReh	Subacute Rehab
#' AReh	Acute Rehab
#' Homx	Home
#'
df$disposition  = as.numeric(factor( 
  df$disposition, levels = c("Homx", "AReh", "SReh", "SNur") )
)

# conposite modeling starts here ------------------------------------------
source("utils_composite.R")
ron0 = f_infer_rankings_simple(df[,-1]) 

#r0 is composite score 
df$r0 = rank(ron0$rs)
compo = df

D = D[,D$Patient_ID %in% df$Patient_ID]
metab = data.frame(Patient_ID = D$Patient_ID, colData(D)[,c("age","sex", "bmi","Time")], as.matrix( assay(D)) %>% t)
metab$sex = as.numeric(metab$sex=="Male")

# metab: metabolomics data
# compo: composite outcomes with other variables
save(file = "data_for_predictive_modeling.Rdata", metab, compo)


