library(magrittr)
library(survival)
library(coxme)

# load dataset 
load(path_data )
# set working directory to location of source code
setwd(glue::glue("{root}/3b_prediction"))

df = metab
df = dplyr::left_join(df, compo[,c("Patient_ID","r0")], by ="Patient_ID" )
# one bmi missing fill with conditional expectation
df[is.na(df$bmi),"bmi"] = 
  predict( lm(bmi~age+sex,df), newdata = df[is.na(df$bmi),c("age","sex"),drop=F] )
# remove NA r0 row 
colSums(is.na(df)) %>% {.[.!=0]}
dim(df)

mets = colnames(df)[6:130]
# remove samples which have no r0
df = df[!is.na(df$r0),]
# scale mets, coluld be important for rx prior ridge, this is different than 
# previous train-split analysis
df[,mets] = scale(df[,mets])

# Patient ID
colnames(df)[1] = "ID"

# # regularized linear mixed effect model
# functions needed to cast problem to be solved with cox(ph)/(me)
source("./model/utils_coxph_sorceries.R")
# ridge regularization 
source("./model/utils_coxme_sorceries.R")
# lasso regularization
source("./model/utils_lasso_sorceries.R")
# utils for lasso regularization path
source("./model/utils_lasso_path_sorceries.R")

# final lasso fit with whole data
fit  <- f_reor4me(df, "r0", mets, c("sex", "age", "bmi"), "Time", 
                  r_FF = r_OR, max_iter = 30, return_objs = F)

# final base model with whole data
fit0  <- r_OR(df, "r0", mets, c("sex","age","bmi"), "Time",  excludeM = T) 

# final lasso path with whole data
fit_path  <- f_reor4me_path(df, "r0", mets, c("sex", "age", "bmi"), "Time", 
                       r_FF = r_OR, N_iter = 550, pmf = 0.99, return_objs = F)


# save final model fitted with whole data  --------------------------------

# base model signature 
B_base = fit0$coefficients

# final model signature
B_final = fit$B1[fit$B1 !=0]

# minimal model signature see supp fig. 2
#' most of the predictive power (d ~ 0.68) was already achieved with
#' the first 6 metabolites + 3 base variables, 
#' after which prediction performance did not significantly improve
B_final_minimal = fit_path$B_path %>% 
  {.[,max(which(colSums(. != 0)==9))]} %>% {.[.!=0]}

fits = list(base = fit0, 
                   lasso = fit,
                   lasso_path = fit_path)

signatures = list(B_base = B_base, 
                  B_final = B_final,
                  B_final_minimal = B_final_minimal)  

save(file = "reor4me_signatures.Rdata", fits, signatures)
