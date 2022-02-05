library(parallel)
library(magrittr)
library(dplyr)
library(survival)
library(coxme)
library(bdsmatrix)

# number of clusters for parallel running 
# number_of_clusters = 37

# path_data = "data_for_predictive_modeling.Rdata"
# root = "~"

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

# scale for spurious intercept effect
colnames(df)[1] = "ID"

# # regularized linear mixed effect model
# functions needed to cast problem to be solved with cox(ph)/(me)
source("./model/utils_coxph_sorceries.R")
# ridge regularization 
source("./model/utils_coxme_sorceries.R")
# lasso regularization
source("./model/utils_lasso_sorceries.R")

# prepare the setting for parallel running
PARPRE <- function(ncl = 3){
  flist = ls(envir = globalenv()) 
  # logfile
  rand_txt <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  }
  
  cl = makeCluster(ncl, outfile = tolower(rand_txt(6)), setup_strategy = "sequential")
  clusterEvalQ(cl, library(survival))
  clusterEvalQ(cl, library(bdsmatrix))
  clusterEvalQ(cl, library(coxme))
  clusterEvalQ(cl, library(magrittr))
  clusterExport(cl, flist)
  cl
}

# variables to be modelled as outcome
vv = c("r0") %>% {names(.)=.;.}

# create the clusters
cl = PARPRE(number_of_clusters)

# loo-cross-validation ----------------------------------------------------
# models with metabolites
nods_loo <-  
  lapply(vv, function(i){
    # # this does not work here no excludeM in reor4me
    # lapply(list(basic= T, fullm = F), function(flag) 
    env = environment()
    clusterExport(cl, "i", envir = env)
    re <- parLapply(cl, unique(df$ID), function(j) try({
      set.seed(j)
      print(j)
      # training set samples 
      trids = df$ID != j
      trids = trids & !is.na(df[,i])
      
      # scale the metabolites
      df_ = df[trids,]
      df_[,mets] = scale(df_[,mets], center  = T, scale = F)
      
      fit <- 
        if(length(unique(df_[,i])) < 3 ) # logistic regression
          f_reor4me(df_, i, mets, c("sex","age","bmi"), "Time",  r_FF = r_LR, max_iter = 30) 
        else # ordinal regression 
          f_reor4me(df_, i, mets, c("sex","age","bmi"), "Time", r_FF = r_OR, max_iter = 30) 
      
      list(fit = fit, trids = trids)
    }))
    print(i)
    re
  })

# base model
nods_base_loo <-  
  lapply(vv, function(i){
    # # this does not work here no excludeM in reor4me
    # lapply(list(basic= T, fullm = F), function(flag) 
    env = environment()
    clusterExport(cl, "i", envir = env)
    re <- parLapply(cl, unique(df$ID), function(j) try({
      set.seed(j)
      print(j)
      # train set samples 
      trids = df$ID != j
      trids = trids & !is.na(df[,i])
      
      # scale the metabolites
      df_ = df[trids,]
      df_[,mets] = scale(df_[,mets], center  = T, scale = F)
      
      fit <- 
        if(i %in% c("intubation","death","renal_replacement") )
          r_LR(df_, i, mets, c("sex","age","bmi"), "Time",  excludeM = T) 
        else 
          r_OR(df_, i, mets, c("sex","age","bmi"), "Time",  excludeM = T) 
      
      list(fit = ff_koni(fit), trids = trids)
    }))
    print(i)
    re
  })
stopCluster(cl)

# save model fits 
save(file = "lasso_loo_fits.Rdata", nods_loo, nods_base_loo, df, compo)



