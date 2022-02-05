# dataset 3 test set 2
# Shen et al --------------------------------------------------------------

root_dir = glue::glue("{root}/4_replication")
# set working directory to location of source code
setwd(root_dir)

library(magrittr)
library(readxl)
library(survival)
library(glmnet)

# download the data  ------------------------------------------------------
# Shen et al paper: https://doi.org/10.1016/j.cell.2020.05.032
url_shen_et_al_clin = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420306279-mmc1.xlsx"
url_shen_et_al_data = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420306279-mmc2.xlsx"
url_shen_et_al_test = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420306279-mmc5.xlsx" 

download.file(url_shen_et_al_clin, destfile = gsub(pattern = ".*/", "" ,url_shen_et_al_clin) )
url_shen_et_al_clin = gsub(pattern = ".*/", "" ,url_shen_et_al_clin)
download.file(url_shen_et_al_data, destfile = gsub(pattern = ".*/", "" ,url_shen_et_al_data) )
url_shen_et_al_data = gsub(pattern = ".*/", "" ,url_shen_et_al_data)
download.file(url_shen_et_al_test, destfile = gsub(pattern = ".*/", "" ,url_shen_et_al_test) )
url_shen_et_al_test = gsub(pattern = ".*/", "" ,url_shen_et_al_test)
# -------------------------------------------------------------------------

# proteomics data
df_prot = read_xlsx(url_shen_et_al_data, "Proteomics_proteins_training", skip = 1)
df_prot = as.data.frame(df_prot)
cn = df_prot$...1
rn = colnames(df_prot)[-(1:2)]
df_prot = t(df_prot[,-seq(2)])
df_prot = as.matrix(df_prot) %>% apply(2, as.numeric)
rownames(df_prot) = rn
colnames(df_prot) = cn
rm(rn, cn)

# metabolomics data 
df_met = read_xlsx(url_shen_et_al_data, "Metabolomics", skip = 1) %>% as.data.frame()
rownames(df_met) = df_met$Metabolites
df_met = t(df_met[,-seq(1)])

# impute missing values as described in Shen et al 
df_prot = apply(df_prot,2,function(x){x[is.na(x)] = min(x,na.rm = T);x})
df_met = apply(df_met,2,function(x){x[is.na(x)] = min(x,na.rm = T);x})

# clinial data 
df_clin = read_xlsx(url_shen_et_al_clin, 2) %>% 
  as.data.frame() %>% dplyr::filter(`MS ID b` != "/" & `Metabolomics ID e` != "/" )

# available clinical data for overlap patients 
ids = intersect(
  match( rownames(df_prot), df_clin$`MS ID b` ),
  match( rownames(df_met), df_clin$`Metabolomics ID e`)
) %>% na.omit %>% unique 
df_clin = df_clin[ids,]
df_prot = df_prot[match( df_clin$`MS ID b`, rownames(df_prot) ),]
df_met = df_met[match( df_clin$`Metabolomics ID e`, rownames(df_met) ),]

# all are in same order
identical(rownames(df_met), df_clin$`Metabolomics ID e`)
identical(rownames(df_prot), df_clin$`MS ID b`)
rownames(df_clin) <- rownames(df_prot) <- rownames(df_met) <- df_clin$`Patient ID a`
df_clin$`COVID-19 patients corhort`
df_clin$status = as.numeric( df_clin$`Date of progression to severe state` != "/" ) 
df_clin = df_clin[,c("Patient ID a", "COVID-19 patients corhort","Age (year)","Sex g","BMI h", "status")]
colnames(df_clin) =  c("PID", "cohort" ,"age", "sex", "bmi", "status")
df_clin$bmi = as.numeric(df_clin$bmi)
# sex=1 -> male # df_clin[df_clin$PID == "XG40",]

# training set 
df_clin = df_clin[df_clin$cohort %in% "Training (C1)", ]
df_met = df_met[df_clin$PID,]
df_prot = df_prot[df_clin$PID,]
rm(ids)

# # dataset 3 -------------------------------------------------------------
df3 = read_xlsx(url_shen_et_al_test, 2, skip = 1)
rn = colnames(df3)[-seq(2)]
cn = df3$...1
df3 = t(df3[,-(1:2)])
rownames(df3) = gsub(pattern = "_", "-", rn )
colnames(df3) = cn
rm(cn,rn)
df3 = scale(log(df3))

# clin data of d3
df_clin3 = read_xlsx(url_shen_et_al_clin, 2) %>% 
  as.data.frame() %>% dplyr::filter(`MS ID b` != "/" & `Metabolomics ID e` != "/" )
df_clin3 = df_clin3[match(rownames(df3), df_clin3$`Patient ID a`),]
identical(df_clin3$`Patient ID a`, rownames(df3))

df_clin3$status = as.numeric( df_clin3$`Date of progression to severe state` != "/" ) 
df_clin3 = df_clin3[,c("Patient ID a", "COVID-19 patients corhort","Age (year)","Sex g","BMI h", "status")]
colnames(df_clin3) =  c("PID", "cohort" ,"age", "sex", "bmi", "status")
df_clin3$bmi = as.numeric(df_clin3$bmi)

# our model and matched metabolites ---------------------------------------
load(glue::glue("{root}/3b_prediction/reor4me_signatures.Rdata"))
B0 = signatures$B_final_minimal
# matched metabo names 
selected_mets0 = names(B0)[-seq(3)]
met_match = readxl::read_xlsx(glue::glue("{root}/DATA/metabolite_annotations.xlsx"),"Shen_et_al")
met_match = met_match[met_match$match, ] %>% {structure(.$target, names = .$us)}
met_match = met_match[selected_mets0]
# --------------------------------------- ---------------------------------

# find the signal coming from our model with the available measured molecules
# reference space in training data
dt_reference = cbind(df_clin[,c("sex","age","bmi")], df_met[,met_match] %>% log %>% scale) %>% as.matrix
# target space in training data
dt_target = cbind( as.matrix(df_clin[,c("sex","age","bmi")]), cbind( df_met, df_prot)[,colnames(df3)] %>% log %>% scale)

# score calculated with our model 
yh = dt_reference %*% B0
# model to predict our score
fit_latent = cv.glmnet(y=yh, x = dt_target , alpha = 0, nfolds = length(yh))

# our score through latent fit in dataset 3
dt3 = cbind(df_clin3[,c("sex","age","bmi")] %>% as.matrix, df3 )

# our model predicted scores
yh_d3 =  predict(fit_latent, dt3) 

# base model predicted scores
yh_d3_base =  as.matrix(df_clin3[,c("sex","age","bmi")]) %*% signatures$B_base 

# test-set2 data extracted from the paper
# https://www.cell.com/action/showPdf?pii=S0092-8674%2820%2930627-9
# figure 2E, and figure 1A
df3_from_figure <- 
  data.frame(
    yh = seq(19), # rank of estimated risk scores
    y  = c(rep(0,9), 1,0,0,rep(1,6),0), # outcomes 
    PID = paste0("X2-",c(11,3,9,21,23,2,8,20,28,18,24,26,13,16,19,12,14,7,22))
  )
# check if all matches 
identical(df3_from_figure$PID %>% sort, df_clin3$PID %>% sort)
identical(df3_from_figure[match(df_clin3$PID, df3_from_figure$PID),"y"], df_clin3$status)

# shen et al model ranks of predicted scores 
yh_d3_shen = df3_from_figure$yh[match( df_clin3$PID, df3_from_figure$PID )]

# collect all prediction into one data.frame
df_shen_test2 = data.frame( 
  df_clin3[,c("PID","cohort","status")],
  yh_base = yh_d3_base, 
  yh_shen = yh_d3_shen,
  yh_mubu = as.numeric(yh_d3)
)


# model performances  
c(
  d_base = concordance(status~yh_base, df_shen_test2)$concordance,
  d_shen = concordance(status~yh_shen, df_shen_test2)$concordance,
  d_mubu = concordance(status~yh_mubu, df_shen_test2)$concordance
) %>% round(3)

save(file = "Shen_et_al_testset2.Rdata",df_shen_test2)




