# Shen et al data1 replication
# Shen et al --------------------------------------------------------------

root_dir = glue::glue("{root}/4_replication")
# set working directory to location of source code
setwd(root_dir)

library(magrittr)
library(readxl)
library(survival)

# download the data  ------------------------------------------------------
# Shen et al paper: https://doi.org/10.1016/j.cell.2020.05.032
url_shen_et_al_clin = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420306279-mmc1.xlsx"
url_shen_et_al_data = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420306279-mmc2.xlsx"

download.file(url_shen_et_al_clin, destfile = gsub(pattern = ".*/", "" ,url_shen_et_al_clin) )
url_shen_et_al_clin = gsub(pattern = ".*/", "" ,url_shen_et_al_clin)
download.file(url_shen_et_al_data, destfile = gsub(pattern = ".*/", "" ,url_shen_et_al_data) )
url_shen_et_al_data = gsub(pattern = ".*/", "" ,url_shen_et_al_data)
# -------------------------------------------------------------------------

# clinical data
df_clin2 = read_xlsx(url_shen_et_al_clin, 2) %>% as.data.frame()
df_clin2$status = as.numeric( df_clin2$`Date of progression to severe state` != "/" ) 
df_clin2$`BMI h` = as.numeric(df_clin2$`BMI h`)

# test-set1 data extracted from the paper
# https://www.cell.com/action/showPdf?pii=S0092-8674%2820%2930627-9
# figure 2D, and figure 1A
df2_from_figure <- 
  data.frame(
    yh  = seq(10), # rank of estimated risk scores
    y   = c(rep(0,4), 1,0,1,0,1,1), # outcomes 
    PID = paste0("XG", c(24,20,23,21,45,22,46,25,44,43))
  )

df_clin2 = df_clin2[df_clin2$`Patient ID a` %in% df2_from_figure$PID,]

# metabolomics data for test set 1
df_met = read_xlsx(url_shen_et_al_data, "Metabolomics", skip = 1) %>% as.data.frame
rownames(df_met) = df_met$Metabolites
df_met = t(df_met[,-seq(1)])
# impute missing values as described in Shen et al 
df_met = apply(df_met,2,function(x){x[is.na(x)] = min(x,na.rm = T);x})

df_met = df_met[df_clin2$`Metabolomics ID e`,] 
rownames(df_met) = df_clin2$`Patient ID a`

df_clin2 = df_clin2[,c("Patient ID a", "COVID-19 patients corhort","Age (year)","Sex g","BMI h", "status")]
colnames(df_clin2) =  c("PID", "cohort" ,"age", "sex", "bmi", "status")

# check if matches with the data in the figure
identical(df2_from_figure$PID %>% sort, df_clin2$PID %>% sort)
identical(df2_from_figure[match(df_clin2$PID, df2_from_figure$PID),"y"], df_clin2$status)

# our model and matched metabolites ---------------------------------------
load(glue::glue("{root}/3b_prediction/reor4me_signatures.Rdata"))
B0 = signatures$B_final_minimal
# matched metabo names 
selected_mets0 = names(B0)[-seq(3)]
met_match = readxl::read_xlsx(glue::glue("{root}/DATA/metabolite_annotations.xlsx"),"Shen_et_al")
met_match = met_match[met_match$match, ] %>% {structure(.$target, names = .$us)}
met_match = met_match[selected_mets0]
# --------------------------------------- ---------------------------------

# predicted score with our model 
yh_d2 = as.matrix( cbind(df_clin2[,names(B0)[1:3]], df_met[,met_match] %>% log %>% scale) ) %*% B0

# predicted score with base model
yh_d2_base =  as.matrix(df_clin2[,c("sex","age","bmi")]) %*% signatures$B_base 

# shen et al model ranks of predicted scores 
yh_d2_shen = df2_from_figure$yh[match( df_clin2$PID,df2_from_figure$PID )]

# collect all prediction into one data.frame
df_shen_test1 = data.frame( 
  df_clin2[,c("PID","cohort","status")],
  yh_base = yh_d2_base, 
  yh_shen = yh_d2_shen,
  yh_mubu = yh_d2
)

# model performances  
c(
  d_base = concordance(status~yh_base, df_shen_test1)$concordance,
  d_shen = concordance(status~yh_shen, df_shen_test1)$concordance,
  d_mubu = concordance(status~yh_mubu, df_shen_test1)$concordance
) %>% round(3)

save(file = "Shen_et_al_testset1.Rdata",df_shen_test1)
