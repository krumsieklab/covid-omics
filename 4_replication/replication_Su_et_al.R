# Su et al --------------------------------------------------------------

root_dir = glue::glue("{root}/4_replication")
# set working directory to location of source code
setwd(root_dir)

library(magrittr)
library(survival)

# download the data  ------------------------------------------------------
# Su et al: https://doi.org/10.1016/j.cell.2020.10.037
url_data = "https://data.mendeley.com/public-files/datasets/tzydswhhb5/files/ead2d262-fb88-4bc7-b697-9b0f2bc7d4db/file_downloaded"
download.file(url_data, destfile = "Su_et_al_TableS1.xlsx")
url_data = "Su_et_al_TableS1.xlsx"
# -------------------------------------------------------------------------

df_su = read_xlsx(url_data, 2)
df_su = df_su[!is.na(df_su$`Who Ordinal Scale`),]
df_su$Sex = as.numeric(df_su$Sex != "Female")
df_su$wos = df_su$`Who Ordinal Scale` %>% {.[. %in% c("1","1 or 2","2")] = "2"; as.numeric(.)}

df_su <- df_su %>% dplyr::group_by(`Study Subject ID`) %>%
  dplyr::summarise( 
    subject_id = unique(`Study Subject ID`),
    age = mean(Age), 
    sex = mean(Sex), 
    bmi = mean(as.numeric(BMI), na.rm = T), 
    dwos = diff(wos[order(`Blood draw time point`)]),
    wost2 = wos[`Blood draw time point` == "T2"]
  ) %>% as.data.frame()

rownames(df_su) = df_su$subject_id

df_met = read_xlsx(url_data, 5) %>% as.data.frame
df_met = df_met %>% dplyr::filter(`Healthy donor sample or COVID19 sample`=="COVID19")
rownames(df_met) = df_met$sample_id
# impute missing values
df_met = as.matrix( df_met[,-(1:2)] ) %>% 
  apply(2,function(x){x[is.na(x)] = min(x,na.rm = T);x})

df_met = df_met[intersect(rownames(df_met), paste0(df_su$subject_id, "-AC")),]
rownames(df_met) = gsub(pattern = "-AC", "", rownames(df_met))

df_su = df_su[rownames(df_met),]

# our model and matched metabolites ---------------------------------------
load(glue::glue("{root}/3b_prediction/reor4me_signatures.Rdata"))
B0 = signatures$B_final_minimal
# matched metabo names 
selected_mets0 = names(B0)[-seq(3)]
met_match = readxl::read_xlsx(glue::glue("{root}/DATA/metabolite_annotations.xlsx"), "Su_et_al")
met_match = met_match[met_match$match, ] %>% {structure(.$target, names = .$us)}
met_match = met_match[selected_mets0]
# --------------------------------------- ---------------------------------

df1 = df_su
# bmi is missing for 1/3 of samples, don't use bmi, instead use age and sex
df1[,"bmi"] = predict(lm(bmi~age+sex,df1), newdata = df1)

# predicted score base
df1$yh_base = as.matrix( df1[,c("sex", "age", "bmi")] ) %*% signatures$B_base
concordance(wost2~yh_base, df1)

# predicted score base+mets
df1$yh_mode = as.matrix(cbind( df1[,c("sex","age","bmi")], df_met[,met_match])) %*%  B0
concordance(wost2~yh_mode, df1)

fsum <- function(dfit){
  d = dfit$concordance %>% unname
  se = sqrt(dfit$var) %>% unname
  c(d = d,se = se)
} 

library(ggplot2)
data.frame(
  rbind(
    concordance(wost2~yh_mode, df1) %>% fsum,
    concordance(wost2~yh_base, df1) %>% fsum
  ), model = c("final", "base") 
)  %>% ggplot(aes(x = model, y = d, color = model)) + 
  geom_pointrange(aes(ymin = d-se, ymax = d+se), lwd = 1.5) +
  geom_point(size = 6, color = "white", fill = "white") +
  geom_point(size = 3) +
  theme_minimal()  


# collect all predictions into one data.frame
df_su = df1[,c("subject_id","wost2","yh_base","yh_mode")]
colnames(df_su)[4] = "yh_mubu"

# model performances  
c(
  d_base = concordance(wost2~yh_base, df_su)$concordance,
  d_mubu = concordance(wost2~yh_mubu, df_su)$concordance
) %>% round(3)

save(file = "Su_et_al.Rdata",df_su)




