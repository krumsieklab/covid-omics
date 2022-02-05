#' 
#' runs replication analysis on data from
#' 
#' Shen et al., https://doi.org/10.1016/j.cell.2020.05.032
#' Su et al., https://doi.org/10.1016/j.cell.2020.10.037
#'
 
rm(list = setdiff(ls(),"root")) 
# root directory 
# root = "~/Documents/R/_github_/covid-omics/"
setwd(root)

# back up paths
objlist = ""; objlist = ls()

# performance of our model on Shen et al data  ----------------------------
if(!file.exists("4_replication/Shen_et_al_testset1.Rdata")){
  # Shen et al testset1 aka C2
  source(glue::glue("{root}/4_replication/replication_Shen_et_al1.R"))
  # clean up for the next analysis 
  rm(list = setdiff(ls(), objlist))
  setwd(root)
  
  # Shen et al testset2 aka C3
  source(glue::glue("{root}/4_replication/replication_Shen_et_al2.R"))
  # clean up for the next analysis 
  rm(list = setdiff(ls(), objlist))
  setwd(root)
  
  # reanalysis of Shen et al model showing that it is no better than age 
  source(glue::glue("{root}/4_replication/reanalysis_Shen_et_al.R"))
  # clean up for the next analysis 
  rm(list = setdiff(ls(), objlist))
  setwd(root)
}else{
  print("Shen et al data already retrieved online and saved in the folder, see:")
  cat("\t Shen_et_al_testset1.Rdata \n")
  cat("\t Shen_et_al_testset2.Rdata \n")
}

# performance of our model on Su et al data  ------------------------------
if(!file.exists("4_replication/Su_et_al.Rdata")){
  source(glue::glue("{root}/4_replication/replication_Su_et_al.R"))
  # clean up for the next analysis 
  rm(list = setdiff(ls(), objlist))
  setwd(root)
}else{
  print("Su et al data already retrieved online and saved in the folder, see:")
  cat("\t Su_et_al.Rdata \n")
}






