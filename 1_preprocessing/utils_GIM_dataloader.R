# GIM data loader functions
GIM_data_loader <- function(dataset, file_path, metadata_path, exclude.intubated = F){
  
  D = f_load(file_path)
  # read sample annotations
  df = readxl::read_xlsx(metadata_path, sheet = 1)
  df = dplyr::left_join( as.data.frame(colData(D)), df, "Patient_ID") 
  colnames(df)[colnames(df)== "intubation"] = "intubation1"
  identical( colnames(assay(D)), df$the_ID )
  colData(D) =  DataFrame(df)
  D$Group = ifelse(is.na(D$Status), NA, paste(D$Status, D$Time, sep = "_"))
  colnames(D) = D$the_ID
  D$age = round(D$age,2)
  D = if(exclude.intubated) D[,!D$already_intubated] else D
  
  if(dataset == "metabo") 
    return( lapply(c(D1="GIM1", D2="GIM2", D3="GIM3"), function(i) D[, D$batch==i] ) )
  
  if(dataset == "proteo") return( D )
  
  return(NULL)
}

# load data from xlsx to SE 
f_load <- function(file_path){
  sheets = readxl::excel_sheets(file_path)
  
  data = readxl::read_xlsx(file_path, sheet = "data")
  data = as.matrix(data[,-1]) %>% {rownames(.)= data$BIOCHEMICAL;.}
  metx = readxl::read_xlsx(file_path, sheet = grep(pattern = "metabolite|protein", x = sheets, value = T) )
  samx = readxl::read_xlsx(file_path, sheet = "sample_annotation")
  
  SummarizedExperiment::SummarizedExperiment(
    list(
      data
    ), 
    rowData = metx, 
    colData = samx 
  )
}


