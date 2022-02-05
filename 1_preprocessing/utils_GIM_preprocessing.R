# GIM preprocessing codes
# previously GIM_1_data_loader.R

# For processing any of the GIM datasets
# It calls function GIM_preprocessing_default for processing dataset1/2/3 (a single SE)
# It calls function GIM_preprocessing_dataset23 for processing dataset23 (a list of 2 SEs)
# @param D: a single SE or a list of two SEs
# @param anchor_method: optional. If input D is a list of two SEs, it can be used to define method of choice for anchor normalization
# @param exclude.pregnant: optional flag to exclude pregnant samples from data
# @returns D: a processed summarized experiment object
GIM_preprocessing <- function(D, exclude.pregnant=T, anchor_method=NULL){
  
  if("SummarizedExperiment" %in% class(D)){
    if (exclude.pregnant) {
      D %<>% mt_modify_filter_samples(filter = (pregnancy=="No"))
    }
    D <- GIM_preprocessing_default (D)
  } else if (("list" %in% class(D)) && length(D)==2) {
    if(("SummarizedExperiment" %in% class(D[[1]])) && ("SummarizedExperiment" %in% class(D[[2]]))){
      D <- GIM_preprocessing_dataset23 (D, exclude.pregnant=exclude.pregnant, 
                                        anchor_method=anchor_method) 
    }
  } else if (("list" %in% class(D)) && length(D)==3) {
    if(("SummarizedExperiment" %in% class(D[[1]])) && ("SummarizedExperiment" %in% class(D[[2]]))
       && ("SummarizedExperiment" %in% class(D[[3]]))){
      D <- GIM_preprocessing_dataset123 (D, exclude.pregnant=exclude.pregnant) 
    }
  } else {stop("Unknown input data format!")}
  
  return(D)
}


# For processing any of the GIM datasets. It is called by function GIM_preprocessing and by function GIM_processing_dataset23
# It calls function normalize_by_deviation if a reference SE is provided
# @param D: SE
# @param refD: optional reference SE to be used for anchor normalization of D
# @param anchor_method: optional method to be used for anchor normalization of D given refD
# @param to_impute: optional flag to set imputation T/F
# @returns D: a proccessed SE
GIM_preprocessing_default <- function(D, refD=NULL, anchor_method=NULL, to_impute=T){
  
  D <- D %>% 
    # header
    mt_reporting_heading("Preprocessing") %>%
    # missing value filtering
    mt_plots_missingness(feat_max = 0.25, plot_features = T, plot_samples = T, plot_data = T) %>%
    mt_pre_filter_missingness(samp_max=0.5) %>%
    mt_plots_missingness(feat_max=0.25, plot_features = T, plot_samples = T, plot_data = T) %>%
    mt_pre_filter_missingness(feat_max=0.25) %>%
    mt_plots_missingness(feat_max=0.25, plot_features = T, plot_samples = T, plot_data = T) %>%
    # normalization
    mt_plots_sample_boxplot(plot_logged=T, color=Group) %>%
    mt_pre_norm_quot(feat_max = 0) %>%
    mt_plots_dilution_factor(in_col = "Group") %>%
    mt_plots_sample_boxplot(plot_logged=T, color=Group) %>%
    # transformation & imputation
    mt_pre_trans_log()
  
  # if there is a reference D switch to anchor based normalization
  if("SummarizedExperiment" %in% class(refD)){
    # filter and quotient normalize reference SE
    refD <- refD %>% mt_pre_filter_missingness(feat_max=0.25) %>% 
      mt_pre_norm_quot(feat_max = 0) %>% mt_pre_trans_log()
    D %<>% normalize_by_deviation (refD, anchor_method=anchor_method)
  } 
  
  if (to_impute){
    D %<>% mt_pre_impute_knn() 
  }
  # pathway annotations, if there is an HMDB column
  if ("HMDB" %in% (D %>% rowData() %>% colnames())) {
    D %<>%  
      # HMDB to KEGG
      mt_anno_hmdb_to_kegg(
        in_col = 'HMDB',
        out_col = 'KEGG_mapped'
      ) %>%
      # pull pathway annotations
      mt_anno_pathways_hmdb(
        in_col = "HMDB", 
        out_col = "kegg_db", 
        pwdb_name = "KEGG") %>% 
      mt_anno_pathways_remove_redundant(feat_col = "HMDB", pw_col = "kegg_db") %>%
      mt_write_pathways(pw_col = "kegg_db", file="results/pw_annos_kegg.xlsx")
  }
  # pathway annotations, if there is an UniProt column
  if ("UniProt" %in% (D %>% rowData() %>% colnames())) {
    D %<>%  
      # pull pathway annotations
      mt_anno_pathways_uniprot(
        in_col = "UniProt", 
        out_col = "kegg_db") %>% 
      mt_anno_pathways_remove_redundant(feat_col = "UniProt", pw_col = "kegg_db") %>%
      mt_write_pathways(pw_col = "kegg_db", file="results/pw_annos_kegg.xlsx")
    
    # average duplicate molecules --> specifically needed for proteins
    D %<>% mt_modify_avg_features(group_col = 'name')  
  }
  
  # final dataset info
  D %<>% mt_reporting_data()
  
  return(D)
}


# For processing GIM dataset23. It is called by function GIM_preprocessing
# It calls function GIM_preprocessing_default for normalization of the two SEs one by one
# @param D: (dataset3) SE 
# @param refD: (dataset2) a second SE to be used for anchor normalization of D
# @param anchor_method: optional, method to be used for anchor normalization of D given refD
# @param exclude.pregnant: optional flag to exclude pregnant samples from data
# @returns D: a processed combined SE
GIM_preprocessing_dataset23 <- function(Ds, exclude.pregnant=T, anchor_method=NULL){
  
  D <- Ds$D
  refD <- Ds$refD
  # exclude pregnant except if it is an anchor
  if (exclude.pregnant) {
    D %<>% mt_modify_filter_samples(filter = (pregnancy=="No" | anchor=="Yes"))
    refD %<>% mt_modify_filter_samples(filter = (pregnancy=="No" | anchor=="Yes"))
  }
  # preprocess D (dataset3)
  D %<>% GIM_preprocessing_default (refD, anchor_method = anchor_method, to_impute = F)
  # preprocess refD (dataset2)
  refD %<>% GIM_preprocessing_default(to_impute = F)
  ##### prepare to merge D AND refD
  # exclude anchor samples from D they will be included from refD
  D %<>% mt_modify_filter_samples(filter = (anchor=="No"))
  # add batch info
  refD %<>% mt_anno_mutate(anno_type = 'samples', col_name='batch', term='GIM2')
  D %<>% mt_anno_mutate(anno_type = 'samples', col_name='batch', term='GIM3')
  # check for common metabolites
  common_mets <- intersect(rowData(refD)$name, rowData(D)$name)
  # check for order of metabolites
  m1 <- match(common_mets, rowData(refD)$name)
  m2 <- match(common_mets, rowData(D)$name)
  # finally combine
  Dcomb <- SummarizedExperiment(
    assay    = cbind(assay(refD)[m1, ], assay(D)[m2, ]),
    colData  = plyr::rbind.fill(data.frame(colData(refD)), data.frame(colData(D))),
    rowData  = rowData(refD)[m1, c('name', 'HMDB')],
    metadata = list()
  )
  # impute the datasets together and flag logged it is lost while combining
  Dcomb %<>% mt_pre_impute_knn() %>% mt_flag_logged()
  # exclude pregnant anchor samples
  if (exclude.pregnant) {
    Dcomb %<>% mt_modify_filter_samples(filter = (pregnancy=="No"))
  }
  # pathway annotations, if there is an HMDB column
  if ("HMDB" %in% (Dcomb %>% rowData() %>% colnames())) {
    Dcomb %<>%  
      # HMDB to KEGG
      mt_anno_hmdb_to_kegg(
        in_col = 'HMDB',
        out_col = 'KEGG_mapped'
      ) %>%
      # pull pathway annotations
      mt_anno_pathways_hmdb(
        in_col = "HMDB", 
        out_col = "kegg_db", 
        pwdb_name = "KEGG") %>% 
      mt_anno_pathways_remove_redundant(feat_col = "HMDB", pw_col = "kegg_db") %>%
      mt_write_pathways(pw_col = "kegg_db", file="results/pw_annos_kegg.xlsx")
  }
  
  # final dataset info
  Dcomb %<>% mt_reporting_data()
  return(Dcomb)
}


# For anchor normalization of dataset3. It is called by function GIM_preprocessing_default. 
# It calls function get_anchor_ratios
# It is hard coded to remove two "outlier metabolites" we observed during optimization
# @param D: SE (dataset3)
# @param refD: SE (dataset2)
# @param anchor_method: optional to set method for anchor norm
# @return D: an anchor normalized SE (dataset3)
normalize_by_deviation <- function(D, refD, anchor_method=NULL){
  
  if(is.null(anchor_method)){
    anchor_method <- "imputed_ratios"
  }
  # remove the outlier metabolites
  to_rem <- c("Phosphorylcholine", "Adenosine monophosphate")
  D %<>% mt_modify_filter_features(name%in%to_rem==F)
  
  # ratio of batch3 / batch2
  ratio_mat <- get_anchor_ratios(D, refD)
  # column averages
  col_avg <- colMeans(ratio_mat, na.rm = T)
  if(anchor_method=="mean_of_mean"){ # mean of mean ratios --> can be removed
    X <- assay(D)
    # subtract the mean of mean ratios
    X <- X - mean(col_avg)
  } else if (anchor_method=="mean"){ # mean ratios --> can be removed
    # subset metabolites similar to batch2
    D %<>% mt_modify_filter_features(name%in%names(col_avg))
    # order metabolites
    D <- D[match(names(col_avg), rowData(D)$name), ]
    X <- assay(D)
    # subtract the ratio avg 
    X <- apply(X, 2, FUN=function(x) x-col_avg)
  } else if (anchor_method=="imputed_ratios"){ # mean ratios no subset
    mm <- mean(col_avg) # mean of means
    # mean of remaining metabolites
    rem_avg <- rep(mm, length(which(rowData(D)$name%in%names(col_avg)==F)))
    names(rem_avg) <- rowData(D)$name[which(rowData(D)$name%in%names(col_avg)==F)]
    all_avg <- c(col_avg, rem_avg)
    # order metabolites
    D <- D[match(names(all_avg), rowData(D)$name), ]
    X <- assay(D)
    # subtract the ratio avg 
    X <- apply(X, 2, FUN=function(x) x - all_avg)
  }
  
  # add status information
  funargs <- maplet:::mti_funargs()
  uuid = uuid::UUIDgenerate()
  
  metadata(D)$results <- maplet:::mti_add_to_list(
    metadata(D)$results, 
    list(fun = funargs$fun, args = funargs$args, logtxt = maplet:::mti_logmsg(""), 
         logshort = "", uuid = uuid, output = NULL, output2 = NULL), 
    oname = paste(paste(funargs$fun, collapse = "_"), uuid, sep = ".")
  )
  
  assay(D) <- X
  return(D)
}


# For anchor normalization of dataset3. It is called by function normalize_by_deviation 
# @param D: SE (dataset3)
# @param refD: SE (dataset2)
# @return ratio_mat: a matrix with ratios of anchor samples (dataset3/dataset2)
get_anchor_ratios <- function(D, refD){
  # fetch Patient_IDs of anchor samples
  anchors <- refD %>% colData() %>% as.tibble() %>% filter(anchor=="Yes") %>% select(Sample_ID) %>% unlist()
  # check for common metabolites
  common_mets <- intersect(rowData(D)$name, rowData(refD)$name)
  # check for order of metabolites
  m1 <- match(common_mets, rowData(D)$name)
  m2 <- match(common_mets, rowData(refD)$name)
  # compute ratios
  ratios <- lapply(1:length(anchors), FUN=function(i){
    x <- assay(D)[m1, which(colData(D)$Sample_ID %in% anchors[i])]
    y <- assay(refD)[m2, which(colData(refD)$Sample_ID %in% anchors[i])]
    r <- x - y
    return(r)
  })
  ratio_mat <- do.call(rbind, ratios)
  rownames(ratio_mat) <- anchors
  colnames(ratio_mat) <- common_mets
  return(ratio_mat)
}


# For processing GIM dataset123. It is called by function GIM_preprocessing
# It calls function GIM_preprocessing_default for normalization of the two SEs one by one
# @param D: (dataset3) SE 
# @param refD: (dataset2) a second SE to be used for anchor normalization of D
# @param anchor_method: optional, method to be used for anchor normalization of D given refD
# @param exclude.pregnant: optional flag to exclude pregnant samples from data
# @returns D: a processed combined SE
GIM_preprocessing_dataset123 <- function(Ds, exclude.pregnant=T){
  
  D1 <- Ds$D1
  D2 <- Ds$D2
  D3 <- Ds$D3
  
  # exclude pregnant except if it is an anchor
  if (exclude.pregnant) {
    D1 %<>% mt_modify_filter_samples(filter = (pregnancy=="No"))
    D2 %<>% mt_modify_filter_samples(filter = (pregnancy=="No" | anchor=='Yes'))
    D3 %<>% mt_modify_filter_samples(filter = (pregnancy=="No" | anchor=="Yes"))
  }
  
  # Preprocess separately 
  D3 %<>% GIM_preprocessing_default (D2, to_impute = F) %>%
    # anchors no longer needed
    mt_modify_filter_samples(filter = (anchor=="No"))
  # pregnant anchor no longer needed 
  D2 %<>% mt_modify_filter_samples(filter = (pregnancy=="No")) %>% 
    GIM_preprocessing_default (to_impute = F)
  D1 %<>% GIM_preprocessing_default (to_impute = F)
  
  # check for common metabolites
  common_mets <- Reduce(f=intersect, x=list(rowData(D1)$name, rowData(D2)$name, rowData(D3)$name))
  
  # check for order of metabolites
  m1 <- match(common_mets, rowData(D1)$name)
  m2 <- match(common_mets, rowData(D2)$name)
  m3 <- match(common_mets, rowData(D3)$name)
  
  # finally merge
  D <- SummarizedExperiment(
    assay    = cbind(assay(D1)[m1, ], assay(D2)[m2, ], assay(D3)[m3, ]),
    colData  = plyr::rbind.fill(data.frame(colData(D1)), data.frame(colData(D2)), data.frame(colData(D3))),
    rowData  = rowData(D1)[m1, c('name', 'HMDB')],
    metadata = list()
  )
  
  D %<>% mt_pre_trans_exp() %>% 
    mt_pre_batch_median(batch_col='batch') %>%
    mt_pre_trans_log() %>% mt_pre_impute_knn()
  
  # pathway annotations, if there is an HMDB column
  if ("HMDB" %in% (D %>% rowData() %>% colnames())) {
    D %<>%  
      # HMDB to KEGG
      mt_anno_hmdb_to_kegg(
        in_col = 'HMDB',
        out_col = 'KEGG_mapped'
      ) %>%
      # pull pathway annotations
      mt_anno_pathways_hmdb(
        in_col = "HMDB", 
        out_col = "kegg_db", 
        pwdb_name = "KEGG") %>% 
      mt_anno_pathways_remove_redundant(feat_col = "HMDB", pw_col = "kegg_db") %>%
      mt_write_pathways(pw_col = "kegg_db", file="results/pw_annos_kegg.xlsx")
  }
  
  # final dataset info
  D %<>% mt_reporting_data()
  
  return(D)
}


# For outcome datatype conversion. It is called in the main script 
# @param D: SE
# @param outcomes: dataframe tiwth outcome names and data types
# @return D: with right datatype of the outcomes
GIM_outcome_type_conversion <- function(D, outcomes){
  # make sure numeric outcomes are numeric
  for (outcome in (outcomes %>% filter(outcomeMode=="numeric") %>% .$outcome)) {
    D %<>% mt_anno_mutate(anno_type = "samples", col_name = outcome, term = as.numeric(!!sym(outcome)))
  }
  for (outcome in (outcomes %>% filter(outcomeType=="binary") %>% .$outcome)) {
    D %<>% mt_anno_mutate(anno_type = "samples", col_name = outcome, term = as.factor(!!sym(outcome)))
  }
  
  D
}
