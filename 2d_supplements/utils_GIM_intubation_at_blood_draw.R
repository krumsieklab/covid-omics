
# workhorse function to run whole analysis pipeline
f_WORKHORSE_int_at_bld_drw <- function(dataset, #COVID19_datapath, 
                        datestamp, 
                        formula_conf, p.adj.cut, 
                        write_results = T, 
                        output_html = T){
  
  if(is.list(dataset)){
    D = dataset$D
    dataset = dataset$dataset_name
  }else{
    D = NULL
  }
  
  mcall = match.call()
  print(dataset)
  # get the list of already created files
  files_created = list.files(path = "results/", recursive = T)
  
  #'
  #' Steps
  #' 1. loading & preprocessing
  #' 2. global analysis
  #' 3. differential analysis
  #' 
  
  if(write_results){
    # sanity check
    # stopifnot(dataset=='dataset1' | dataset=="dataset23" | dataset=="dataset1_proteo") # ugly, but necessary due to copy/pasted code
    # construct correction string
    # corradd <- paste0("_corr_", formula_conf %>% strsplit("\\w*\\+\\w*") %>% .[[1]] %>% trimws() %>% paste0(collapse = "_"))
    # file name suffix for excluded samples?
    intubadd <- "intubation_at_blood_draw"
    # base path for output files
    outbase <- sprintf("results/GIM_%s%s", dataset, intubadd) #, datestamp
    # make sure output directory exists
    dir.create("results", showWarnings = F, recursive = T)
  }
  

  #### global analysis ----
  D %<>%
    # header
    mt_reporting_heading("Global analysis") %>%
    # heatmap
    mt_plots_heatmap(scale_data = T, 
                      annotation_col = "Group", 
                      clustering_method = "ward.D2", 
                      fontsize = 5, 
                      cluster_rows = T, 
                      color=gplots::bluered(101), 
                      sym_zero = T) %>%
    
    # PCA & UMAP
    mt_plots_pca(scale_data = T, color = Group) %>%
    mt_plots_pca(scale_data = T, color = batch) %>%
    mt_plots_umap(scale_data = T, color = batch) %>%
    mt_plots_umap(scale_data = T, color = Group) %>%
    {.}
  
  #### differential analysis, all COVID samples pooled ----
  compname = "int_at_bld_drw"
  D %<>%
    # header
    mt_reporting_heading("Intubation at blood draw") %>%
    mt_anno_reorder_factor("int_at_bld_drw", c('No', 'Yes')) %>%
    # analyze
    mt_stats_univ_lm(
      formula      = as.formula(glue("~ int_at_bld_drw + {formula_conf} + (1|Patient_ID)")),
      stat_name    = compname
    ) %>% 
    mt_post_fold_change(stat_name = compname) %>%
    mt_post_multtest(stat_name = compname, method = "BH") %>%
    mt_reporting_stats(stat_name     = compname, stat_filter = p.adj < !!enquo(p.adj.cut )) %>%
    mt_plots_volcano(stat_name     = compname,
                     feat_filter = p.adj < !!enquo(p.adj.cut),
                     colour       =  p.adj < !!enquo(p.adj.cut) ) %>%
    # add boxplot
    mt_plots_box_scatter(stat_name           = compname,
                             plot_type = 'box',
                             x                  = Status,
                             fill               = Status,
                             feat_sort         = p.value,
                             feat_filter       = p.adj<p.adj.cut,
                             annotation         = "{sprintf('P-value: %.1e', p.value)}\nadjusted: {sprintf('%.1e', p.adj)}",
                             # rows               = 2,
                             # cols               = 2,
                             restrict_to_used_samples = T)
  
  
  if(write_results){
    # export all
    D %<>% 
      mt_write_stats(file = sprintf("%s.xlsx", outbase), sort_by_p = T, feat_col= 'name'
                     #, output.dir = T
                     ) %>%
      {.}
  }
  
  if(output_html){
    #### output to HTML and save ----
    D %>% 
      mt_reporting_html(sprintf("%s.html", outbase), use_plotly = F, number_sections = T) %>%
      #mt_stripresults(strip = "plots") %>%
      mt_clean_remove_results(remove = "plots") %>%
      mt_write_se_rds(file = sprintf("%s.rds", outbase))
  }
  
 
  # check important warnings do.call(warnings, args = list(), envir = globalenv())
  important_warnings =  attributes(warnings())$names %>% 
    unique %>% {.[-grep(.,pattern = "deprecated|non-finite")]}  #%>% cat(sep = "\n")
  
  files_created = setdiff( list.files(path = "results", recursive = T), files_created )
  
  list(D = D, important_warnings = important_warnings, call = mcall, files_created = files_created)
}

