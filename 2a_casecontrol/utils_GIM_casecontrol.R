
# workhorse function to run whole analysis pipeline
f_WORKHORSE <- function(dataset, #COVID19_datapath, 
                        datestamp, exclude.intubated, formula_conf, p.adj.cut, 
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
    intubadd <- if (exclude.intubated)""else"_withintub"
    # base path for output files
    outbase <- sprintf("results/GIM_Case-Control_%s", dataset) #, datestamp
    # make sure output directory exists
    dir.create("results", showWarnings = F, recursive = T)
  }
  
  if(is.null(D)){
    #### loading & preprocessing ----
    D <- 
      # load data
      GIM_data_loader(exclude.intubated=exclude.intubated, dataset=dataset) %>%
      # preprocess
      GIM_preprocessing() %>% 
      {.}
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
  compname = "casecontrol"
  D %<>%
    # header
    mt_reporting_heading("Case/control") %>%
    mt_anno_reorder_factor("Status", c('Control', 'COVID')) %>%
    # analyze
    mt_stats_univ_lm(
      formula      = as.formula(glue("~ Status + {formula_conf} + (1|Patient_ID)")),
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
  
  # stats barplots (if kegg_db given)
  if ("kegg_db" %in% colnames(rowData(D))){
    D %<>% mt_plots_stats_pathway_bar(stat_list = compname,
                                 group_col = "kegg_db",
                                 feat_filter = p.adj<p.adj.cut, 
                                 y_scale = "count",
                                 sort = F,
                                 assoc_sign_col = "fc", 
                                 outfile=if(write_results) sprintf("%s_%s.xlsx", outbase, 'pathway_stat') else NULL)
  }else{
    print("stats barplots is not produced  because kegg_db is NOT given!")
  }
  
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
  
  gg_equ <-
    get_equalizer(D=D, comparison = "casecontrol", pwvar = "kegg_db", filename = if(write_results) 
      sprintf("results/Equalizer_casecontrol_%s_%s", dataset, datestamp) else NULL)
  
  # check important warnings do.call(warnings, args = list(), envir = globalenv())
  important_warnings =  attributes(warnings())$names %>% 
    unique %>% {.[-grep(.,pattern = "deprecated|non-finite")]}  #%>% cat(sep = "\n")
  
  files_created = setdiff( list.files(path = "results", recursive = T), files_created )
  
  list(D = D, gg_equalizer = gg_equ, important_warnings = important_warnings, call = mcall, files_created = files_created)
}
# equalizer plot
get_equalizer <- function(D, comparison, pwvar, filename = NULL, pathway_type='kegg'){
  if(pathway_type=='kegg'){
    pp <- metadata(D)$pathways[[pwvar]] %>% 
      dplyr::select(ID, pathway_name)
    
    # get pathway annotations
    rd <- rowData(D) %>% as.data.frame %>%
      dplyr::mutate(var = rownames(D))
    
    # add pathway annotations to results
    res0 <- 
      maplet::mtm_get_stat_by_name(D, comparison) %>%
      dplyr::left_join(rd, by="var") %>%
      dplyr::mutate(x = sign(statistic)*log10(p.adj)) 
    res0$kegg_db[which(sapply(res0$kegg_db, is.null))] <- "unmapped"
    
    res <- res0 %>%
      tidyr::unnest(!!sym(pwvar)) %>%
      dplyr::left_join(pp , by=c("kegg_db"="ID"))
    res$pathway_name[res$kegg_db=="unmapped"] <- "Unmapped Metabolites"
    
  } else if(pathway_type=='metabolon'){
    # get pathway annotations
    rd <- rowData(D) %>% as.data.frame %>%
      dplyr::mutate(var = rownames(D), pathway_name=SUB_PATHWAY)
    
    # add pathway annotations to results
    res0 <- 
      maplet::mtm_get_stat_by_name(D, comparison) %>%
      dplyr::left_join(rd, by="var") %>%
      dplyr::mutate(x = sign(statistic)*log10(p.adj)) 
    
    res0$pathway_name[which(sapply(res0$pathway_name, is.null))] <- "unmapped"
    
    res <- res0 %>%
      tidyr::unnest(!!sym(pwvar))
    
    res$pathway_name[res$pathway_name=="unmapped"] <- "Unmapped Metabolites"
    
  } else if (pathway_type=='lipids') {
    # get pathway annotations
    rd <- rowData(D) %>% as.data.frame %>%
      dplyr::mutate(var = rownames(D), pathway_name=Class)
    
    # add pathway annotations to results
    res0 <- 
      maplet::mtm_get_stat_by_name(D, comparison) %>%
      dplyr::left_join(rd, by="var") %>%
      dplyr::mutate(x = sign(statistic)*log10(p.adj)) 
    
    res0$pathway_name[which(sapply(res0$pathway_name, is.null))] <- "unmapped"
    
    res <- res0 %>%
      tidyr::unnest(!!sym(pwvar))
    
    res$pathway_name[res$pathway_name=="unmapped"] <- "Unmapped Metabolites"
  }
  # colors
  clrs <- c("#9494FF","red")
  
  # compute multiple testing correction line
  sel <- res %>% 
    dplyr::filter(p.adj < p.adj.cut) 
  if(nrow(sel)>0){
    xfine <- sel %>% .$p.adj %>% max(., na.rm = T)
  } else {
    xfine <- Inf
  }
  
  # create equalizer plots
  
  # select metabolites in pw y
  df <- res %>%
    # dplyr::filter(!!sym(pwvar)==y) %>%
    # dplyr::mutate(x = fc) %>%
    dplyr::mutate(x = -sign(statistic)*log10(p.adj)) %>%
    dplyr::mutate(col = ifelse(p.adj<p.adj.cut, "significant", "non significant"))
  # x axis limits
  a = max(abs(df$x))+0.3
  
  eqp <- ggplot(df, aes(x = x, y = name, fill=col)) +
    geom_vline(xintercept = 0, color ="gray") +
    (if(!is.infinite(xfine)){geom_vline(xintercept = c(-log10(xfine),log10(xfine)), color="red", alpha=0.4)}) +
    (if(!is.infinite(xfine)){ggtitle(sprintf("Differential Metabolites at alpha %.2f",p.adj.cut))}else{ggtitle(sprintf("No significant results at alpha %.2f", p.adj.cut))}) +
    geom_point(pch = 22, color = "black", size = 3) +
    facet_grid(as.formula(sprintf("pathway_name~.")), scales = "free_y", space = "free_y") +
    theme(strip.background =element_rect(fill=NA),
          strip.text = element_text(colour = 'black', face = "bold"),
          strip.text.y = element_text(angle = 0, hjust = 0),
          panel.grid.major.y = element_line(color ="gray"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill=NA, color ="black")) +
    ylab("") +
    xlab("Directed log10(p.adj)") +
    scale_x_continuous(limits = c(-a,a)) 
  
  if(!is.null(filename)){
    pdf(sprintf("%s.pdf",filename), height = 200, width = 10)
    print(eqp)
    dev.off()
  }
  eqp
  
}
