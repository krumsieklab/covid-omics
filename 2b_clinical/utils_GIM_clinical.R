
# workhorse function to run whole analysis pipeline
f_WORKHORSE <- function(dataset, #COVID19_datapath, 
                        outcomes,
                        datestamp, 
                        case_only,
                        exclude.intubated, formula_conf, p.adj.cut, 
                        do.boxplots, pboxcut,
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
  files_created = list.files(path = "results", recursive = T)
  
  # intubation
  intubadd <- if (exclude.intubated)"" else "_withintub"
  
  #'
  #' Steps
  #' 1. loading & preprocessing
  #' 2. global analysis
  #' 3. differential analysis
  #' 
  
  if(write_results){
    # file name suffix for excluded samples?
    intubadd <- if (exclude.intubated)""else"_withintub"
    # base path for output files
    outbase <- sprintf("results/GIM_clinicals_%s%s", dataset, intubadd) # , datestamp ageadd, ratioadd
    # make sure output directory exists
    dir.create("results", showWarnings = F, recursive = T)
  }
  
  #### re-usable analysis function ----
  p_analyze <- function(D, outcome_info, compname, write_results=T ) {
    
    # lm for binary and numeric outcomes
    if (outcome_info$outcomeType %in% c("binary", "numeric")) {
      # lm
      D %<>%
        mt_stats_univ_lm(
          formula = as.formula(glue("~{outcome_info$outcome} + {formula_conf} + (1|Patient_ID)")), 
          stat_name         = compname
        ) 
    } else if (outcome_info$outcomeType ==  "composite") {
      # lm
      D %<>%
        mt_stats_univ_cs(
          in_col = outcome_info$outcome,
          id_col = "Patient_ID", 
          stat_name = compname
        ) 
    }
    else {
      stop("only composite, binary and numeric outcomes implemented so far")
    }
    
    # pval QQ plot
    D %<>% mt_plots_pval_qq(stat_name = compname)
    
    # fold changes, multiple testing correction
    D <- D %>%
      {if (outcome_info$outcomeType=="binary"){mt_post_fold_change(., stat_name = compname)}else{.}} %>%
      mt_post_multtest(stat_name = compname, method = "BH") %>% {.}
        
    # stats
    D %<>% 
      mt_reporting_stats(stat_name = compname, stat_filter =  p.adj < !!enquo(p.adj.cut))
    
    # binary: Volcano plot and boxplots
    # numeric: ordered p-value plot and scatter 
    if (outcome_info$outcomeType == "binary") {
      D %<>% 
        mt_plots_volcano(stat_name     = compname,
                         feat_filter = p.adj < !!enquo(p.adj.cut),
                         colour       =  p.adj < !!enquo(p.adj.cut))
      
      if (do.boxplots) {
        D %<>%
          mt_plots_box_scatter(stat_name           = compname,
                                   plot_type = 'box',
                                   x                  = !!sym(outcome_info$outcome),
                                   fill               = !!sym(outcome_info$outcome),
                                   feat_sort         = p.value,
                                   feat_filter       = p.adj<!!enquo(pboxcut),
                                   annotation         = "{sprintf('P-value: %.1e', p.value)}\nadjusted: {sprintf('%.1e', p.adj)}",
                                   #rows               = 2,
                                   #cols               = 3,
                                   restrict_to_used_samples = T)
      }
      # stats barplots (if kegg_db given)
      if ("kegg_db" %in% colnames(rowData(D))) {
        D %<>% 
          mt_plots_stats_pathway_bar( 
            stat_list = compname,
            group_col = "kegg_db",
            feat_filter = p.adj<!!enquo(p.adj.cut), 
            y_scale = "count",
            sort_by_y = F,
            assoc_sign = "fc",
            outfile = if(write_results) sprintf("%s_%s_%s_%s.xlsx",outbase, dataset, compname, 'pathway_stat') else NULL)
      }
      
      
    } else if (outcome_info$outcomeType == "numeric" | outcome_info$outcomeType == "composite") {
      D %<>%
        mt_plots_volcano(stat_name     = compname,
                         x = statistic,
                         feat_filter = p.adj < !!enquo(p.adj.cut),
                         colour       =  p.adj < !!enquo(p.adj.cut))
      if (do.boxplots) {
        D %<>% 
          mt_plots_box_scatter(stat_name           = compname,
                                   plot_type = 'scatter',
                                   x                  = !!sym(outcome_info$outcome),
                                   feat_sort         = p.value,
                                   feat_filter       = p.adj<!!enquo(pboxcut),
                                   annotation         = "{sprintf('P-value: %.1e', p.value)}\nadjusted: {sprintf('%.1e', p.adj)}",
                                   #rows               = 2,
                                   #cols               = 3,
                                   restrict_to_used_samples = T) 
      }
      # stats barplots (if kegg_db given)
      if ("kegg_db" %in% colnames(rowData(D))) {
        D %<>% 
          mt_plots_stats_pathway_bar(
            stat_list = compname,
            group_col =  "kegg_db",
            feat_filter = p.adj<!!enquo(p.adj.cut), 
            y_scale = "count",
            assoc_sign = "statistic",
            sort_by_y = F, 
            outfile = if(write_results) sprintf("%s_%s_%s_%s.xlsx",outbase, dataset, compname, 'pathway_stat') else NULL)
      }
      
      
    }
    
    
    # return
    D
  }
  
  #### loading & preprocessing ----
  if(is.null(D)){
    D <- 
      # load data
      GIM_data_loader(exclude.intubated=exclude.intubated, dataset=dataset) %>%
      # preprocess
      GIM_preprocessing() %>% 
      # convert outcomes into numeric / factors 
      GIM_outcome_type_conversion (outcomes)  %>%
      {.}
  }
  
  # filter down outcomes to those that are actually present in the dataset
  outcomes <- outcomes %>% dplyr::filter(outcome %in% colnames(colData(D)))
  
  # if case only
  if(case_only){
    D %<>% mt_modify_filter_samples(Status=="COVID")
  }
  
  #### global analysis ----
  D %<>%
    # header
    mt_reporting_heading("Global analysis") %>%
    
    # heatmap
    mt_plots_heatmap(scale_data = T, annotation_col = "Group", clustering_method = "ward.D2", fontsize = 5, 
                      cluster_rows = T, color=gplots::bluered(101), sym_zero = T) %>%
    
    # PCA & UMAP
    mt_plots_pca(scale_data = T, color = Group) %>%
    mt_plots_umap(scale_data = T, color = Group) %>%
    
    # for coding convenience (empty statement to terminate %>% pipe)
    {.}
  
  # loop over outcomes
  for (i in 1:nrow(outcomes)) {
    D %<>% 
      mt_reporting_heading(heading = sprintf("%s", outcomes$outcome[[i]])) %>%
      p_analyze(outcome_info = outcomes[i,], 
                compname = outcomes$outcome[[i]], 
                write_results = write_results) 
  }
  
  # create equalizer plots
  ggs_equ<- lapply( structure(outcomes$outcome, names = outcomes$outcome), function(outcome){
    get_equalizer(D=D,comparison = outcome, pwvar = "kegg_db",  filename = if(write_results)  
      sprintf("results/Equalizer_%s_%s",dataset, outcome) else NULL)
  })
  
  if(write_results){
    # export all, write out all statistics
    D %>% mt_write_stats(file=sprintf("%s_allstats.xlsx", outbase), feat_col = 'name', sort_by_p=T)
  }
  
  #### add overview heatmap ----
  # gather names of comparisons first
  comps <- D %>% 
    maplet::mtm_res_get_entries("stats") %>%
    map("output") %>% map("name") %>% unlist()
  
  D %<>% mt_reporting_heading(heading = "Overview heatmaps") 
  # first version, leave out sex,age
  D %<>%
    mt_plots_multstats_heatmap(
      stat_list = setdiff(comps, c("age","sex")),
      color_signif = T, 
      cutoff = 0.05, 
      signif_mark = "", 
      cluster_rows = T, 
      cluster_cols = T,
      filter_signif = F
    )
  
  if(output_html){
    #### output to HTML and save ----
    D %>%
      mt_reporting_html(sprintf("%s.html", outbase), use_plotly = F, number_sections = T) %>%
      #mt_stripresults(strip = "plots") %>%
      mt_clean_remove_results(remove = "plots") %>%
      mt_write_se_rds(file = sprintf("%s.rds", outbase))
  }
  
  #### assemble all results into one view ----
  # will produce PDF file in results/ folder
  library(stringr)
  library(RColorBrewer)

  if(write_results){
    graphics.off(); 
    pdf(sprintf("%s_overview.pdf", outbase), width=4, height=4)
  }
  
  # extractor helper function
  get.all <- function(D){
    allstats <- D %>% 
      maplet::mtm_res_get_entries("stats") %>%
      map("output") %>% map("table")
    names(allstats) <- D %>% 
      maplet::mtm_res_get_entries("stats") %>%
      map("output") %>% map("name")
    allstats
  }
  
  # converts a matrix of lists into a matrix of atomics
  # terrible code... no idea how else to do this.... 
  unlist.mat <- function(m) {
    m2 <- matrix(NA, nrow=nrow(m), ncol=ncol(m))
    for (i in 1:nrow(m2)) {
      for (j in 1:ncol(m2)) {
        m2[i,j] <- m[[i,j]]
      }
    }
    rownames(m2) <- rownames(m)
    colnames(m2) <- colnames(m)
    m2
  }
  
  ## extract all tests 
  allres <- get.all(D)
  # generate color scale by max hits
  m1 <- allres %>% map(~.x$p.adj < p.adj.cut) %>% map(sum) %>% unlist() %>% max()
  maxhits <- max(m1)
  # pheatmap parameters
  breaksList = seq(0, maxhits, by = 1)
  color = colorRampPalette(rev(brewer.pal(n=7,name = "Spectral")))(length(breaksList))
  
  ## main run
  # store objects in matrix form
  mat <- as.list(rep(NA, nrow(outcomes)))
  for (i in 1:length(allres)) {
    # find outcome and day
    s <- names(allres)[i]
    #day = unlist(strsplit(s, "  "))[2]
    # write into matrix
    mat[[which(outcomes$outcome==s)]] <- allres[[i]]
  }
  # label matrix
  names(mat) <- outcomes$outcome
  
  # get number of significant p-values
  mat.hits <- sapply(mat, function(tab){sum(tab$p.adj<p.adj.cut)})
  
  # visualize
  phemp<-
  pheatmap::pheatmap(
    as.data.frame(mat.hits),
    cluster_cols = F, cluster_rows = F, 
    display_numbers = T, number_color = "#000000", fontsize_number = 12, angle_col = 0,
    main = sprintf("p.adj< %f", p.adj.cut),
    color = color, breaks = breaksList,
    number_format = "%d"
  )
  
  if(write_results) dev.off()
  
  # check important warnings do.call(warnings, args = list(), envir = globalenv())
  important_warnings =  attributes(warnings())$names %>% 
    unique %>% {.[-grep(.,pattern = "deprecated|non-finite")]}  #%>% cat(sep = "\n")
  
  files_created = setdiff( list.files(path = "results", recursive = T), files_created )
  
  list(D = D, ggs_equalizer = ggs_equ, results_heatmap = phemp,important_warnings = important_warnings, call = mcall, files_created = files_created)
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

