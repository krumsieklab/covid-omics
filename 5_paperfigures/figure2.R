
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(scales)

# load the results of case-control analysis 
load(glue::glue("{root}/2a_casecontrol/results/GIM_Case-Control_metabo.rds"))
Ds = list(metabo = D)
load(glue::glue("{root}/2a_casecontrol/results/GIM_Case-Control_proteo.rds"))
Ds$proteo = D
rm(D)

# reverse log transform
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# results of association analysis 
res = lapply(Ds, function(D){
  re = maplet::mtm_res_get_entries(D, "stats")[[1]]$output$table
  re$name = make.names(rowData(D)$name)
  re
})


# volcano plots -----------------------------------------------------------
ggs_volcano <-
lapply(res, function(df){
  df$p.adj[df$p.adj>0.05] = NA
  df$name[is.na(df$p.adj)] = NA 
  
  gg = ggplot(df, aes(x = fc, y = p.value, 
                      color =  c("down", "up")[(p.adj < 0.05) * (fc>0)+1]  )) + 
    geom_point() + #aes(alpha = -log10(p.value))  
    scale_y_continuous(
      trans = reverselog_trans(10), 
      breaks = scales::trans_breaks("log10", function(x) 10^x), 
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    labs(color = "", x = "fold change") +
    geom_text_repel(aes(x = fc, y = p.value, label = name ))+
    theme_minimal() +
    scale_color_manual(values =c("#E69F00", "#56B4E9"), na.value = "gray70")
  xlims = ggplot_build(gg)$layout$panel_params[[1]]$x.range 
  gg + 
    annotate("text", x = xlims[1], y = 1, label ="atop(bold(`Contol high`))", 
             hjust = 0, vjust = 1, color = "#E69F00", parse = TRUE)+
    annotate("text", x = xlims[2], y = 1, label = "atop(bold(`COVID19 high`))", 
             hjust = 1, vjust = 1, color = "#56B4E9", parse = TRUE) +
    scale_alpha_continuous(range = c(0.3,1)) + 
    guides(alpha = "none")
})


# boxplots of top hits ----------------------------------------------------
ggs_boxplots<-
structure(names(res), names = names(res)) %>% lapply(function(idata){
  # idata = names(res)[2]
  mets = 
    res[[idata]] %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(p = unique(p.adj)) %>% # p.adj, fc
    dplyr::slice_max(-abs(p), n=40) %>% 
    { as.character( .$name ) }  
  mets = mets[c(1:2,3)]
  
  x = assay(Ds[[idata]]) %>% t
  colnames(x) = make.names(rowData(Ds[[idata]])$name)
  
  df = data.frame(x[,mets], check.names = F,
                  colData(Ds[[idata]])[, c("Time", "Patient_ID", "Status")])
  
  # scale the data based on control
  for(i in mets){
    sm = unlist( attributes(scale(df[[i]][df$Status == "Control"]))[c("scaled:center","scaled:scale")])
    df[[i]] = scale(df[[i]], center = sm[1], scale = sm[2]) %>% as.numeric()
  }
  rm(sm, i)
  
  # melt the data 
  df = reshape2::melt(df[,setdiff(colnames(df),"Time")], id = c("Patient_ID", "Status"))
  
  # summarize multiple data points of each patients
  mdf <-  
    df %>% 
    dplyr::group_by(Patient_ID, variable, Status) %>% 
    dplyr::summarise(value = mean(value), Status = unique(Status))
  
  mdf$Status = factor(mdf$Status, labels = c("Control", "COVID19"))
  
  gg_met <-
    ggplot(mdf, aes(y=log10(value), x = Status, color = Status, fill = Status)) +
    #geom_jitter(width = 0.2, alpha = 0.20, pch =21, color = "white", size = 2) + 
    geom_boxplot(fill = NA, outlier.color = NA)+#, notch = T) +
    facet_wrap(~variable, nrow = 1, scales = "free_x", switch = "y") +
    theme_minimal() +
    ggsci::scale_color_aaas() +
    theme(strip.background = element_rect(fill = "wheat", colour = "black"), 
          panel.background = element_rect(colour = "black")) +
    ylim(c(-2.5,NA)) + labs(x="", y = "centered expression")
  #scale_y_continuous(limits = c(NA, 5), oob = scales::squish)
  
  (gg_met = 
      gg_met + 
      theme(legend.position = c(0.87, 0.2),
            axis.text.x = element_text(angle = -30, hjust = 0.5, vjust = 0)) + 
      labs(color = "", y = expression(paste("Normalized ",log[2], " fold change"))))
  
  # svg(filename = "5_paperfigures/boxplot_metabo.svg",width=4, height=3.6)
  # gg_met
  # dev.off()
  gg_met
})


# pathway barplots --------------------------------------------------------
library(SummarizedExperiment)
library(readxl)
library(dplyr)
library(magrittr)
library(ggplot2)

# Figure 2 (result overview), pathway plots 
# list of curated pathways
pws <- read_excel(glue::glue("{root}/DATA/KEGGpw_Hierarchy.xlsx"), sheet="ManuallyCuratedList")
# remove disease pathways
pws = setdiff(pws$Name,pws[grep("disease", pws$Group ),]$Name  %>% {.[-grep("COVID",.)]})
pws = setdiff(pws, grep("Chagas|infection|cancer|arthritis|sclerosis|Influenza|Chemokine|cytokine", pws,value = T))

# files where pathway results in 
resfiles = list( 
  metabo = glue::glue("{root}/2a_casecontrol/results/GIM_Case-Control_metabo_pathway_stat.xlsx"),
  proteo = glue::glue("{root}/2a_casecontrol/results/GIM_Case-Control_proteo_pathway_stat.xlsx")
)

ggs_pwbarplots<- resfiles %>% lapply(function(input_file){
  
  # authors: RB, JK
  label = ifelse( length(grep("metabo", input_file))>0, "metabolites", "proteins" )
  print(label)
  xlabel=sprintf('# of %s', label) # x-axis label
  grp_name1='Down in COVID-19' # statistics < 0
  grp_name2='Up in COVID-19' # statistics > 0
  path_order='max2min' # order in the barplot
  topN = 15
  add_pw_size = FALSE
  legend_title = ""
  filter_to = if( label == "proteins" ) pws else NULL
  
  # load both sheets
  tmp <- read_xlsx(input_file, sheet='AggregatedPathways') %>%
    select(pathway_name, label)
  data_mat <- read_xlsx(input_file, sheet='IndividualResults') %>%
    left_join(tmp, by='pathway_name') %>% unique()  
  
  # filter to non-disease pathways 
  if (!is.null(filter_to)) {
    data_mat = data_mat %>% filter(pathway_name  %in% filter_to)
  }
  
  # get data into the right shape
  plot_mat <- dplyr::left_join(data_mat, data_mat %>% dplyr::group_by(label) %>%
                                 dplyr::count(name) %>% dplyr::group_by(label) %>% 
                                 dplyr::count(label), by='label') %>% 
    dplyr::rename(path_wt=n) %>%
    dplyr::left_join(data_mat %>% group_by(label) %>% 
                       dplyr::count(statistic<0) %>% 
                       filter(`statistic < 0` ==TRUE) %>% 
                       select(- `statistic < 0`), by='label') %>% 
    dplyr::rename(grp1=n) %>%
    dplyr::left_join(data_mat %>% group_by(label) %>% 
                       dplyr::count(statistic>0) %>% 
                       filter(`statistic > 0` ==TRUE) %>% 
                       select(- `statistic > 0`), by='label') %>% 
    dplyr::rename(grp2=n) %>% 
    select(label, path_wt, grp1, grp2) %>% unique()
  
  # filter top N?
  if(is.finite(topN) && (topN <nrow(plot_mat))) {
    plot_mat <- plot_mat[order(plot_mat$path_wt, decreasing = T), ]
    plot_mat <- plot_mat[1:topN, ]
  }
  
  # order up/down or down/up
  if(path_order=='max2min'){
    # max at the top
    path_order <- plot_mat$label[order(plot_mat$path_wt, decreasing = F)]  
  } else if(path_order=='min2max'){
    # max at the bottom
    path_order <- plot_mat$label[order(plot_mat$path_wt, decreasing = T)]
  } else {
    # bad input
    print('Unrecognized order... Plotting max2min!')
    path_order <- plot_mat$label[order(plot_mat$path_wt, decreasing=F)]  
  }
  
  # prepare matrix for ggplot
  plot_mat <- reshape2::melt(plot_mat, id=c('label')) %>% 
    filter(variable!='path_wt')
  plot_mat$grp <- sub('..*_', '', plot_mat$variable)
  plot_mat$label <- factor(plot_mat$label, levels=path_order)
  
  if(!add_pw_size){
    levels(plot_mat$label) = gsub(" \\[.*\\]","", levels(plot_mat$label)) 
  }
  
  # interesting pathways for proteins
  if(label == "proteins"){
    plot_mat = plot_mat[-seq(3),]
  }
  
  # ggplot
  p <- ggplot(plot_mat, aes(fill=variable, x=label, y=value)) +
    geom_bar(position='stack', stat='identity') + theme_bw()+
    xlab('') + ylab(xlabel) + ggtitle('')+
    theme(text=element_text(size=12))+
    coord_flip() +
    scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                      name= legend_title,
                      breaks=c("grp1", "grp2"),
                      labels=c(grp_name1, grp_name2))
  return(p)
})

save(file = "figure2.fig", ggs_boxplots, ggs_pwbarplots, ggs_volcano)

# # combine volcano and boxplots --------------------------------------------
# library(patchwork)
# cw<-
#   cowplot::plot_grid(
#     cowplot::plot_grid(ggs_volcano$metabo + guides(color = F), 
#                        ggs_volcano$proteo + guides(color = F), ncol = 1),
#     cowplot::plot_grid(ggs_boxplots$metabo, 
#                        ggs_boxplots$proteo, ncol = 1), 
#     rel_widths = c(3,2), ncol = 1 )
# 
# 
# 
# # svg(filename = "5_paperfigures/combined.svg",width=10.5, height=7)
# # cw
# # dev.off()
# 
# cw2<-
# cowplot::plot_grid(ggs_pwbarplots$metabo + guides(color = F), 
#                    ggs_pwbarplots$proteo + guides(color = F), ncol = 1)
# 
