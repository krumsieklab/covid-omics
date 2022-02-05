# GGM given the preprocessed data
# this script produced figure 3A and 3B networks 

# load the preprocessed data 
# load(paste0(root,"/1_preprocessing/data_preprocessed_for_GGM.Rdata"))
load(path_preprocessed_GGM_data)

# set working directory to location of source code
setwd(glue::glue("{root}/2c_networks"))

# GGM ---------------------------------------------------------------------
library(gplots)
library(GeneNet)
colors <- gplots::bluered(100)

datestamp <- format(Sys.time(),"%Y%m%d")


# molecule types
mol_types <- lapply(Ds, FUN=function(d) {
  res <- d %>% rowData() %>% as_tibble() %>% select(name)
  res %<>% mutate(unique_name=make.names(name, unique=T))
  return(res)
})

# samples common in 2 omics
com_samples <- lapply(Ds, FUN=function(D) D %>% colData() %>% 
                        as_tibble() %>% select(Sample_ID) %>% 
                        unlist()) %>% Reduce(intersect, .)

# join the molecules
data_mat <- lapply(Ds, FUN=function(D) {
  res <- D %>% assay() %>% t() %>% as_tibble %>% filter(colData(D)$Sample_ID%in%com_samples) %>% data.frame()
  colnames(res) <- D %>% rowData() %>% as_tibble() %>% select(name) %>% unlist() %>% make.names(unique = T)
  
  patIDs <- D %>% colData() %>% as_tibble %>% filter(colData(D)$Sample_ID%in%com_samples) %>%
    pull(Patient_ID)
  res %<>% mutate(Patient_ID = sub("*.+_", "", patIDs))
  res <- res %>% # Specify data frame
    group_by(Patient_ID) %>% # Specify group indicator
    summarise_at(vars(-group_cols()), # Specify columns
                 mean) # Specify function mean per patient
  return(res)
})

# Correlation matrices ----
cor_mat <- full_join(data_mat$Proteins, data_mat$Metabolites, by='Patient_ID')
cor_mat %<>% select(-Patient_ID)

# GGM - Biomics ----
this_mat <- full_join(data_mat$Proteins, data_mat$Metabolites, by='Patient_ID') %>% 
  select(-Patient_ID)
pcor_mat <- ggm.estimate.pcor(as.matrix(this_mat), method = "dynamic", verbose = F)
pval_mat <- network.test.edges(pcor_mat, plot = F, verbose = F)
pval_mat$p.adj.bh <- p.adjust(pval_mat$pval, method="BH")
pval_mat$p.adj.bon <- p.adjust(pval_mat$pval, method="bonferroni")

# # you may want to save the data 
# # make sure result directory exists
# dir.create("results/GGM", showWarnings = T, recursive = T)
# save(pcor_mat, pval_mat, mol_types, data_mat, file=paste0("results/GGM/",datestamp, "_GIM_Plasma_Biomics_GGM.RData"))


# get stats from case control analysis ------------------------------------
get_compiled_stats <- function(sub_set, outcome_names=NULL) {
  
  stat_list<- lapply(sub_set, FUN=function(x){
    # reading metadata
    res <- x %>%
      excel_sheets() %>%
      purrr::set_names() %>%
      map(read_excel, path = x)
    return(res)
  })
  names(stat_list) <- names(sub_set)
  stat_data <- lapply(1:length(stat_list[[1]]), FUN=function(i){
    res <- do.call(rbind, list(stat_list[[1]][[i]], stat_list[[2]][[i]]))
    res %<>% mutate (pval_score = -1*log10(p.adj),
                     asso_type = case_when(statistic > 0 ~ 1, 
                                           statistic < 0 ~ -1),
                     sig = case_when(pval_score >= 1.3 ~ "sig", 
                                     pval_score < 1.3 ~ "nsig"),
                     pval_sign = pval_score * asso_type
    )
    return(res)
  })
  names(stat_data) <- names(stat_list[[1]])
  cytoscape_data <- lapply(1:length(stat_data), FUN=function(x){
    res <- stat_data[[x]][order(stat_data[[x]]$var), ] %>% 
      select(pval_score, asso_type, sig, pval_sign)
    names(res) <- paste0(gsub("-", "_", names(stat_data)[x]), "_", names(res))
    res <- bind_cols(res, stat_data[[x]][order(stat_data[[x]]$var), ] %>% 
                       select(var, feat_col))
    return(res)
  }) %>% plyr::join_all(by=c('var', 'feat_col'), type='left')
  cytoscape_data %<>% dplyr::rename(mol_uname=var, mol_name=feat_col)
  return(cytoscape_data)
}

get_stats_list <- function(sub_set, outcome_names) {
  
  stat_list<- lapply(sub_set, FUN=function(x){
    # reading metadata
    res <- x %>%
      excel_sheets() %>%
      purrr::set_names() %>%
      map(read_excel, path = x)
    
    sheet_ids <- sapply(outcome_names, FUN=function(x) 
      grep(paste0('^', x),names(res)))%>% unlist()
    
    res <- res [sheet_ids]
    return(res)
  })
  return(stat_list)
}

library(tidyverse)
library(magrittr) # %<>%
library(readxl) # reading all sheets at once
library(openxlsx) # convertToDate

# all stats files 
stat_files.caco <- list.files(paste0(root, "/2a_casecontrol/results/"), 
                              pattern="^GIM_..*.xlsx", 
                              full.names = T, 
                              all.files = F)  

# GIM groups/ case-control
sub_set = stat_files.caco[ -grep(pattern = "pathway", stat_files.caco) ]

cytoscape_data <- get_compiled_stats(sub_set) %>%
  mutate(mol_uname = make.names(mol_name, unique=T))


# # you may want to save cytoscape data 
# # save cytoscape data 
# # create stats folder for GGM analysis and save the data there
# dir.create("results/GGM/stats", showWarnings = T, recursive = T)
# save(cytoscape_data, 
#      file = paste0('results/GGM/stats/', datestamp, '_GIM_Case-Control_Compiled_Stat_Data.RData'), 
#      version = 2) 

# figure 3A ---------------------------------------------------------------
library(igraph)

# function to compile network data minimally, given the FDR threshold
load_gim_network_data_minimally <- function(
  p_th_for_ggm = NA
){
  # pvalue threshold
  ggm_thresh <- p_th_for_ggm
  
  # transform the pcor mat into edge list
  tmp <- pcor_mat %>% graph_from_adjacency_matrix(mode='undirected', weighted = T) %>% igraph::simplify()
  ggm_edges <- cbind.data.frame(get.edgelist(tmp), edge_attr(tmp)$weight)
  names(ggm_edges) <- c("source", "target", "pcor_val")
  
  # select edges based on threshold
  ggm_edges %<>% filter(abs(pcor_val)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bh<=ggm_thresh]))) #%>% 
  # remove within omic edges
  #filter(!(source%in%mol_types$Proteins$unique_name & target%in%mol_types$Proteins$unique_name)) %>% 
  #filter(!(source%in%mol_types$Metabolites$unique_name & target%in%mol_types$Metabolites$unique_name))
  
  # add edge_type based on sign of pcor
  ggm_edges %<>% mutate(edge_type = case_when(pcor_val > 0 ~ "pos", 
                                              pcor_val < 0 ~ "neg"))%>%
    data.frame(.,stringsAsFactors=FALSE)
  # build a node attribute data frame
  node_attributes <- data.frame(node_name=unique(c(as.matrix(ggm_edges$source), as.matrix(ggm_edges$target))), 
                                node_type="Metabolite") %>% 
    mutate(node_type = case_when(node_name%in%mol_types$Proteins$unique_name ~ "Protein", 
                                 TRUE ~ "Metabolite")) %>% 
    left_join(cytoscape_data, by=c("node_name"= "mol_uname"))
  
  node_attributes %<>% mutate(id=as.character(1:nrow(node_attributes)))
  
  ggm_edges <- node_attributes %>% dplyr::select(id, node_name) %>% 
    left_join(ggm_edges,.,by=c("source"= "node_name")) %>% dplyr::rename(from=source, source=id)
  ggm_edges <- node_attributes %>% dplyr::select(id, node_name) %>% 
    left_join(ggm_edges,.,by=c("target"= "node_name")) %>% dplyr::rename(to=target, target=id)
  
  ggm_edges <- data.frame(ggm_edges, stringsAsFactors=FALSE)
  node_attributes <- data.frame(node_attributes, row.names = node_attributes$id, stringsAsFactors = F)
  
  # add node types 
  ggm_edges$from_type = c("protein", "metabolite")[2-ggm_edges$from %in% mol_types$Proteins$unique_name]
  ggm_edges$to_type = c("protein", "metabolite")[2-ggm_edges$to %in% mol_types$Proteins$unique_name]
  
  as.list(environment())
}

# with fdr threshold 0.2
nod = load_gim_network_data_minimally(0.2)

# ggm data and stat associations 
df_ggm = nod$ggm_edges[,c("from","to","pcor_val","edge_type","from_type","to_type")]
df_nod = nod$node_attributes[,c("node_name","node_type","casecontrol_pval_score",
                                "casecontrol_asso_type","casecontrol_sig","mol_name")]

g0 = igraph::graph_from_data_frame( df_ggm, directed = F )
g = mst( g0 )
g=g0
gl = mst( g0 )

# layout based on MST of whole network
set.seed(42)
LO = layout_with_fr( mst( g0 ) )

vshape = (structure( df_nod$node_type, names = df_nod$node_name) == "Protein")+1
V(g)$shape = factor(vshape[ names( V(g) ) ], labels = c("metabolite","protein")) 

vcolor= c("white","red")[((structure( df_nod$casecontrol_sig, names = df_nod$node_name) == "sig")+1)[names( V(g))]]
V(g)$color = vcolor

# store node attributes
xnod = nod$node_attributes
rownames(xnod) = xnod$node_name
for(i in colnames(xnod)){
  print(i)
  V(g) = igraph:::`$<-.igraph.vs`(V(g), i,  xnod[names(V(g)),i] )
}
rm(i,xnod)

# node clusters based on louven 
grps = cluster_louvain(g)$membership %>% {lapply(seq(max(.)), function(i) which(.==i))}

# # to save
# svg(filename = "network3a.svg",width=15, height=10)

# no names
plot(
  g, 
  edge.color = adjustcolor( "black", alpha.f = 0.23),
  #mark.groups = grps,  
  mark.col = adjustcolor( "black", alpha.f = 0.05),
  mark.border = NA,
  vertex.border = NA,
  vertex.shape = (shapes()[c(1,9)])[as.numeric(V(g)$shape)],
  layout = LO,
  vertex.color = V(g)$color,
  vertex.label =NA, 
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  vertex.size = 2
)

# dev.off()

# g is the figure 3A
net1 = g 
LO_net1 = LO

# figure 3B ---------------------------------------------------------------

nod = load_gim_network_data_minimally( 0.2)

df_ggm = nod$ggm_edges[, c("from",
                           "to",
                           "pcor_val",
                           "edge_type",
                           "from_type",
                           "to_type") ]

df_nod = nod$node_attributes[, c("node_name",
                                 "node_type",
                                 "casecontrol_pval_score",
                                 "casecontrol_asso_type",
                                 "casecontrol_sig",
                                 "mol_name") ]

# base graph
g0 = igraph::graph_from_data_frame( df_ggm, directed = F )
g = mst( g0 )

vshape = (structure( df_nod$node_type, names = df_nod$node_name) == "Protein")+1
V(g)$shape = factor( vshape[ names( V(g) ) ], labels = c("metabolite","protein")) 

vcolor= c("white","red")[((structure( df_nod$casecontrol_sig, names = df_nod$node_name) == "sig")+1)[names( V(g))]]
V(g)$color = vcolor

grps = cluster_louvain(g)$membership %>% {lapply(seq(max(.)), function(i) which(.==i))}

# MST of hits 
ii = df_nod$node_name[ df_nod$casecontrol_sig == "sig" ]
dit = igraph::distances(g0)[ii,ii]
g3 = mst(igraph::graph_from_adjacency_matrix(dit,weighted = TRUE ,mode = "undirected"))

grps = cluster_louvain(g3)$membership %>% {lapply(seq(max(.)), function(i) which(.==i))}
vshape = (structure( df_nod$node_type, names = df_nod$node_name) == "Protein")+1
V(g3)$shape = vshape[ names( V(g3) ) ] 
vcolor= c("white","red")[((structure( df_nod$casecontrol_sig, names = df_nod$node_name) == "sig")+1)[names( V(g3))]]
V(g3)$color = vcolor

# load proteomics data
D<-
  (function(){
    load(paste0(root,"/1_preprocessing/data_preprocessed_proteo.Rdata"))
    return(D)
  })()

# we need panel information 
aa = rowData(D)
aa = structure( gsub("Olink | III| II| I", "", aa$Panel), names = make.names(aa$name))

df_nod$panel = aa[df_nod$node_name]
df_nod$panel[is.na(df_nod$panel)] = df_nod$node_type[is.na(df_nod$panel)]

# 1 cardio, 2 inflam, 3 Metabolites 
vcolor= c("tomato", "orange", "cadetblue3" )[ structure( df_nod$panel, names = df_nod$node_name)[names( V(g3))] %>% factor %>% as.numeric ]
V(g3)$color = vcolor

# store node attributes
xnod = nod$node_attributes
rownames(xnod) = xnod$node_name
for(i in colnames(xnod)){
  print(i)
  V(g3) = igraph:::`$<-.igraph.vs`(V(g3), i,  xnod[names(V(g3)),i] )
}
rm(i,xnod)
# panel info
V(g3)$panel = unname( c(tomato = "Cardiovascular", orange = "Inflammation", cadetblue3 = "Metabolites" )[V(g3)$color] )

# shortest path for mst 
E(g3)$d_shortest_path = apply(get.edgelist(g3), 1, function(i) dit[i[1], i[2]])


k = 55 #755
set.seed(k)
# E(g3)$weight = 1
LO = layout_with_fr(g3)

# # to save
# svg(filename = "network_nonames.svg",width=9, height=9)

# no names
plot(g3, 
     edge.color = adjustcolor( "black", alpha.f = 0.7),
     vertex.border = NA,
     vertex.shape = (shapes()[c(1,9)])[V(g3)$shape],
     vertex.color = V(g3)$color,
     mark.groups = grps,  
     mark.col = adjustcolor( "black", alpha.f = 0.05),
     mark.border = NA,
     layout = LO,
     vertex.label = NA,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     vertex.size=2 )
# dev.off()

net2= g3
LO_net2 = LO
# rm(list = setdiff(ls(), c("net1", "net2")))

# save networks
dir.create("results")
save(file = "results/igraph_networks.Rdata", net1, net2, LO_net1, LO_net2)
 

# # to export to cytoscape --------------------------------------------------
# library(igraph)
# library(tidyverse)
# library(RCy3) # make sure cytoscape is running
# library(magrittr)
# 
# 
# # function to extract edge and vertex data from a network 
# f_net <- function(net){
#   dfe = as.data.frame(get.edge.attribute(net))
#   dfv = as.data.frame(get.vertex.attribute(net))
#   list(ed = dfe, ve = dfv)
# }
# 
# nett = net2
# nod <- f_net(nett)
# 
# this_collection <- sprintf('%s_%.2f_OR_Lower', "BH", 0.2)
# this_net <- "casecontrol_simplified"
# outcome_name = "casecontrol"
# 
# # load the network to cytoscape
# createNetworkFromDataFrames(edges=get.edgelist(nett) %>% {colnames(.)= c("source", "target");. } %>% data.frame, 
#                             title=this_net, 
#                             collection=this_collection)
# 
# # workaround for a bug 
# ggm_edges = cbind(nod$ed, get.edgelist(nett) %>% {colnames(.)= c("source", "target");. } %>% data.frame)
# ggm_edges <- ggm_edges %>% mutate(key=paste(source, "(interacts with)", target))
# 
# # get SUID for edges and match with edge attributes
# cy_edges <- getTableColumns(table = 'edge')
# cy_edges <- cy_edges[order(cy_edges$name), ]
# ggm_edges <- ggm_edges[order(ggm_edges$key), ]
# ggm_edges$cpSUID <- cy_edges$SUID
# 
# # add edge attributes
# loadTableData(subset(ggm_edges, select = -c(source, target)), 
#               data.key.column = 'cpSUID', table.key.column = 'SUID',
#               table = 'edge')
# 
# # add node attributes
# node_attributes = nod$ve
# loadTableData(node_attributes, data.key.column = 'node_name', table = 'node')
# 
# #then prepare style variables
# style.name = this_net
# nodeLabels <- mapVisualProperty('node label','mol_name','p')
# 
# nodeShapes <- mapVisualProperty('node shape','node_type','d',
#                                 c('Protein',  'Metabolite'),
#                                 c('DIAMOND', 'ELLIPSE'))
# edgeStyles <- mapVisualProperty('edge line type', 'edge_type', 'd', 
#                                 c("neg", "pos"), c("LONG_DASH", "SOLID"))
# 
# 
# 
# #and then create the style
# createVisualStyle(style.name=style.name, base.url = 'http://localhost:1234/v1',
#                   mappings = list(nodeLabels,nodeShapes,edgeStyles)) 
# 
# #finish by applying the style
# setVisualStyle(style.name=style.name, base.url = 'http://localhost:1234/v1',
#                network = this_net)
# max_val <- max(abs(node_attributes[(sprintf('%s_pval_sign',outcome_name))]), na.rm=T)
# setNodeColorMapping(
#   table.column = as.name(sprintf('%s_pval_sign', outcome_name)),
#   table.column.values = c(-1*max_val, 0, max_val),
#   mapping.type = "c",
#   colors= c('#E74C3C' ,'#F0F3F4', '#3498DB'),
#   style.name=style.name, base.url = 'http://localhost:1234/v1',
#   network = this_net)












