#' 
#' generates and plots paper figures
#' 

rm(list = setdiff(ls(),"root")) 
library(ggplot2)
library(ggrepel)
library(magrittr)
library(igraph)

# root directory 
# root = "~/Documents/R/_github_/covid-omics/"

# set folder to paper figures
wd = glue::glue("{root}/5_paperfigures")
setwd(wd)
      
# figure 2
if(!file.exists("figure2.fig")){
  source(glue::glue("{root}/5_paperfigures/figure2.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}

# figure 3
if(!file.exists("figure3.fig")){
  source(glue::glue("{root}/5_paperfigures/figure3.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}

# figure 4
if(!file.exists("figure4.fig")){
  source(glue::glue("{root}/5_paperfigures/figure4.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}

# figure 5
if(!file.exists("figure5.fig")){
  source(glue::glue("{root}/5_paperfigures/figure5.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}
# performance of the  model on other datasets
if(!file.exists("figure5replication.fig")){
  source(glue::glue("{root}/5_paperfigures/figure5replication.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}
# supplement i.e. lasso path: parsimony vs accuracy tradeoff 
if(!file.exists("figure5Supplement.fig")){
  source(glue::glue("{root}/5_paperfigures/figure5Supplement.R"))
  rm(list = setdiff(ls(), c("root","wd")))
  setwd(wd)
}



#### output generated figures as PDFs

path <- glue::glue("{root}/5_paperfigures/")

# Figure 2
load(glue::glue("{path}/figure2.fig"))
pdf(file = glue::glue("{path}/figure2.pdf")); 
sapply(ggs_volcano, plot)
sapply(ggs_boxplots, plot); 
sapply(ggs_pwbarplots, plot); 
dev.off()

# Figure 3
load(glue::glue("{path}/figure3.fig"))
library(igraph)
library(magrittr)
pdf(file = glue::glue("{path}/figure3.pdf"), width = 14, height = 7)
plot_networks() # this is an internally defined function that plots the networks
dev.off()

# Figure 4
load(glue::glue("{path}/figure4.fig"))
pdf(file = glue::glue("{path}/figure4.pdf"))

plot(gg_lollipops)

# pheatmap objects 
plot.new()
hitheatmaps$metabo
plot.new()
hitheatmaps$proteo

# volcano plots
sapply(ggs_volcano_clins$metabo, plot)
sapply(ggs_volcano_clins$proteo, plot)
dev.off()

# Figure 5 + replication
load(glue::glue("{path}/figure5.fig"))
load(glue::glue("{path}/figure5replication.fig"))
pdf(file = glue::glue("{path}/figure5.pdf"))
plot(gg_figure5)
plot(gg_su)
plot(gg_shen)
ggplot() + 
  annotate("text",x = 1,y = 1, size = 8, label = "Note: Figure 5C was manually created.") + 
  theme_void()
dev.off()

# Figure 5 supplement
load(glue::glue("{path}/figure5Supplement.fig"))
pdf(file = glue::glue("{path}/figure5_supplement.pdf"))
plot(gg_lpath)
dev.off()

setwd(root)

