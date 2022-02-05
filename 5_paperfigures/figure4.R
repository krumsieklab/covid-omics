
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(scales)
library(pheatmap)

# load the results of case-control analysis 
load(glue::glue("{root}/2b_clinical/results/GIM_clinicals_metabo.rds"))
Ds = list(metabo = D)
load(glue::glue("{root}/2b_clinical/results/GIM_clinicals_proteo.rds"))
Ds$proteo = D
rm(D)

# annotations of outcomes
df_outc = 
  rbind( 
    data.frame(
      outc = c("lymphocyte", "platelet", "ferritin_level", "crp", "d_dimer"), 
      type = "lab", stringsAsFactors = F),
    data.frame(
      outc = c("ckd_or_esrd", "dm", "htn", "max_SOFA"), 
      type = "predisposition", stringsAsFactors = F),
    data.frame(
      outc = c("age", "sex", "bmi"), 
      type = "demographics", stringsAsFactors = F),
    data.frame(
      outc = c("death", "intubation"), 
      type = "events", stringsAsFactors = F)
    
  )

# results of association analysis 
res = lapply(Ds, function(D){
  re = maplet::mtm_res_get_entries(D, "stats")
  nm = sapply(re, function(x) x$output$table$term[1])
  names(re) = nm
  varnames = make.names(rowData(D)$name)
  lapply(re, function(x){
    x =  x$output$table
    x$name = varnames
    x
  })
})

# volcano plots -----------------------------------------------------------
# reverse log transform
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

ggs_volcano_clins <-
  lapply( res, lapply, function(df){
    
    df$p.adj[df$p.adj>0.05] = NA
    df$name[is.na(df$p.adj)] = NA
    if(is.null(df$fc)){
      df$fc = df$statistic
      xlab = "statistics"
    }else{
      xlab = "fold change"
    }
    
    gg = ggplot(df, aes(x = fc, y = p.value, 
                        color =  c("down", "up")[(p.adj < 0.05) * (fc>0)+1]  )) + 
      geom_point() + #aes(alpha = -log10(p.value))  
      scale_y_continuous(
        trans = reverselog_trans(10), 
        breaks = scales::trans_breaks("log10", function(x) 10^x), 
        labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      labs(color = "", x = xlab) +
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


# hit heatmaps  -----------------------------------------------------------
# wrapper function to collect results and return in an usable format for pheatmap
wrap_get_data <- function(D){
  
  # extract all the results data 
  get_all <- function(D){
    re = maplet::mtm_res_get_entries(D, "stats")
    nm = sapply(re, function(x) x$output$table$term[1])
    names(re) = nm
    varnames = make.names(rowData(D)$name)
    lapply(re, function(x){
      x =  x$output$table
      x$name = varnames
      x
    })
  }
  
  res = get_all(D)
  # check if order is correct 
  sapply(seq(length(res)-1), function(i) identical(res[[i]]$var, res[[i+1]]$var))
  
  mm = sapply(res, function(x) x$p.adj) %>% {rownames(.) = res[[1]]$var;. }
  colnames(mm)[2] = "intubation"
  
  # annotations of outcomes
  df_outc = 
    rbind( 
      data.frame(
        outc = c("lymphocyte", "platelet", "ferritin_level", "crp", "d_dimer"), 
        type = "lab", stringsAsFactors = F),
      data.frame(
        outc = c("ckd_or_esrd", "dm", "htn", "max_SOFA"), 
        type = "predisposition", stringsAsFactors = F),
      data.frame(
        outc = c("age", "sex", "bmi"), 
        type = "demographics", stringsAsFactors = F),
      data.frame(
        outc = c("death", "intubation"), 
        type = "events", stringsAsFactors = F)
      
    )
  
  df_outc$outc
  m1 = mm[,df_outc$outc ] 
  
  # m1: what we want
  # mm: all results including ones that need to be removed
  # df_outc: annotation of outcomes
  list(m1 = m1, mm = mm, df_outc = df_outc)
}

nods = lapply(Ds, wrap_get_data)
# distance based on significance yes no, given pvalues, x
fd <- function(x, p_th = 0.05) dist(1*(x< p_th), method = "binary")


# metabolomics ------------------------------------------------------------
m1 = nods$met$m1
df_outc = nods$met$df_outc
df_outc$type = factor(df_outc$type, levels = unique(df_outc$type))

# 7 metabolites has no hits, delete them
m1 = m1[rowSums(m1 < 0.05)>0,] 

X = t(m1)
df_row = df_outc[,-1, drop = F] %>% {rownames(.) = df_outc$outc;.}
identical(rownames(df_row), rownames(X))
identical(df_outc$outc, rownames(X))
new_names = structure( 
  c("Lymphocyte", "Platelet", "Ferritin", "CRP**", "D-dimer", "Kidney disease", 
    "DM*", "Hypertension", "SOFA", "Age", "Sex", "BMI", "Death", "ARDS***"),
  names = rownames(df_row))
rownames(df_row) <- rownames(X) <- new_names[rownames(df_row)]

df_outc$outc = rownames(df_row)
ann_color <- list(
  type = structure( hue_pal()(8)[-seq(4)], names =  levels(df_outc$type))
)

# back up for barplot
Xmet = X

# some configurations for metabolomics
(p1 = pheatmap(log10(X) %>% abs, 
               breaks = c(0, -round(log10( 0.05),1),4,8,16,32), 
               clustering_distance_cols = fd(m1),
               clustering_distance_rows = fd(t(m1)),
               clustering_method = "ward.D2",
               cutree_cols = 4,
               cutree_rows = 2,
               legend_breaks =  c(0, -round(log10( 0.05),1),4,8,16,32),
               annotation_row = df_outc[,-1, drop = F] %>% {rownames(.) = df_outc$outc;.},
               annotation_colors =  ann_color,
               legend_labels =  c("", "0.05", "1e-4", "1e-8", "1e-16", "1e-32"), 
               scale = "none", fontsize_col = 6,show_colnames = F,
               color = c("floralwhite",rev(heat.colors(10))[-c(1:3,9:10)])))

# proteomics --------------------------------------------------------------
m1 = nods$pro$m1
df_outc = nods$pro$df_outc
df_outc$type = factor(df_outc$type, levels = unique(df_outc$type))

# prot info
dp = Ds$pro %>% rowData()
identical(dp$feature_id, rownames(m1))
rownames(m1) = dp$name
df_prot = data.frame("Olink Panel" = dp$Panel %>% gsub(pattern = "Olink| II| III", replacement = "" ) %>% 
                       tolower, check.names = F)
rownames(df_prot) = dp$name


# 126 proteins has no hits, delete them
m1 = m1[rowSums(m1 < 0.05)>0,] 

X = t(m1)
df_row = df_outc[,-1, drop = F] %>% {rownames(.) = df_outc$outc;.}
identical(rownames(X), rownames(df_row))
rownames(X) <- rownames(df_row) <- new_names[rownames(X)]

df_outc$outc = rownames(df_row)
ann_color <- list(
  type = structure( hue_pal()(8)[-seq(4)], names =  levels(df_outc$type))
)
ann_color$`Olink Panel` = structure( c("darkslategrey", "gray"), names =  levels(factor(df_prot$`Olink Panel`)))

# back up for barplot
Xpro = X

# manual annotation of protein panels
ids <-
  c(
    "CCL3",
    "CXCL1",
    "FGF-21",
    "FGF-23",
    "IL18",
    "IL6",
    "SCF",
    "MCP-1",
    "OPG",
    "uPA"
  )

df_prot[ids, ] <- c(" inflammation", " cardiovascular", " cardiovascular",
                    " cardiovascular", " cardiovascular", " inflammation",
                    " cardiovascular", " inflammation", " cardiovascular", " cardiovascular" )
#---
(p2 = # some configurations for proteomics
    pheatmap(log10(X) %>% abs, 
             breaks = c(0, -round(log10( 0.05),1),2,3,4,6), 
             legend_breaks =  c(0, -round(log10( 0.05),1),2,3,4,6),
             annotation_row = df_outc[,-1, drop = F] %>% {rownames(.) = df_outc$outc;.}, 
             annotation_col = df_prot, 
             legend_labels =  c("", "0.05", "0.01", "0.001", "1e-4", "1e-6"), 
             clustering_distance_cols = fd(m1),
             clustering_distance_rows = fd(t(m1)),
             clustering_method = "ward.D2",
             cutree_cols = 4,
             cutree_rows = 2,
             annotation_colors =  ann_color,
             scale = "none", show_colnames = F, fontsize_col = 6,
             color = c("floralwhite",rev(heat.colors(10))[-c(1:3,9:10)])))


# lollipop plots ----------------------------------------------------------
# get the data
mm = sapply(nods, function(x) colSums( x$m1 < 0.05 ))
identical( rownames(Xmet) %>% names, rownames(mm) )
rownames(mm) <- rownames(Xmet)

plot(mm %>% apply(2, function(x) x/max(x)) )
abline(0,1,lty=2, col ="red")

ddf = reshape2::melt(mm)
ddf$Var1 =  factor(ddf$Var1, levels = rev(c(levels(ddf$Var1)[7], levels(ddf$Var1)[-7]) ))
ddf$Var2 = factor(ddf$Var2, labels = c("Metabolites","Proteins"))

gg_lollipops <-
  ddf %>% 
  ggplot(aes(y = Var1, x = value, color  = Var2)) + 
  geom_segment(aes(xend = 0, yend = Var1, group = Var1)) + #, lwd = 2, alpha = 0.5
  geom_point(size =  2)+
  # geom_point(size =  0.3, color = "red")+
  geom_line()+
  facet_wrap(~Var2, nrow = 1, scales = "free_x") + 
  theme_minimal() +
  theme(legend.position = "n", 
        strip.background = element_rect(fill = "wheat"),
        panel.background = element_rect(colour = "black"),
        axis.ticks.x = element_line(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  scale_color_manual(values = c("blue", "firebrick")) + 
  labs(y = "", x = "number of significant features")

#
# gg2 = gg_lollipops 
# # save(file = "lollipops.Rdata",gg_lollipops)
# 
# # svg(filename = "hits_comparison_counts.svg",width=8, height=4)
# gg2
# # dev.off()
# 
# # svg(filename = "hits_comparison_counts_met.svg",width=4, height=4)
# gg_met = gg2 %>% {.$data = gg2$data[gg2$data$Var2 == "Metabolites",];.} + scale_color_manual(values = "blue")
# # dev.off()
# 
# # svg(filename = "hits_comparison_counts_pro.svg",width=4, height=4)
# gg_pro = gg2 %>% {.$data = gg2$data[gg2$data$Var2 != "Metabolites",];.} + scale_color_manual(values = "firebrick")
# # dev.off()




# death, kidney disease,  c-reactive protein
ggs_volcano_clins$metabo[c("death", "ckd_or_esrd","crp")]

# death, platelet count, ferritin
ggs_volcano_clins$proteo[c("death", "platelet", "ferritin_level")]


# lollipop plots 
gg_lollipops
# pheatmap 1
p1
# pheatmap 2
p2


# # save the plots --------------------------------------------------------
hitheatmaps = list(metabo = p1, proteo = p2)

# set the names
ggs_volcano_clins$metabo = ggs_volcano_clins$metabo[c("death", "ckd_or_esrd","crp")]
names(ggs_volcano_clins$metabo) = c("death", "kidney disease",  "c-reactive protein")

# set the names
ggs_volcano_clins$proteo = ggs_volcano_clins$proteo[c("death", "platelet", "ferritin_level")]
names(ggs_volcano_clins$proteo) = c("death", "platelet",  "ferritin")

# add title 
for(i in names(ggs_volcano_clins$metabo)){
  ggs_volcano_clins$metabo[[i]] =  ggs_volcano_clins$metabo[[i]] + ggtitle(i)
}
for(i in names(ggs_volcano_clins$proteo)){
  ggs_volcano_clins$proteo[[i]] =  ggs_volcano_clins$proteo[[i]] + ggtitle(i)
}

# save the plots
save(file = "figure4.fig", gg_lollipops, hitheatmaps, ggs_volcano_clins)















