
# set wd to here
setwd(glue::glue("{root}/2d_supplements"))
library(SummarizedExperiment)
library(magrittr)

# load the results of case-control analysis 
res_paths = grep(".rds", list.files(glue::glue("{root}/2a_casecontrol/results/")), value = T)
load(paste0(glue::glue("{root}/2a_casecontrol/results/"), res_paths[1]))
Ds = list(metabo = D)
load(paste0(glue::glue("{root}/2a_casecontrol/results/"), res_paths[2]))
Ds$proteo = D
rm(D)


# PCAs --------------------------------------------------------------------
ggs_pca<-
  lapply(Ds, function(D){
    # data 
    x = D %>% assay %>% t %>% as.matrix() %>% scale
    
    df= colData(D)[,c("Group", "Status", "batch", "Time")] %>% as.data.frame
    dfx = data.frame(df, prcomp(x)$x)
    
    list(
      status = # PCA for status
        ggplot(dfx, aes(x = PC1, y = PC2, color  = Status)) + 
        geom_point() + 
        theme_minimal(20),
      batch = #  batch effect 
        ggplot(dfx, aes(x = PC1, y = PC2, color = batch)) + 
        geom_point()+
        theme_minimal(20)
    )
  })



# Time analysis -----------------------------------------------------------
res = lapply(Ds, function(D){
  X = D %>% assay %>% t %>% as.matrix() %>% scale
  colnames(X) = make.names( rowData(D)$name )
  
  df= colData(D)[,c("Patient_ID","Group", "Status", "batch", "Time")] %>% as.data.frame
  is_proteo = if(length(table(df$batch))==1) TRUE else FALSE
  
  df$Patient_ID %<>% factor
  ps<-
    apply(X, 2, function(x){
      #sm = lm(x~Time+Status, data.frame(x, df)) %>% summary
      m1 = if(!is_proteo) x~Time+Status+(1|Patient_ID)+(1|batch) else x~Time+Status+(1|Patient_ID)
      m0 = if(!is_proteo) x~Status+(1|Patient_ID)+(1|batch) else x~Status+(1|Patient_ID)
      
      fit1 = lmerTest::lmer(m1, data.frame(x=x, df))
      sm = fit1 %>% summary
      sm = sm$coefficients[c("Timed2","Timed3"),c("Estimate","Pr(>|t|)")]
      ifelse(sign(sm[1,"Estimate"]) == sign(sm[2,"Estimate"]),{
        #min(sm[,"Pr(>|t|)"])
        anova(lmerTest::lmer(m0, data.frame(x=x, df)), fit1)$`Pr(>Chisq)`[2]
      }, 1)
    })
  
  # slack to Jan most extreme one 
  table( p.adjust(ps) < 0.05 ) %>% print
  
  rm(D)
  pplot <- function(){
    p.adjust(ps) %>% sort %>% 
      plot(log ="y", ylab = "adjusted p-value", pch = "*", 
           xlab = if(!is_proteo) "metabolites" else "proteins")
    abline(v = sum(p.adjust(ps) < 0.05) ,h = 0.05, col = "red")
    mtext(side=3, line=3, at=-0.07, adj=0, cex=1, 
          if(!is_proteo) "effect of time to metabolite measurements"
          else "effect of time to protein measurements" )
    mtext(side=3, line=2, at=-0.07, adj=0, cex=0.7, 
          "adjusted p-values for Time in x~Time+Status+(1|Patient_ID)")
  }
  env = environment(pplot)
  
  gg4 = 
    order(ps)[1:4] %>% lapply(function(i)
      data.frame(x = X[,i], df[,c("Time","Status","Patient_ID")], met = colnames(X)[i])
    ) %>% do.call(what = rbind) %>% 
    ggplot(aes(y = x, x = paste(Time, Status), color = Time)) + 
    geom_jitter(width = 0.075, height = 0, pch = 21) +
    geom_line(aes(group = Patient_ID),  alpha = 0.35, color = "gray") +
    geom_boxplot(aes(color = Time), fill = NA, width = 0.38, outlier.color = NA) + 
    facet_wrap(met~.,ncol = 2,scales = "free") +
    theme_minimal(16) +
    theme(strip.background = element_rect(fill = "wheat"), 
          panel.background = element_rect())
  
  
  # all metabolites are plotted
  dfm <-
    order(ps) %>% lapply(function(i)
      data.frame(x = X[,i], df[,c("Time","Status","Patient_ID")], met = colnames(X)[i])
    ) %>% do.call(what = rbind) 
  
  dfm$met = factor(dfm$met, levels = colnames(X)[order(ps)] )
  
  gg_all = 
    ggplot(dfm[dfm$Status=="COVID",], aes(y = x, x = paste(Time, Status), color = Time)) + 
    geom_jitter(width = 0.1, height = 0, pch = 21, size = 0.3, alpha = 0.35) +
    geom_line(aes(group = Patient_ID),  alpha = 0.15, color = "gray") +
    geom_boxplot(aes(color = Time), fill = NA, width = 0.38, outlier.color = NA,notch = T) + 
    facet_wrap(met~.,ncol = 14,scales = "free") +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "wheat"), 
          panel.background = element_rect(),
          strip.text = element_text(size = 5),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5)) +
    labs(x = "Time", 
         y = if(!is_proteo) "scaled metabolite abundance"
         else "scaled protein abundance")
  list(ps=ps, pplot = pplot, gg_most4=gg4, gg_all = gg_all )
})

rm(Ds)

# # 1000 x 755
# # 800 X 600
# 
# # metabolomics
# res$metabo[-4]
# res$metabo$pplot()
# res$metabo$gg_most4
# # proteomics
# res$proteo[-4]
# res$proteo$pplot()
# res$proteo$gg_most4

# save only scatter plot and example boxplots 
ggs_time = res
ggs_time$metabo[4] <- ggs_time$proteo[4] <- NULL
dir.create("results")
save(file = "results/supplement_PCA_Time.fig", ggs_time, ggs_pca)
