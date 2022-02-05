library(magrittr)
library(survival)
library(coxme)

wd = getwd()
# set working directory to location of source code
setwd(glue::glue("{root}/3b_prediction"))
# load pre-calculated model fits and data
load("lasso_loo_fits.Rdata")

# # regularized linear mixed effect model
# functions needed to cast problem to be solved with cox(ph)/(me)
source("./model/utils_coxph_sorceries.R")
# ridge regularization 
source("./model/utils_coxme_sorceries.R")
# lasso regularization
source("./model/utils_lasso_sorceries.R")
# utility function to calculate C-index for loo
source("./model/utils_loo_prediction.R")

# concatanete constituents of composite outcome
colnames(compo)[ colnames(compo)%in%c("intubation1", "max_aki_stage" ) ] = c("intubation", "max_aki")
df = dplyr::left_join(df, compo[,-ncol(compo)], by = c("ID"="Patient_ID")) 

# Hypoxia: order on O2_device, then order based on Fi_O2 within each device
Hypoxia = rank( nrow(df)^as.numeric(factor(df$O2_device)) + 
                  as.numeric(factor(df$Fi_O2)), na.last = "keep" )
# # check it out 
# par(mfrow = c(1,2))
# plot(Hypoxia~df$O2_device)
# plot(Hypoxia~df$Fi_O2)
df$Hypoxia = Hypoxia

# resulting loo fits
res_nodes =  list(base = nods_base_loo, `base + metabolites`= nods_loo )

# name namer for named list to be returned 
fnm <- function(x) structure(x, names = x)

# predicted scores 
# all prediction based on r0(composite score) model
# results collected as list
resl <-
  lapply(names(res_nodes) %>% fnm, function(setting){ 
    dset = "test"
    re<-
      sapply( c("r0", "death", "intubation", "max_aki", "Hypoxia", "hospitalization") %>% fnm, function(outc){
        j = "r0"
        cat("\n", setting, " - ", j," - " ,outc, "\n")
        
        fit = calculate_d_loo(res_nodes[[setting]][[j]], df, outc, dset,f = max, return_data = T)
        list(d= fit$fit$concordance, se = sqrt(fit$fit$var), nod = fit)
        
      }) %>% t %>% data.frame
    # names(re$nod) <- "r0"
    re
  })

# plot the results 
library(dplyr)
library(ggplot2)
# rownames
rn = resl$`base + metabolites`$nod$r0$dfx %>% rownames
# retrieve composite outcome for per patient
dfy = 
  data.frame(id = resl$`base + metabolites`$nod$r0$dfx$id, 
             sapply(resl$`base + metabolites`$nod, function(x) x$dfx[rn, "y"])) %>% 
  group_by(id) %>% summarise_all( list(~max(.,na.rm = T))) %>% 
  as.matrix() %>% {.[is.infinite(.)]=NA;.} %>% data.frame() 


# (base)-model prediction
dfyh0 = 
  data.frame(id = resl$base$nod$r0$dfx$id, 
             sapply(resl$base$nod, function(x) x$dfx[rn, "yh"])) %>%
  group_by(id) %>% summarise_all( list(~max(.,na.rm = T))) %>% 
  as.matrix() %>% {.[is.infinite(.)]=NA;.} %>% 
  data.frame()


# (base+metabolites)-model prediction
dfyh1 = 
  data.frame(id = resl$`base + metabolites`$nod$r0$dfx$id, 
             sapply(resl$`base + metabolites`$nod, function(x) x$dfx[rn, "yh"])) %>%
  group_by(id) %>% summarise_all( list(~max(.,na.rm = T))) %>% 
  as.matrix() %>% {.[is.infinite(.)]=NA;.} %>% 
  data.frame()

# compile results as summary table
sm <- sapply( 
  structure( colnames(dfy)[c(-1,-8)], names = colnames(dfy)[c(-1,-8)]), 
  function(i){
    print(i)
    d0 = lm(y~yh+0, data.frame(y = dfy[[i]], yh = dfyh0[[i]] ))
    d1 = lm(y~yh+0, data.frame(y = dfy[[i]], yh = dfyh1[[i]] ))
    
    ctest <- concordance(d0, d1)
    contr <- c(-1, 1)
    dtest <- contr %*% coef(ctest)
    dvar <- contr %*% vcov(ctest) %*% contr
    d0 = concordance(d0)
    d1 = concordance(d1)
    c(contrast=dtest, sd=sqrt(dvar), z=dtest/sqrt(dvar), 
      d0 = d0$concordance, se0 = sqrt(d0$var), 
      d1 = d1$concordance, se1 = sqrt(d1$var) )
  }
) %>% t

sm = data.frame( v = rownames(sm), sm, 
                 p = format.pval( pnorm(abs(sm[,"z"]), lower.tail = F), digits = 1) )
sm$v = factor(sm$v, 
              levels = c("r0", "death", "intubation", "max_aki", "Hypoxia", "hospitalization"),
              labels = c("Composite", "Death", "Intubation", "Aki", "Hypoxia", "Hospitalization"))
sm


# plot the results --------------------------------------------------------

# bootsrapping confidence intervals
boo <- lapply( 
  structure( colnames(dfy)[c(-1,-8)], names = colnames(dfy)[c(-1,-8)]), 
  function(i){
    set.seed(42)
    replicate(200,{
      j = sample(length(dfy[[i]]), replace = T)
      d0 = concordance(y~yh, data.frame(y = dfy[[i]], yh = dfyh0[[i]])[j,])$concordance
      d1 = concordance(y~yh, data.frame(y = dfy[[i]], yh = dfyh1[[i]])[j,])$concordance
      c(d0=d0, d1=d1)
    }) %>% t
  }
) 
ps = seq(0.025,0.475,by = 0.05)
boo = 
  lapply(c(upper = 1, lower = 0), function(k){
    lapply(names(boo), function(i){
      data.frame( apply(boo[[i]],2,quantile, probs = abs(k-ps)), p = seq(ps), v = i)
    }) %>% 
      do.call(what = rbind) %>% 
      reshape2::melt(id = c("p","v") )
  }) %>% {
    mm = .$upper;
    colnames(mm)[4] = "dmax";
    mm$dmin = .$lower$value ;
    mm
  }

colnames(boo)[3] = "M"
boo$M = factor(boo$M, labels = c("base", "base + metabolites"))
boo$v = factor(boo$v, 
               levels = c("r0", "death", "intubation", "max_aki", "Hypoxia", "hospitalization"),
               labels = c("Composite", "Death", "Intubation", "Aki", "Hypoxia", "Hospitalization"))

library(ggplot2)
gg0<-
  rbind(
    sm[,c("v","contrast","sd","z","d0","se0","p")] %>% 
      {colnames(.)[5:6] = c("d","se"); cbind(., M = "base")},
    sm[,c("v","contrast","sd","z","d1","se1","p")] %>% 
      {colnames(.)[5:6] = c("d","se"); cbind(., M = "base + metabolites")}
  )  %>% 
  ggplot( aes(x = M, y= d, color = M) ) + 
  facet_grid( ~v ) + 
  geom_point( size = 2.5 ) +
  geom_point( size = 9, shape = "_") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "wheat", colour = "black"), 
    panel.background = element_rect(colour = "black") ) + 
  scale_y_continuous(limits = c(0.5, NA), oob = scales::squish) +
  labs(color = "Model", y = "concordance-index", x ="") + 
  scale_color_manual( values = c("hotpink3", "forestgreen") )

gg0<-gg0 + 
  geom_linerange( 
    data = boo, inherit.aes = F,
    aes(ymin = dmin, ymax = dmax, size = p, alpha = p, x= M, color = M)
  ) + 
  scale_alpha_continuous(range = c(0.025, 0.15)) +
  guides(alpha = "none", size = "none")


# difference 
smi = sm[,c("v","contrast","sd","p")]
smi["Hypoxia", "contrast"] = sm["Hypoxia", c("d1","d0")] %>% unlist %>% {.[.<0.5]=0.5;.} %>% diff %>% abs
colnames(smi)[2] = "d"

gg_figure5  = 
gg0 + geom_text(inherit.aes = F, 
            data = smi , 
            aes(y = 0.5025, x = 0.8, label = paste0("p.val = ", p)), hjust = -0.1, size = 2.5)

setwd(wd)
# figure 5
save(file = "figure5.fig", gg_figure5)
