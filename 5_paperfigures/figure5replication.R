library(magrittr)
library(survival)
library(ggplot2)

wd = getwd()
# set wd to here
setwd(glue::glue("{root}/4_replication"))

# load results Shen et al test set1
load("Shen_et_al_testset1.Rdata")
# load results Shen et al test set2
load("Shen_et_al_testset2.Rdata")
# load results Su et al
load("Su_et_al.Rdata")

# function to calculate bootstrapping confidence interval
#' y :outcome
#' yh: predicted scores
#' nBoo: number of boo samples 
#' seed: random seed 
#' probs: confidence envelopes
f_CI <- function(y, yh, nBoo = 200, seed = 42,
                 probs = seq(0.025, 0.475, by = 0.05)){
  
  n = length(y)
  df = data.frame(y=y, yh=yh)
  set.seed(seed)
  dh = replicate(nBoo,{
    inds = sample(n, replace = T)
    concordance(y~yh, df[inds,])$concordance
  })
  dmin = quantile(na.omit(dh), probs)
  dmax = quantile(na.omit(dh), 1-probs)
  
  data.frame(ci = 1-2*probs ,dmax = dmax, dmin = dmin)
}


# Shen et al --------------------------------------------------------------
df1 = lapply(list(df_shen_test1, df_shen_test2), function(df_){
  # calculate CIs and get ggplottable data
  df = lapply(grep("yh", colnames(df_), value = T), function(i){
    y = df_$status
    yh = df_[[i]]
    re = f_CI(y, yh)
    data.frame( re, model = i, cohort = paste0("Shen ", df_$cohort[1]),
                d0 = concordance(y~yh)$concordance )
  }) %>% do.call(what = rbind)
  
  df$model = factor(df$model, 
                    levels = c("yh_shen", "yh_base","yh_mubu"),
                    labels = c("Shen", "base", "our_met6" ) )
  df
}) %>% do.call(what = rbind) 

gg_shen = 
ggplot(df1) +
  geom_linerange( data = df1,
                  aes(x=model, ymin=dmin, ymax=dmax,
                      color=model, size = -ci, alpha = -ci) ) +
  geom_point(data = df1, size = 2.5,
             aes(x=model, y=d0, color=model)) +
  geom_point(data = df1, size = 9, shape = "_",
             aes(x=model, y=d0, color=model)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "wheat", colour = "black"), 
    panel.background = element_rect(colour = "black") ) + 
  labs(color = "Model", y = "concordance-index", x ="") + 
  scale_color_manual( values = c("black","hotpink3", "forestgreen") )+
  scale_alpha_continuous(range = c(0.025, 0.15)) +
  guides(alpha = "none", size = "none") +
  facet_grid(~cohort)



# Su et al  ---------------------------------------------------------------

# calculate CIs and get ggplottable data
df_ = df_su
df_$cohort = "Su"
df2 = lapply(grep("yh", colnames(df_), value = T), function(i){
  y = df_$wost2
  yh = df_[[i]]
  re = f_CI(y, yh)
  data.frame( re, model = i, cohort = df_$cohort[1],
              d0 = concordance(y~yh)$concordance )
}) %>% do.call(what = rbind)

df2$model = factor(df2$model, 
                  levels = c("yh_base","yh_mubu"),
                  labels = c("base", "our_met6" ) )

gg_su = 
ggplot(df2) +
  geom_linerange( data = df2,
                  aes(x=model, ymin=dmin, ymax=dmax,
                      color=model, size = -ci, alpha = -ci) ) +
  geom_point(data = df2, size = 2.5,
             aes(x=model, y=d0, color=model)) +
  geom_point(data = df2, size = 9, shape = "_",
             aes(x=model, y=d0, color=model)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "wheat", colour = "black"), 
    panel.background = element_rect(colour = "black") ) + 
  labs(color = "Model", y = "concordance-index", x ="") + 
  scale_color_manual( values = c("hotpink3", "forestgreen") )+
  scale_alpha_continuous(range = c(0.025, 0.15)) +
  guides(alpha = "none", size = "none") +
  facet_grid(~cohort)


# # load figure 5 data
# load(glue::glue("{root}/3b_prediction/fig5.gg"))
# 
# library(patchwork)
# 
# gg0 = gg0 + scale_y_continuous(limits = c(0.5, 0.8), oob = scales::squish) 
# gg_su = gg_su + scale_y_continuous(limits = c(0.675, 0.875), oob = scales::squish) 
# gg_shen = gg_shen + scale_y_continuous(limits = c(0.5, 1), oob = scales::squish) 
# 
# # pdf 9 X 9
# (gg0 + coord_fixed(ratio = 15))/(
#   (gg_su + coord_fixed(ratio = 22.5)) +
#     (gg_shen + coord_fixed(ratio = 12))
# )

setwd(wd)
# figure 5replication
save(file = "figure5replication.fig", gg_su, gg_shen)

