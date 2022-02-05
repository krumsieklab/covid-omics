# this script so as to create figure5 supplement based on 
# precalculated lasso path objects 
library(magrittr)
library(survival)
library(coxme)


wd = getwd()
# set working directory to location of source code
setwd(glue::glue("{root}/3b_prediction"))
load("lasso_path_loo_fits.Rdata")

# # regularized linear mixed effect model
# functions needed to cast problem to be solved with cox(ph)/(me)
source("./model/utils_coxph_sorceries.R")
# ridge regularization 
source("./model/utils_coxme_sorceries.R")
# lasso regularization
source("./model/utils_lasso_sorceries.R")
# utils for lasso regularization path
source("./model/utils_lasso_path_sorceries.R")
# utility function to calculate c-index for loo
source("./model/utils_loo_prediction.R")

# concatanete constituents of composite outcome
colnames(compo)[ colnames(compo)%in%c("intubation1", "max_aki_stage" ) ] = c("intubation", "max_aki")
df = dplyr::left_join(df, compo[,-ncol(compo)], by = c("ID"="Patient_ID")) 

# which n feas can be feasible to show 
sapply(nods_lpath, function(x) unique(colSums(x$fit$B_path!=0))) %>% 
  unlist %>% table %>% {as.numeric(names(.[.>50]))} %>% sort

# collect prediction performance for each point on the path 
d_train = calculate_d_lpath_loo(nods_lpath, df, "r0", "training", 0)
d_test = calculate_d_lpath_loo(nods_lpath, df, "r0", "test", 0)

# compile all path results per outcome ------------------------------------
d_all<- cbind(rbind(cbind(d_train, dset = "training"), 
                    cbind(d_test, dset = "test")), 
              outc = "Composite Outcome")

df_all = d_all 

gg_lpath<-
  ggplot(df_all, aes(x = nf-3, y = d, 
                     ymin = d_lb, 
                     ymax = d_ub, 
                     color = outc,
                     fill = outc,
                     lty = dset)) +
  geom_ribbon(color = NA, alpha = 0.1) +
  geom_line() + 
  # geom_vline(xintercept = 32, lty = 2) +
  facet_wrap(outc~., nrow = 1, dir = "v", )+
  theme_minimal() + 
  theme(strip.background = element_rect(fill = "wheat", colour = "wheat")) + 
  labs(color = "outcome",fill = "outcome", 
       y = "C-index", 
       x ="number of metabolites selected", 
       lty = "dataset") +
  theme(legend.position = "bottom"#, 
        #strip.text = element_text(size = 12),
        #axis.text.x.bottom =element_text(angle = 45, hjust = 1, vjust = 0.5)
  ) +
  scale_color_viridis_d(option="magma", end = 0.85) +
  scale_fill_viridis_d(option="magma", end = 0.85)
# save(file = "gg_lpath.Rdata", gg_lpath)

gg_lpath = 
gg_lpath + 
  geom_vline(xintercept = 32, color = "red", lty =2 , lwd = 0.5)+
  geom_vline(xintercept = 6, color = "blue", lty =2 , lwd = 0.5) +
  ylab("c-index")

gg_lpath + 
  geom_vline(xintercept = 32, color = "red", lty =2 , lwd = 0.5)+
  geom_vline(xintercept = 6, color = "blue", lty =2 , lwd = 0.5)+
  xlim(c(0, 40)) + theme(legend.position = "n") + ylab("C-index")

setwd(wd)
# figure 5replication
save(file = "figure5Supplement.fig", gg_lpath)

