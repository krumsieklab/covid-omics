library(survival)
library(magrittr)

# test-set1 data extracted from the paper
# https://www.cell.com/action/showPdf?pii=S0092-8674%2820%2930627-9
# figure 2D, and figure 1A
df1 <- 
  data.frame(
    yh  = seq(10), # rank of estimated risk scores
    y   = c(rep(0,4), 1,0,1,0,1,1), # outcomes 
    age = c(56,60,44,45,62,43,59,43,65,55),
    sample_id = paste0("XG", c(24,20,23,21,45,22,46,25,44,43))
  )

# test-set2 data extracted from the paper
# https://www.cell.com/action/showPdf?pii=S0092-8674%2820%2930627-9
# figure 2E, and figure 1A
df2 <- 
  data.frame(
    yh = seq(19), # rank of estimated risk scores
    y  = c(rep(0,9), 1,0,0,rep(1,6),0), # outcomes 
    age = c(37,18,54,48,37,55,36,51,47,68,46,39,51,33,77,44,56,65,66),
    sample_id = paste0("X2-",c(11,3,9,21,23,2,8,20,28,18,24,26,13,16,19,12,14,7,22))
  )

# bootstrappings based CI
f_CIboo <- function(y, yh,  n = 200){
  set.seed(42)
  n2 = length(y)
  df2 = data.frame(y = y, yh = yh)
  dq2_boo = replicate(concordance(y~yh, df2[sample(n2, replace = T),])$concordance, n = n)
  quantile( na.omit(dq2_boo), probs = c(0.025, 1-0.025))
}

# CIs of model and age for test sets
ci1 = f_CIboo(df1$y, df1$yh)
ci1_age = f_CIboo(df1$y, df1$age)

ci2 = f_CIboo(df2$y, df2$yh)
ci2_age = f_CIboo(df2$y, df2$age)

# gather all stats
dff<-
  data.frame(
    d = c(
      concordance(y~yh,  df1)$concordance,
      concordance(y~age, df1)$concordance,
      concordance(y~yh,  df2)$concordance,
      concordance(y~age, df2)$concordance
    ),
    rbind(ci1, ci1_age, ci2, ci2_age) %>% {colnames(.)=c("lb","ub");.},
    model = c("met&pro","age","met&pro","age"),
    test_set = c("test1(n=10)","test1(n=10)", "test2(n=19)", "test2(n=19)") 
  )

library(ggplot2)
dff %>%
  ggplot(aes(color = model,y = model, x = d)) +
  geom_errorbar(aes( xmin = lb, xmax = ub ), size = 1.25) +
  geom_point(size = 5) + 
  theme_minimal(14) + 
  labs(y = "model", x = "AUC") +
  facet_wrap(test_set~.,ncol = 1) +
  theme(strip.background = element_rect(fill = "wheat"), 
        panel.background = element_rect(colour = "black")) +
  labs(title = "Models developed in cell paper predicts the risk not better than age predicts",
       subtitle = "95% CIs are estimated with bootstrapping")


dff$test_set = factor(dff$test_set, labels = c("C2(n=10)", "C3(n=19)"))

dff %>%
  ggplot(aes(color = model,y = model, x = d)) +
  geom_errorbar(aes( xmin = lb, xmax = ub ), size = 1.25) +
  geom_point(size = 5) + 
  theme_minimal(24) + 
  labs(y = "model", x = "AUC") +
  facet_wrap(test_set~.,ncol = 1) +
  theme(strip.background = element_rect(fill = "wheat"), 
        panel.background = element_rect(colour = "black")) #+
# labs(title = "Models developed in cell paper predicts the risk not better than age predicts",
#      subtitle = "95% CIs are estimated with bootstrapping")

dff %>% {colnames(.)[1:3] = c("AUC", "0.025 CI", "0.975 CI"); . } %>% 
  {.[1:4,1:3] = round(dff[1:4,1:3],digits = 4);.[,c(4,1:3,5)]} 




