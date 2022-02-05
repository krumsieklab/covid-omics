

# # ridge-regularized-linear-mixed-effect ordinal and logistic regression --------

# this is for ordinal outcome
# roof function to run ordinal regression with coxme model i.e.
# regularized linear mixed effect prop odds ordinal regression
# outc: outcome names, ranking: death high ranking, mild low ranking
# mets: metabolites to be ridge regularized 
# vars: fixed effect vars, no-regularization
# day: variable to adjust for random effect intercept
# if excludeM, then only vars, ridge part is discarded
# vfixed: fixed variance platforms 
r_OR<- function(df, outc, mets, vars, day, excludeM = F, vfixed = NULL){
  
  nod = f2NoTies(df[,outc], df[,c(mets,vars,day)])
  
  # construct the data frame
  adf = data.frame(S= nod$Sc,
                   Time = nod$tdf[,day],
                   nod$tdf[, vars, drop = F])
  
  # part to be ridge regularized
  adf$M <- as.matrix(nod$tdf[,mets,drop=F])
  
  fr <- if(!excludeM) paste(c("S ~ (M|1) + (1|Time)", vars),collapse = " + ") else 
    paste(c("S ~ (1|Time)", vars),collapse = " + ") 
  # print(fr)
  
  cx = coxme(as.formula(fr), adf, vfixed = vfixed)
  cx$excludeM = excludeM
  class(cx) = c("sorcery",class(cx)) # sorcery object
  cx
}


# this is for binary outcome 
# roof function to run LR with coxme model i.e.
# regularized linear mixed effect logistic regression
# outc: outcome names
# mets: metabolites to be ridge regularized 
# vars: fixed effect vars, no-regularization
# day: variable to adjust for random effect intercept
# if excludeM, then only vars, ridge part is discarded
r_LR<- function(df, outc, mets, vars, day, excludeM = F){
  
  dt = df[,c(vars, day, mets)] %>% 
    {.[,day] =as.numeric(factor(.[,day]));. } %>% 
    as.matrix 
  
  # outcome
  y = df[,outc]
  
  nod = f2LogReg(y, dt)
  
  # construct the data frame
  adf = data.frame(S= nod$S,
                   Time = nod$dd[, day, drop = F],
                   z =  nod$stz,  # latent strata for coxph
                   nod$dd[, vars, drop = F])
  
  # part to be ridge regularized
  adf$M <- as.matrix(nod$dd[,mets,drop=F])
  
  # fr <- paste(c("S ~ strata(z) + (M|1) + (1|Time)", vars),collapse = " + ")
  fr <- if(!excludeM) paste(c("S ~ strata(z) + (M|1) + (1|Time)", vars),collapse = " + ") else 
    paste(c("S ~ strata(z) + (1|Time)", vars),collapse = " + ") 
  # print(fr)
  
  cx = coxme(as.formula(fr), adf)
  cx$excludeM = excludeM
  class(cx) = c("sorcery",class(cx)) # sorcery object
  cx
}

coef.sorcery <- function(fit) 
  if(fit$excludeM) fit$coefficients else
    c(fit$coefficients, fit$frail$`1` %>% 
        {names(.)=substring(text = names(.),first = 2);. })

predict.sorcery <- function(fit, newx){
  B = coef.sorcery(fit)
  xh = as.matrix(newx[, names(B)])
  xh %*% B
}


