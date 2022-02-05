# training logistic regression with S in coxph function --------------------
# create log-reg compatible Surv outcome for coxph 
f_lr <- function(n, y){
  t0 = seq(n)/2
  t1 = t0+0.2
  # zero classes
  zrs = cbind( t0, t1+0.05*(1-y), y)
  # surrogate 1 scores: 1/(1+e^x)
  zrs1 = cbind( t0-0.2, t1+0.05*(y), 1-y)
  rbind(zrs,zrs1)
}

# given classes y and data dt, creates log-reg compatible 
# data and Surv outcome for coxph
f2LogReg <- function(y, dt, dkeep=NULL){
  i0 = y==0 # zero class
  y0 =  f_lr(sum(i0), 0) # create cout type survival outcome
  d0 = rbind(dt[i0,], dt[i0,]*0) # shape data accordingly 1/(1+e^x)
  if(!is.null(dkeep)) d0 = cbind(d0, rbind(dkeep[i0,,drop=F], dkeep[i0,,drop=F]))
  
  i1 = y==1 # one class
  y1 =  f_lr(sum(i1), 1)
  d1 = rbind(dt[i1,], dt[i1,]*0) # e^x/(1+e^x)
  if(!is.null(dkeep)) d1 = cbind(d1, rbind(dkeep[i1,,drop=F], dkeep[i1,,drop=F]))
  
  yy = rbind(y0, y1)
  dd = rbind(d0, d1)
  zz = c(y0[,1]*0, y1[,1]*0+1)
  list(S = Surv(yy[,1], yy[,2], yy[,3]), dd = dd, stz = zz)
}


# removing ties to get concordance/ sommers'D optimizer directly ----------

# function to melt ties
# t0: base start
f_mt <- function(n, t0=0){
  # if no ties
  if(n<2) return( cbind(0, t0+0.4, 1) )
  
  k = rep(0.4+t0, n)
  # first all be in until melting time comes 
  bas =  cbind(0, k, 0)
  # batir cikar, each comes in and out per time 
  sug = 
    cbind(t0+seq(k)/2,
          t0+seq(k)/2 + 0.4,1)
  rbind(bas, sug)
}

# function, removing the ties for given 
# r: survival times ranking, 
#    high(first low risk) to low(last high risk) (will be reverted)
# df: data to be reconstructed accordingly
f2NoTies <- function(r, df){
  
  # ranks, death first, and space allocation for melting
  r = rank(-r, ties.method = "min", na.last = "keep")
  df = df[order(r),] 
  r = r[order(r)] 
  
  tb = table(r)
  rr = cumsum(c(0,tb)) %>% {.[-length(tb)]} 
  names(rr)<- names(tb)
  # count type survival data
  Sc = apply( cbind(start = rr,n = as.numeric(tb)),1, function(x) f_mt(x[2],x[1])) %>%
    do.call(what = rbind)
  
  Sc = Surv(Sc[,1], Sc[,2], Sc[,3])
  # df without ties
  tdf<-
    lapply(unique(r), function(i){
      if(sum(r==i)<2) df[r==i, ,drop = F] else
        rbind( df[r==i, ,drop = F], df[r==i, ,drop = F])
    }) %>% do.call(what = rbind)
  
  list(Sc = Sc, tdf = tdf)
}
