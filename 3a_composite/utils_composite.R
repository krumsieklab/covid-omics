
# function to infer rankings given the variables of interests as sm 
# columns are of importance order of variables 

f_infer_rankings_simple <- function(sm){
  # conposite modeling starts here
  A = t(combn(nrow(sm),2))
  adf =  sign( sm[A[,1],] - sm[A[,2],] )
  # how many resolution we can get per comparison
  dresolution = rowSums( abs(adf), na.rm = T) %>% table
  
  # deciding order based on preference
  # assuming variables are ordered based on the importance
  ff <- function(x){
    x = na.omit(x)
    x = x[x!=0]
    if(length(x)<1) return(0)
    x[1]
  }
  
  # resolve the pairwise rankings based on preference order
  # latent outcome for optimization
  a = apply(adf,1,ff)
  
  # which pairs coming from where
  w = apply(adf,1,function(x){
    x = na.omit(x)
    x = x[x!=0]
    if(length(x)<1) return(NA)
    names(x)[1]
  })
  w = factor(w, levels = rev(colnames(sm))) 
  
  
  # estimated contribution of each variable to final comparison scheme
  contribution = factor(w,levels = colnames(adf)) %>% 
    {table(.,useNA = "ifany")} %>% 
    sqrt %>% {./sum(.)} %>% 
    {round(.,4)} %>% {cbind(`%`=.)}
  
  # create surrogate data to find patients rankings
  x = matrix(0, nrow = nrow(sm), ncol = nrow(sm))
  diag(x) = 1
  x = x[A[,1],] - x[A[,2],]
  
  list(rs =  colSums( sweep(x,1,a,`*`)[a!=0,] ), vars_contributions = contribution, dresolution = dresolution)
}

# function to return rankings from preferred order of variables 
# accepts missing values and optimizes the pairwise -> ranking 
#
f_infer_rankings_from_preferences <- function(sm, weighting = F){
  # conposite modeling starts here
  A = t(combn(nrow(sm),2))
  adf =  sign( sm[A[,1],] - sm[A[,2],] )
  # how many resolution we can get per comparison
  dresolution = rowSums( abs(adf), na.rm = T) %>% table
  
  # deciding order based on preference
  # assuming variables are ordered based on the importance
  ff <- function(x){
    x = na.omit(x)
    x = x[x!=0]
    if(length(x)<1) return(0)
    x[1]
  }
  
  # resolve the pairwise rankings based on preference order
  # latent outcome for optimization
  a = apply(adf,1,ff)
  
  # which pairs coming from where
  w = apply(adf,1,function(x){
    x = na.omit(x)
    x = x[x!=0]
    if(length(x)<1) return(NA)
    names(x)[1]
  })
  w = factor(w, levels = rev(colnames(sm))) 
  
  if(weighting){
    ww = as.numeric( w)
    table(ww,useNA = "ifany")
    # preference weights * upsampling, because of missing values
    ww = 2^(ww-1) * sqrt(min(table(w))/table(w)[as.character(w)]) 
  }
  
  # estimated contribution of each variable to final comparison scheme
  contribution = factor(w,levels = colnames(adf)) %>% 
    {table(.,useNA = "ifany")} %>% 
    sqrt %>% {./sum(.)} %>% 
    {round(.,4)} %>% {cbind(`%`=.)}
  
  # create surrogate data to find patients rankings
  x = matrix(0, nrow = nrow(sm), ncol = nrow(sm))
  diag(x) = 1
  x = x[A[,1],] - x[A[,2],]
  
  fitb= glm(y~.+0, family = binomial(),
            data = data.frame(y=a*0, sweep(x,1,a,`*`))[a!=0,], 
            weights = if(weighting) ww[a!=0] else NULL)
  
  list(fit = fitb, vars_contributions = contribution, dresolution = dresolution)
  
}

# sharpen pairwise ranking to ranking by expected ordering 
# conditioned on one before 
f_sharpen_p2r <- function(sm0, yh){
  
  # set missing values
  library(rms)
  # sm = df[,-1]
  # rownames(sm) = df$sample_id
  # sm = sm[sm$death<1,-c(1,2)]
  
  # keep only needed variables
  sm = cbind(sm0, yh = yh)

  if(any(is.na( sm[,"yh"] ))){
    # variables to be used to fill yh in 
    ii = names(which(colSums(is.na(sm[is.na(sm[,"yh"]), -ncol(sm),drop = F])) == 0 ))
    jj = names(which(colSums(!is.na(sm[!is.na(sm[,"yh"]), -ncol(sm),drop = F]))>10))
    ij = intersect(ii, jj)
    print(ij)
    
    dk0 = data.frame(yh = sm[,"yh"],sm[,ij])
    # browser()
    return(dk0)
    fit = orm(yh~., dk0)
    # browser()
    # sapply( na.omit(dk0), var ) %>% print
    sm[is.na(sm[,"yh"]), "yh"] = round( predict(fit, dk0[is.na(dk0$yh),,drop = F], type = "mean"),2 )
    
    rm(dk0)
  }
  
  for( i in rev(colnames(sm))[-1]){
    # i = "pof_support"  #"intubation"  #"renal_replacement" #"O2_device"  
    # "Fi_O2"  #"max_aki" #"hospitalization" #"disposition"
    vars = rev(colnames(sm))
    print(i)
    
    j = vars[which(vars == i)-1] #vars[seq(which(vars == i)-1)]
    x = sm
    xdf = data.frame(y = as.numeric(sm[,i]), x)
    fm =  as.formula(paste0("y~yh+", paste(j,collapse = "+")))
    # fm =  as.formula(y~yh)
    print(fm)
    
    if(length(unique(na.omit(sm[,i])))>2){
      # input based on death and aki
      if(any(is.na(sm[,i]))){
        # print(round( predict(orm(fm,xdf),xdf[is.na(xdf$y), ], type = "mean") ))
        try({sm[is.na(sm[,i]), i] = round(predict(orm(fm,xdf),xdf[is.na(xdf$y), ], type = "mean"),2) })
      }
    }else{
      if(any(is.na(sm[,i]))){
        # print(round( predict(orm(fm,xdf),xdf[is.na(xdf$y), ], type = "mean") ))
        # browser()
        try({
          sm[is.na(sm[,i]), i] = round( 
            predict(glm(fm,xdf[,c("y",j)], family = binomial()), xdf[is.na(xdf$y), ], 
                    type = "response"),3)
        })
      }
    }
  }
  
  # order with preference
  f_op <- function(x){
    x = data.frame(apply(x,2,rank))
    # cardinalities 
    ks = sapply(x, function(a) length(unique(a)))
    d = max(ks)+1
    # number of variables 
    n = length(ks)
    
    r0 = x[,1]
    x = x[,-1]
    while(1){
      r0 = rank(r0*d + x[,1])
      if(ncol(x)==1) break
      x = x[,-1, drop = F]
    }
    r0
    # 
    # structure( rowSums( sapply(seq(n), function(i) d^(n-i)*x[,i] )), d=d, n=n)
  }
  
  
  f_op(sm)
  
}

