# simplify multiple measurement risk score estimation
# alternative to weighting, getting max risk of repeated measures
frog <- function(id, y, yh, f = max, add_id = F){
  bb = sapply( unique(id), function(k) c( (y[id==k])[1], f(yh[id==k]) ) )
  dd = data.frame( y = bb[1,], yh=bb[2,])
  dd$id = if(add_id) unique(id) else NULL
  dd
}

get_c_repeated <- function(id, y, yh, f = max, return_fit = F){
  fit = concordance(y~yh, frog(id,y,yh,f))
  if(return_fit) return(fit)
  fit$concordance
}

# get prediction accuracy for lasso path, this is utility function for loo runs
calculate_d_lpath_loo <- function(nods_lpath, df, yname="r0", dset = c("training","test"), nf_th =0, f= max){
  if(dset == "training"){
    # get all predictions
    mdf = seq(nods_lpath) %>% lapply(function(i){
      yh = predict(nods_lpath[[i]]$fit, df[nods_lpath[[i]]$trids,]) 
      # keep 1 prediction for each number of features
      nf = c( length(nods_lpath[[i]]$fit$B2), colSums( nods_lpath[[i]]$fit$B_path != 0 ))
      j = duplicated(nf)
      yh = yh[,!j]
      nf = nf[!j]
      # we are only interested in training set here
      # yh[!nods_lpath[[i]]$trids,] = NA # see above
      r0 = df[[yname]][nods_lpath[[i]]$trids]
      isNA = is.na(r0)
      list( df = data.frame(i=i, id = df$ID[nods_lpath[[i]]$trids], y = r0)[!isNA, ], yh = yh[!isNA,, drop = F], nf = nf)
    })
    
    # nfs where lasso path is going to be plotted
    rnf = lapply(mdf, `[[`, "nf") %>% unlist %>% table %>% {as.numeric(names(.[.>nf_th]))} %>% sort
    
    # collect all of them 
    rops = lapply(mdf, function(rop){
      # rop = mdf[[2]]
      # print("a")
      colnames(rop$yh) = rop$nf
      names(rop$nf) = rop$nf
      nf = structure( rop$nf[as.character(rnf)], names = as.character(rnf))
      yh = matrix(nrow = nrow(rop$yh), ncol = length(rnf), dimnames = list(NULL,rnf))
      # results available for nfs
      ii = intersect(colnames(rop$yh), colnames(yh))
      yh[,ii] = rop$yh[,ii]
      list(df= rop$df, yh = yh)
    })
    
    # check rnf
    rowSums( sapply(rops, function(x) colSums( is.na( x$yh ) ) ==0 )*1  )
    
    # df for concordance
    ds = structure(as.character(rnf), names = rnf) %>% sapply(function(j){
      # j = as.character(rnf[5])
      dfx = 
        do.call(rbind, lapply(rops, function(x){
          # x = rops[[1]]
          # j = as.character(rnf[5])
          dfj = frog( id = x$df$id, y = x$df$y, yh = x$yh[,j],add_id = T, f = f)
          dfj$i = x$df$i[1]
          dfj
        }))
      # remove of no results for some nfs for some training set
      dfx = na.omit(dfx)
      
      # # it is working perfectly thumbs up
      # concordance(y~yh+strata(i)+cluster(id), dfx)
      # concordance(y~yh+strata(i)+cluster(id), dfx[dfx$i==1,])
      
      dfit = concordance(rank(y)~yh+strata(i)+cluster(id), dfx)
      c(d = dfit$concordance, se = sqrt(dfit$var))
    }) 
    
    return(
      data.frame(d_lb = ds["d",]+qnorm(0.025)*ds["se",],
                 d = ds["d",],
                 d_ub = ds["d",]+qnorm(0.975)*ds["se",],
                 nf = rnf)
    )
  }
  
  if(dset == "test"){
    # get all predictions for test
    mdf_test = seq(nods_lpath) %>% lapply(function(i){
      yh = predict( nods_lpath[[i]]$fit, df ) 
      # keep 1 prediction for each number of features
      nf = c( length(nods_lpath[[i]]$fit$B2), colSums( nods_lpath[[i]]$fit$B_path != 0 ))
      j = duplicated(nf)
      yh = yh[,!j]
      nf = nf[!j]
      # we are only interested in training set here
      # yh[!nods_lpath[[i]]$trids,] = NA # see above
      r0 = df[[yname]]
      
      list( df = data.frame(i=i, id = df$ID, y = r0), 
            yh = yh, nf = nf, 
            trids = nods_lpath[[i]]$trids )
    })
    
    # nfs where lasso path is going to be plotted
    rnf = lapply(mdf_test, `[[`, "nf") %>% unlist %>% table %>% {as.numeric(names(.[.>nf_th]))} %>% sort
    
    # collect all of them 
    rops_test = lapply( mdf_test, function(rop){
      # rop = mdf[[2]]
      # print("a")
      colnames(rop$yh) = rop$nf
      names(rop$nf) = rop$nf
      nf = structure( rop$nf[as.character(rnf)], names = as.character(rnf))
      yh = matrix(nrow = nrow(rop$yh), ncol = length(rnf), dimnames = list(NULL,rnf))
      # results available for nfs
      ii = intersect(colnames(rop$yh), colnames(yh))
      yh[,ii] = rop$yh[,ii]
      list(df= rop$df[!rop$trids,,drop = F], yh = yh[!rop$trids,,drop = F])
    })
    
    rops_test = list( 
      df = do.call(rbind, lapply(rops_test, function(x) x$df)),
      yh = do.call(rbind, lapply(rops_test, function(x) x$yh)) %>% # smooth spline for missing values
        apply(1, function(x) zoo::na.spline(zoo::zoo(x= x,order.by = rnf)) %>% as.numeric() ) %>% t
    )
    colnames(rops_test$yh) = rnf
    sapply(rops_test, dim)
    
    # test ds 
    ds_test = structure( as.character(rnf), names = rnf ) %>% sapply( function(j){
      
      re<- try({
        dfx = frog( id = rops_test$df$id, y = rops_test$df$y, yh = rops_test$yh[,j], add_id = T, f = f)
        
        # # remove of no results for some nfs for some training set
        # dfx = na.omit(dfx)
        
        dfit = concordance(rank(y)~yh+cluster(id), dfx)
        c(d = dfit$concordance, se = sqrt(dfit$var))
      })
      if(inherits(re,"try-error")) return(c(d = NA, se = NA))
      re
    }) 
    
    
    return(
      data.frame(d_lb = ds_test["d",]+qnorm(0.025)*ds_test["se",],
                 d = ds_test["d",],
                 d_ub = ds_test["d",]+qnorm(0.975)*ds_test["se",],
                 nf = rnf)
    )
  }
  
}

# this can be use non-path other predictions as well
# old name for this function: calculate_d_base_loo
calculate_d_loo <- function(nods_base_loo, df, yname="r0", dset = c("training","test"), f = max, return_data = F){
  
  if(dset == "training"){
    # base train set 
    dfx_train = do.call(rbind, lapply(seq(nods_base_loo), function(i){
      x= nods_base_loo[[i]]
      aa = data.frame(id = df$ID, y = df[[yname]], yh = predict( x$fit, df))
      aa = aa[x$trids,]
      data.frame(frog(id = aa$id, y = aa$y, yh = aa$yh, add_id = T, f = f),i=i)
    })) %>% na.omit
    
    dim(dfx_train)
    
    d0_train = concordance(y~yh+strata(i)+cluster(id), dfx_train)
    if(return_data) return(list( dfx = dfx_train ,fit = d0_train))
    return(
      data.frame(d_lb = d0_train$concordance + qnorm(0.025)*sqrt(d0_train$var),
                 d = d0_train$concordance,
                 d_ub = d0_train$concordance+qnorm(0.975)*sqrt(d0_train$var))#,
      #nf = rnf)
    )
  }
  
  if(dset == "test"){
    # browser()
    dfx = do.call(rbind,lapply(nods_base_loo, function(x){
      aa = data.frame(id = df$ID, y = df[[yname]], yh = predict( x$fit, df))
      aa[x$trids,] = NA
      aa
    })) %>% na.omit
    
    
    d0 = get_c_repeated(dfx$id, dfx$y, dfx$yh, return_fit = T, f = f)
    if(return_data) return(list( dfx = dfx ,fit = d0))
    return(
      data.frame(d_lb = d0$concordance + qnorm(0.025)*sqrt(d0$var),
                 d = d0$concordance,
                 d_ub = d0$concordance+qnorm(0.975)*sqrt(d0$var))#,
      #nf = rnf)
    )
  }
}


