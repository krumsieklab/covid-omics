# obtain lasso solution with double ridge, with full path
# pmf: penalty multiplication factor, to get sparser solution
f_doubre_ridge_as_lasso_path <- function(df, outc, mets, vars, day, r_FF, max_iter, pmf, vf0 =  NULL){
  
  # variance factor to get doble ridge as lasso
  get_vf <- function(fit, pmf){
    vf = as.list( unname( unlist( VarCorr(fit) ) ) )
    # print(vf %>% unlist)
    vf[[1]] = vf[[1]]*pmf
    # print(vf %>% unlist)
    vf
  }
  
  # recursively double ridge 
  ff <- function(n, df1, sweepB, flist, vf){
    # if iterations are done 
    if(n > max_iter) return(flist)
    
    df0 = df1
    df1[,mets] = sweep(df1[,mets], 2, abs(sweepB), `*`)
    ror1 = r_FF(df1, outc, mets, vars, day, excludeM = F, vfixed = vf)
    
    # keep same space estimation to compare later 
    ror1$B_samespace = ror1$frail$`1`*abs(sweepB)
    
    ror1 = ff_koni(ror1)
    flist[[as.character(n)]] <- ror1
    
    # cat("\n iteration ", n, ": \n")
    # print( unlist( VarCorr(ror1) ) )
    # plot(ror0$frail$`1`, ror1$frail$`1`)
    # points(ror0$frail$`1`, ror1$B_samespace, col ="blue")
    
    Bsweep_i = sapply(flist,function(x) x$frail$`1`) %>% abs %>% log %>% rowSums %>% exp
    ff(n+1, df0, Bsweep_i, flist, vf = get_vf(ror1, pmf))
  }
  
  # base model i.e. ridge
  ror_ridge = if(!is.list(vf0)) r_FF( df, outc, mets, vars, day, excludeM = F ) else
    r_FF( df, outc, mets, vars, day, excludeM = F, vfixed = vf0 )
  
  ror_ridge = ff_koni(ror_ridge)
  ror_ridge$B_samespace = ror_ridge$frail$`1`
  
  ff(0, df, ror_ridge$frail$`1`, list(Base = ror_ridge))
  
}

f_reor4me_path <- function(df, outc, mets, vars, day, r_FF, N_iter = 20, pmf =0.9, vf0 = NULL,  return_objs = F){
  # stop("this function is under development")
  # # implement some convenience functions like chisq stats, coefs, delta change etc.
  # # ...
  mcall = match.call()
  nods <- f_doubre_ridge_as_lasso_path(df, outc, mets, vars, day, r_FF, N_iter, pmf, vf0)
  
  # chisq stat 
  xsq = sapply( nods, function(x)
    apply(f_cxs(x), 1, function(z) 
      pchisq(z["Chisq"], z["df"], lower.tail = F, log.p = T) )
  )
  
  # # converged fit coefficients
  # fit = nods[[length(nods)]]
  # # class(fit) = c("l1_sorcery",class(fit))
  # 
  # fit$converged <- fit$ldelta <= delta 
  
  coef_l1 <- function(fit){
    B = round(fit$B_samespace, 3)
    B[B!=0] = fit$B_samespace[B!=0]
    # coefficients for penalized and unpenalized together for Lasso
    c(fit$coefficients, B %>% {names(.)=substring(text = names(.),first = 2);. })
  }
  
  # # lasso
  # B1 = coef_l1(fit)
  # ridge
  B2 = coef(nods[[1]])
  
  # pseudo path for optimized lasso
  B_pseudopath = sapply(nods[-1], coef_l1)
  
  # # VarCorr, interpretible effect in coxme
  lvarcorr = sapply(nods, function(x) unlist(VarCorr(x)))
  
  re <- list(
    B_path = B_pseudopath, 
    B2 = B2, 
    xchisq = xsq,
    lvarcorr = lvarcorr,
    nods = if(return_objs) nods else NULL,
    #converged = fit$converged, 
    #n_iter = fit$n_iter, 
    call = mcall 
  )
  
  class(re) <- c("sorcery_paths", class(re)) # lasso optimum parameter solution
  
  re
} 

# keep only, memory efficient, necessary information for a coxme object
ff_koni <- function(acx){
  acx$hmat <- NA
  attr(acx$formulaList$fixed, ".Environment") <- NULL
  # acx$y <- NA # it is fine, just %10 overhead
  acx$terms = NA
  acx
}

# get coxme stats
f_cxs <-  function(x, rcoef=FALSE, digits=options()$digits, ...) {
  # cat("Cox mixed-effects model fit by maximum likelihood\n")
  # if (!is.null(x$call$data)) 
  #   cat("  Data:", deparse(x$call$data))
  # if(!is.null(x$call$subset)) {
  #   cat(";  Subset:", deparse(x$call$subset), sep="\n")
  # }
  # else cat("\n")
  
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail<- nrow(x$var) - nvar
  
  omit <- x$na.action
  # cat("  events, n = ", x$n[1], ', ', x$n[2], sep='')
  if(length(omit))
    cat(" (", naprint(omit), ")", sep = "")
  loglik <- x$loglik + c(0,0, x$penalty)
  temp <- matrix(loglik, nrow=1)
  # cat("\n  Iterations=", x$iter, "\n")
  dimnames(temp) <- list("Log-likelihood", 
                         c("NULL", "Integrated", "Fitted"))
  # print(temp)
  # cat("\n")
  chi1 <- 2*diff(x$loglik[c(1,2)]) 
  
  
  chi1 <- 2*diff(loglik[1:2]) 
  chi2 <- 2*diff(loglik[c(1,3)])
  temp <- rbind(c(round(chi1,2), round(x$df[1],2),
                  signif(1- pchisq(chi1,x$df[1]),5),
                  round(chi1- 2*x$df[1],2),
                  round(chi1- log(x$n[1])*x$df[1],2)),
                c(round(chi2,2), round(x$df[2],2),
                  signif(1- pchisq(chi2,x$df[2]),5),
                  round(chi2- 2*x$df[2],2),
                  round(chi2- log(x$n[1])*x$df[2],2)))
  dimnames(temp) <- list(c("Integrated loglik", " Penalized loglik"),
                         c("Chisq", "df", "p", "AIC", "BIC"))
  temp
}

coef.sorcery_paths <- function(fit) cbind(Base = fit$B2, fit$B_path)

predict.sorcery_paths <- function(fit, newx){
  Bpath = coef.sorcery_paths(fit)
  xh = as.matrix(newx[, rownames(Bpath)])
  xh %*% Bpath
}

