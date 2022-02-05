# lasso solution is obtained by double-ridge as described by Tibshirani

# keep only, memory efficient, necessary information for a coxme object
ff_koni <- function(acx){
  acx$hmat <- NA
  attr(acx$formulaList$fixed, ".Environment") <- NULL
  # acx$y <- NA # it is fine, just %10 overhead
  acx$terms = NA
  acx
}

# calculate delta between two fits, based on coefs
f_delta <- function(fit0, fit1, dg = 3){
  sqrt(sum((round(fit0$B_samespace, digits = dg) - 
              round(fit1$B_samespace, digits = dg))^2))
}

# obtain lasso solution with double ridge 
# max_iter, previously N, is max of iteration to converge
# delta: eps for convergence
# r_FF: function to decide which model, ordinal regression or logistic regression i,e r_OR or r_LR
f_doubre_ridge_as_lasso_converges <- function(df, outc, mets, vars, day, r_FF, max_iter = 20, delta = 0.01){
  
  # base model 
  ror0 = r_FF(df, outc, mets, vars, day, excludeM = F)
  ror0$B_samespace = ror0$frail$`1`
  # keep memory efficient version
  ror0 = ff_koni(ror0)
  
  # recursively double ridge 
  ff <- function(n, df1, sweepB, flist){
    # if iterations are done 
    if( n > max_iter ) return(flist)
    
    df0 = df1
    df1[,mets] = sweep(df1[,mets], 2, abs(sweepB), `*`)
    ror1 = r_FF(df1, outc, mets, vars, day, excludeM = F)
    # keep memory efficient version
    ror1 = ff_koni(ror1)
    
    # keep same space estimation to compare later 
    ror1$B_samespace = ror1$frail$`1`*abs(sweepB)
    ror1$n_iter = n+1
    
    # cat("\n iteration ", n, ": \n")
    flist[[as.character(n)]] <- ror1
    
    # # verbose
    # print( unlist( VarCorr(ror1) ) )
    # plot(ror0$frail$`1`, ror1$frail$`1`)
    # points(ror0$frail$`1`, ror1$B_samespace, col ="blue")
    
    # if delta achived stop
    di =  f_delta(flist[[length(flist)-1]],  ror1)
    # cat("\n delta: ", di, "\n")
    flist[[length(flist)]]$ldelta = di
    # converged and returned
    if(di <= delta) return(flist)
    
    # factoring into origincal space
    Bsweep_i = exp( rowSums( log( abs( sapply(flist,function(x) x$frail$`1`) ) ) ) )
    ff(n+1, df0, Bsweep_i, flist)
  }
  ror0$n_iter = 0
  ff(0, df, ror0$frail$`1`, list(Base = ror0))
}

# linear mixed effect ordinal regression with L1 penalty 
# model obtained by dobly regularizing ridge of coxme which 
# is optimized with integrated partial likelihood over theta and beta
# l1(lasso) regularized ordinal regression for mixed effect model 
# df: dataframe including variables
# mets: names of metabolites to be shrunk
# vars: fixed effect variables
# day: random effect variable 
# max_iter: max iteration for convergence
# delta: eps for convergence
# return_objs: if return intermediate objects
# r_FF: r_OR or r_LR
f_reor4me <- function(df, outc, mets, vars, day, r_FF, max_iter = 20, delta = 0.01, return_objs = F){
  # stop("this function is under development")
  # # implement some convenience functions like chisq stats, coefs, delta change etc.
  # # ...
  mcall = match.call()
  nods <- f_doubre_ridge_as_lasso_converges(df, outc, mets, vars, day, r_FF = r_FF, max_iter = max_iter, delta = delta)
  
  # chisq stat 
  xsq = sapply( nods, function(x)
    apply(f_cxs(x), 1, function(z) 
      pchisq(z["Chisq"], z["df"], lower.tail = F, log.p = T) )
  )
  
  # converged fit coefficients
  fit = nods[[length(nods)]]
  # class(fit) = c("l1_sorcery",class(fit))
  
  fit$converged <- fit$ldelta <= delta 
  
  coef_l1 <- function(fit){
    B = round(fit$B_samespace, 3)
    B[B!=0] = fit$B_samespace[B!=0]
    # coefficients for penalized and unpenalized together for Lasso
    c(fit$coefficients, B %>% {names(.)=substring(text = names(.),first = 2);. })
  }
  
  # lasso
  B1 = coef_l1(fit)
  # ridge
  B2 = coef(nods[[1]])
  
  # pseudo path for optimized lasso
  B_pseudopath = sapply(nods[-1], coef_l1)
  
  # VarCorr, interpretible effect in coxme
  varcorr = rbind(B2 = unlist(VarCorr(nods[[1]])), 
                  B1 = unlist(VarCorr(fit)))
  
  re <- list(B1 = B1, B2 = B2, B_hpath = B_pseudopath, varcorr = varcorr, xchisq = xsq, 
             nods = if(return_objs) nods else NULL, converged = fit$converged, n_iter = fit$n_iter, call = mcall )
  
  class(re) <- c("sorcery_l12", class(re)) # lasso optimum parameter solution
  
  re
} 

coef.sorcery_l12 <- function(fit, l = 1) if(l==1) fit$B1 else fit$B2 

predict.sorcery_l12 <- function(fit, newx, l = 1){
  B = coef.sorcery_l12(fit, l = l)
  xh = as.matrix(newx[, names(B)])
  xh %*% B
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
