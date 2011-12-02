testPoissonTweedie <- function(x, group, saveModel=FALSE, ...){

  if (missing(group))
   stop("'group' argument should be provided")

  if (length(x) != length(group))
    stop("'x' and 'group' must have equal length")
  
  groups <- unique(group)
  ngroups <- length(groups)


  if (!ngroups>=2)
   stop("two groups are required")

  if (ngroups>2)
   stop("more than two groups is not yet supported")


  y <- list() 
  mm <- rep(NA, ngroups)
  for (i in 1:ngroups)
   {
    y[[i]] <- x[group==groups[i]]
    mm[i] <- mean(y[[i]])
   } 
  names(mm) <- c("meanA", "meanB")


  mod <- n <- list()
  for (i in 1:ngroups)
   {
    if (!all(y[[i]]==0))
     mod[[i]] <- mlePoissonTweedie(y[[i]], ...)
    else
     mod[[i]] <- NA
   }
 
   if (any(is.na(mod)))
    {
     non.all.zero <- c(1:ngroups)[!is.na(mod)]
     estad <- (mod[[non.all.zero]]$par[1] - 0) / mod[[non.all.zero]]$se[1]
    }
   else
    {
     if (mod[[1]]$convergence==0 && mod[[2]]$convergence==0)
      {
 #       estad <- abs(mod[[1]]$par[1] - mod[[2]]$par[1])/ sqrt(mod[[1]]$se[1]^2 + mod[[2]]$se[1]^2)
       num <- log(mod[[1]]$par[1]) - log(mod[[2]]$par[1])
       den <-  sqrt(mod[[1]]$se[1]^2/mod[[1]]$par[1]^2 + mod[[2]]$se[1]^2/mod[[2]]$par[1]^2)
       estad <- num/den
      }
     else
       estad <- NA
    }
  pval <- 2*pnorm(abs(estad), lower=FALSE)
    
  names(pval) <-  "pval"

  ans <- list(mean=mm, stat=estad, pvalue=pval)

   if (saveModel)
    attr(ans, "model") <- mod
   ans
}


