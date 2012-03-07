logLik.glmPT <- function(object, ...)
  object$value

AIC.glmPT <- function(object, ...){
  df <- object$df
  ans <- -2*logLik(object) + 2*df
  ans
}
                      


print.glmPT <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  coef <- x$par
  coef <- round(coef,digits=digits)
  if (length(coef)) {
    d <- length(coef)-2
    coef[d+1] <- paste(c(rep(" ",10),coef[d+1]),collapse="")
    coef[d+2] <- paste(c(rep(" ",10),coef[d+2]),collapse="")
    cat("Coefficients:\n")
    print.default(format(coef[1:d], digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\nPoisson-Tweedie parameters:\n")
    print.default(format(coef[(d+1):length(x$par)], digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  invisible(x)
}

summary.glmPT <- function(object, ...){
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(object$par)) {
    cat("Coefficients:\n")
    d <- length(object$par)-2
    coefs <- object$par[1:d]
    t <- coefs/object$se[1:d]
    pv <- 1 - pnorm(t)
    signif.stars <- any(pv < 0.1)
    Signif <- symnum(pv, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    Cf <- array("", dim = c(d,5), dimnames = list(names(object$par)[1:d],c("Estimate","Std.Error","t value","Pr(>|t|)","")))
    Cf[,1] <- round(coefs,4)
    Cf[,2] <- round(object$se[1:d],4)
    Cf[,3] <- round(t,4)
    Cf[,4] <- format.pval(pv)
    Cf[,5] <- Signif
    print.default(Cf, quote = FALSE, right = TRUE, ...)
    cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    cat("\nPoisson-Tweedie parameters\n")
    PTCf <- array("", dim = c(2,2), dimnames = list(c("c","a"),c("Estimate","Std.Error")))
    PTCf[,1] <- round(object$par[(d+1):(d+2)],4)
    PTCf[,2] <- round(object$se[(d+1):(d+2)],4)
    print.default(PTCf,quote=FALSE,right=TRUE, ...)
    invisible(list(model.coef=Cf, pt.coef=PTCf))
  }
  else
    cat("No coefficients\n\n")
}

anova.glmPT <- function(object, modelNull, ...){
  if(missing(modelNull))
    modelNull <- update(object, ~1)
  D <- -2*modelNull$value +2*object$value
  df <- object$df - modelNull$df
  if(D>=0)
    pval <- pchisq(q = D, df = df, lower.tail = FALSE)
  else
    pval <- 1
  pval
}


glmPT.fit <- function(X, Y, offset=NULL, allFactors=FALSE, a = NULL, ...){
  ncov <- ncol(X)
  if (is.null(a)){
    par.ini <- c(log(mean(Y)), rep(0,ncov-1), 0.9, 0)
    lower <- c(0,rep(-Inf,ncov-1),1e-5, -Inf)
    upper <- c(rep(Inf,ncov),1 - 1e-3, 1)
    mle <- tryCatch(optim(par.ini, loglikGlmPT, X=X, Y=Y, offset=offset, allFactors=allFactors, ..., method="L-BFGS-B", lower=lower, upper=upper, hessian=TRUE, control = list(fnscale = -1, maxit=1e3)), warning = function(w) w)
    mle$ncov <- ncov
  }
  else{
    par.ini <- c(log(mean(Y)), rep(0,ncov-1), 0.9)
    lower <- c(0, rep(-Inf,ncov-1), 1e-5)
    upper <- c(rep(Inf,ncov), 1-1e-3)
    mle <- tryCatch(optim(par.ini, loglikGlmPT, X=X, Y=Y, offset=offset, allFactors=allFactors, a=a, ..., method="L-BFGS-B", lower=lower, upper=upper, hessian=TRUE, control = list(fnscale = -1, maxit=1e3)), warning = function(w) w)
    mle$ncov <- ncov
  }
  mle$df <- length(mle$par)
  mle
}
  

glmPT <- function(formula, data, offset=NULL, a=NULL, ...)
 {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf, contrasts)
  ncov <- ncol(X)
  datClass <- attr(mt,"dataClasses")
  covClass <- datClass[2:length(datClass)]
  if(any(covClass!="factor"))
    allFactors = FALSE
  else
    allFactors = TRUE
  mle <- glmPT.fit(X,Y,offset=offset, allFactors=allFactors, a, ...)
  se <- try(sqrt(-diag(solve(mle$hessian))), TRUE)
  if( !inherits(se, "try-error") )
    mle$se <- se
  else
    mle$se <- rep(NA, nrow(mle$hessian))
  mle$call <- cl
  mle$contrasts = attr(X, "contrasts")
  mle$fitted.values <- as.numeric(exp(X%*%mle$par[1:ncov]))
  mle$residuals <- Y - mle$fitted.values
  mle$ncov <- ncov - 1
#  mle$df <- length(mle$par)
  if(length(mle$par)!=ncov+2){
    mle$par[length(mle$par)+1] <- a
    mle$se[length(mle$se)+1] <- NA
  }
  names(mle$par) <- c(colnames(X), colnames(Y),"c","a")
  class(mle) <- "glmPT"
  mle
}
