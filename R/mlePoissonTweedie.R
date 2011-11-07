mlePoissonTweedie <- function (x, a, D.ini, a.ini, maxit = 100, loglik=TRUE, maxCount=20000, ...) 
{
    x <- as.numeric(x)
    if (all(x == 0)) 
        warning("All data are 0's")
    if (any(is.na(x))) {
        cat("There are NA's. They have been removed")
        x <- x[!is.na(x)]
    }
    mu <- mean(x)

    if (max(x)<=maxCount) {
      if (!missing(a)) {
        if (missing(D.ini)) 
            p.ini <- c(var(x)/mean(x))
        else p.ini <- c(D.ini)

        if (p.ini[1]<1)
          p.ini[1] <- 1
        
        MLE <- optim(p.ini, loglikPoissonTweedie2, x = x, a = a, 
            method = "L-BFGS-B", mu = mu, lower = c(1), upper = c(Inf), 
            control = list(fnscale = -1, maxit = maxit), hessian = TRUE, 
            ...)
        p <- MLE$par
        D <- p[1]
        c <- (D - 1)/(D - a)
        b <- (mu * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a))
        se <- c(sqrt(-diag(solve(MLE$hessian))), NA)
      }
      else {

            if (missing(D.ini) || missing(a.ini))
                p.ini <- c(var(x)/mean(x), 0)
            else p.ini <- c(D.ini, a.ini)

            if (p.ini[1]<1)
              p.ini[1] <- 1

            MLE <- try(optim(p.ini, loglikPoissonTweedie, x = x,
                method = "L-BFGS-B", mu = mu, lower = c(1, -Inf),
                upper = c(Inf, 1 - (1e-09)), control = list(fnscale = -1,
                  maxit = maxit), hessian = TRUE, ...), TRUE)
            if (!inherits(MLE, "try-error"))  {
                  p <- MLE$par
                  D <- p[1]
                  a <- p[2]
            }
            else  {
                MLE <- optim(p.ini, loglikPoissonTweedie2, x = x,
                  a = 0, method = "L-BFGS-B", mu = mu, lower = c(1),
                  upper = c(Inf), control = list(fnscale = -1,
                    maxit = maxit), hessian = TRUE, ...)

                 p <- MLE$par
                 D <- p[1]
                 a <- 0
            }
            c <- (D - 1)/(D - a)
              
            b <- (mu * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a))
            se <- try(sqrt(-diag(solve(MLE$hessian))), TRUE)
            if (inherits(se, "try-error"))
              se <- rep(NA, nrow(MLE$hessian))
      }

     se <- c(sqrt((D * mu)/length(x)), se)
     loglik <- MLE$value
     iter <- MLE$counts[1]

     convergence <- MLE$convergence
     message <- MLE$message
     if (convergence!=0)
      {
       warning(paste("algorithm did not converge:", message))

          param <- momentEstimates(x)
          mu <- param[1]
          D <- param[2]
          a <- param[4]
          if (a > 1 || is.na(a))
             a <-NA

#           {
#               # degenerate distribution
#               MLE <- optim(p.ini, loglikPoissonTweedie2, x = x,
#                  a = -1, method = "L-BFGS-B", mu = mu, lower = c(0),
#                  upper = c(Inf), control = list(fnscale = -1,
#                    maxit = maxit), hessian = TRUE, ...)
#               a <- -1
#               D <- MLE$par[1]
#             }    
   
          c <- (D-1)/(D-a)
          b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))
          se <- c(sqrt((D*mu)/length(x)), NA, NA)
 
          
          if (loglik && (!any(is.na(c(mu, D, a)))) && D!=0)
           loglik <- loglikPoissonTweedie(c(D,a), x, mu=mu)
          else
           loglik <- NA

          iter <- message <- NA
          convergence <- 0


      }
     method <- "MLE"
    }
    else {
          param <- momentEstimates(x)
          mu <- param[1] 
          D <- param[2]
          a <- param[4]
          if (a > 1 || is.na(a))
           a <- NA
          c <- (D-1)/(D-a)
          b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))
          se <- c(sqrt((D*mu)/length(x)), NA, NA)
          if (loglik && (!any(is.na(c(mu, D, a)))))
           loglik <- loglikPoissonTweedie(c(D,a), x, mu=mu) 
          else
           loglik <- NA

          iter <- message <- NA 
          convergence <- 0
          method <- "Moments"
    }
   
    parameters <- c(mu, D, a)
    gamma.num <- a^2 * c^2 - 3 * a * c + c + 1
    gamma.den <- sqrt(b * c * (1 - c)^a * (1 - a * c)^3)
    gamma <- gamma.num/gamma.den

    names(parameters) <- names(se) <- c("mu", "D", "a")
   
    paramZhu <- c(a, b, c)
    names(paramZhu) <- c("a", "b", "c")
    names(gamma) <- "skewness"
    alpha <- a
    delta <- b * (c)^a
    theta <- (1 - c)/c
    paramHou <- c(alpha, delta, theta)
    names(paramHou) <- c("alpha", "delta", "theta")
    ans <- list(par = parameters, se = se, loglik = loglik, iter = iter, 
        paramZhu = paramZhu, paramHou = paramHou, skewness = gamma, 
        x = x, convergence = convergence, message=message, method=method)
    class(ans) <- "mlePT"
    ans
}
