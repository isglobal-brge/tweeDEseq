momentEstimates <- function(x, w = NULL){
   if(is.null(w)){
     mu <- mean(x)
     d <- var(x)/mean(x)
     kappa <- kappa3(x)
     k <- k(d, kappa)
     a <- (k-2)/(k-1)
     ans <- c(mu, d, k, a)
     names(ans) <- c("mu", "D", "k", "a")
   }
   else{
     if(length(x)!=length(w))
       stop("'x' and 'w' have different lengths")
     X <- list(as.numeric(x))
     Y <- as.numeric(w)
     p <- 1L
     n <- as.integer(length(x))
     ansC <- .Call("momentEstimates_wt_C", X, Y, n, p)
     ans <- c("mu" = ansC[1], "D" = ansC[3], "k" = ansC[4], "a" = ansC[5])
   }
   ans
 }

k3<-function(x)
 {
  mm <- mean(x)
  num <- (x-mm)^3
  ans <- sum(num)/length(x)
  ans
 }

kappa3 <- function(x)
 {
  ans <- k3(x)/mean(x) - 1
  ans
 }

k <- function(d, k3.est)
 {
   ans <- (k3.est - 3*(d-1))/(d-1)^2
   ans 
 }

kvector <- function(x)
 {
   k(x[1], x[2])
 }
