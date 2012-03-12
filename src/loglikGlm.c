#include "tweeDEseq.h"

SEXP loglikGlm(SEXP NOBS, SEXP NCOV, SEXP A, SEXP C, SEXP PAR, SEXP X, SEXP Y, 
	       SEXP TOL, SEXP OFFSET, SEXP MAXCOUNT){
  
  int *nobs = INTEGER(NOBS), *ncov = INTEGER(NCOV), *y = INTEGER(Y), *maxCount = INTEGER(MAXCOUNT);
  double *a = REAL(A), *c = REAL(C), *par = REAL(PAR), *offset = REAL(OFFSET);

  int i = 0, j = 0;
  double m, b, probs, *loglik, *x;
  SEXP LOGLIK;
  
  PROTECT(LOGLIK = allocVector(REALSXP, 1));
  loglik = REAL(LOGLIK);
  *loglik = 0;

  for(i = 0; i < *nobs; i++){
    x = REAL(VECTOR_ELT(X, i));
    m = offset[i];
    for(j = 0; j < *ncov; j++)
      m += x[j]*par[j];
    m = exp(m);
    b = (m * pow(1 - *c, 1 - *a)) / (*c);
    if(b<0.001)
      b = 0.001;
    if((*a<0.001)&&(-(*a)<0.001)){
      b = (m * (1-*c))/(*c);
      probs = dnbinom_mu((double)y[i], b, m, 1);
    }
    else if(*a<=0.999){
      if(y[i]>*maxCount){
	b = (m * (1-*c))/(*c);
	probs = dnbinom_mu((double)y[i], b, m, 1);
      }
      else{
	probs = zhuprobs2(y[i], A, b, C, TOL);
	probs = log(probs);
      }
    }
    else
      probs = dpois((double)y[i], m, 1);

    if(probs!=probs)
      *loglik -= 1000;
    else
      *loglik += probs;
    //    Rprintf("loglik=%g, probsi=%g\n", *loglik, probs);
  }

  if(*loglik==*loglik+1)
    *loglik = -pow(10,50);

  //  Rprintf("%g with a=%g, m=%g, b=%g, c=%g\n", *loglik, *a, m, b, *c);

  UNPROTECT(1);
  return(LOGLIK);
}
