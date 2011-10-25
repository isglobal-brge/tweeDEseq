#include "tweeDEseq.h"

SEXP zhuprobs(SEXP N, SEXP A, SEXP B, SEXP C, SEXP TOL){
  int *n = INTEGER(N), i, j, nbreak=*n+1;
  double *a = REAL(A), *b = REAL(B), *c = REAL(C), *tol = REAL(TOL), *res, *r, aux;
  SEXP RES;

  PROTECT(RES = allocVector(REALSXP,*n+1));
  res = REAL(RES);

  if(*a==0)
    res[0] = pow(1-*c,*b);
  else
    res[0] =  exp((*b)*(pow(1-*c,*a)-1)/(*a));

  if(*n!=0){
    aux = (*b)*(*c);
    r = (double *)malloc((*n)*sizeof(double));
    r[0] = (1-*a)*(*c);
    for(i=1;i<*n;i++)
      r[i] = (*c)*r[i-1]*(i-1+*a)/(i+1);

    res[1] = aux*res[0];

    for(i=1;i<*n;i++){
      res[i+1] = aux*res[i];
      for(j=1;j<=i;j++)
        res[i+1] += j*r[i-j]*res[j];
      res[i+1] /= (i+1);
      if((res[i+1]<=*tol)&&(res[i+1]<res[i])){
        nbreak = i;
        break;
      }
    }

    for(i=nbreak+1;i<=*n;i++)
      res[i] = 0;
  free(r);
  }

  UNPROTECT(1);
  return RES;
}
