#include "tweeDEseq.h"

double **logcfactor(int n, double alpha){
  int i, j; 
  double **C;
  
  C = (double **)malloc(n*sizeof(double *));
  for(i=0;i<n;i++)
    C[i] = (double *)malloc((i+1)*sizeof(double));

  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      if(i==j)
	C[i][j] = 0;
      else if(j==0)
	C[i][j] = lgamma(i+1-alpha)-lgamma(1-alpha);
      else
	C[i][j] = logspace_add(C[i-1][j-1],C[i-1][j]+log(i-(j+1)*alpha));
    }
  return C;
}
	  

SEXP nprobs(SEXP N, SEXP ALPHA, SEXP DELTA, SEXP THETA){
  int *n = INTEGER(N), k, i;
  double *alpha = REAL(ALPHA), *delta = REAL(DELTA), *theta = REAL(THETA), *res, aux;
  SEXP RES;

  PROTECT(RES = allocVector(REALSXP,*n+1));
  res = REAL(RES);

  double **C = logcfactor(*n,*alpha);

  double p0 = exp(-(*delta)*(pow(*theta+1,*alpha) - pow(*theta,*alpha))/(*alpha));
  res[0] = p0;

  for(k=1;k<=*n;k++){
    aux = 0;
    for(i=1;i<=k;i++)  
      aux += exp(C[k-1][i-1]-lgamma(k+1))*pow(*delta,i)*pow(*theta+1,i*(*alpha)-k);  
    res[k] = aux*p0;
  }

  for(i=0;i<*n;i++)
    free(C[i]);
  free(C);

  UNPROTECT(1);
  return RES;
}

      
    
