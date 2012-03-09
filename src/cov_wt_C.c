#include "tweeDEseq.h"

double *vecProd(int n, double* x, double* y){
  int i = 0;
  double *res;  
  res = (double *)malloc(n*sizeof(double));
  
  for (i = 0; i<n; i++)
    res[i] = x[i]*y[i];
  
  return(res);
}

double vecSum(int n, double* x){
  int i = 0;
  double res=0;

  for (i = 0; i<n; i++)
    res += x[i];
  
  return(res);
}

double weightedMean(int n, double* x, double* w){
  double *aux1, aux2;
  aux1 = vecProd(n, x, w);
  aux2 = vecSum(n, aux1);
  free(aux1);
  return(aux2);
}

double *weightedVar(int n, double* x, double* w, double h1){
  int i = 0;
  double *xc, *resWV, *aux;
  xc = (double *)malloc((n)*sizeof(double));
  resWV = (double *)malloc(2*sizeof(double));
  resWV[0] = weightedMean(n, x, w);
  for(i = 0; i<n; i++)
    xc[i] = x[i] - resWV[0];
  aux = vecProd(n, xc, xc);
  resWV[1] = weightedMean(n, aux, w);
  resWV[1] *= h1;
  free(xc), free(aux);
  return(resWV);
}

SEXP cov_wt_C(SEXP X, SEXP Y, SEXP N, SEXP P){
  int i=0, *n = INTEGER(N), *p = INTEGER(P);
  double *x, *z = REAL(Y), *res, *aux, *wvar, *y, sumY, h1;
  SEXP RES;

  y = (double *)malloc((*n)*sizeof(double));
  res = (double *)malloc((*p*2)*sizeof(double));

  // Check whether sum(Y)=1
  sumY = vecSum(*n, z);
  if(sumY != 1){
    // If not, make it happen
    for (i = 0; i<*n; i++)
      y[i] = z[i]/sumY;
  }
  else
    y[i] = z[i];

  aux = vecProd(*n, y, y);
  h1 = 1/(1 - vecSum(*n, aux));

  //  PROTECT(RES = allocVector(VECSXP, *p));
  PROTECT(RES = allocMatrix(REALSXP, 2, *p));
  res = REAL(RES);

  for(i = 0; i<*p; i++){
    x = REAL(VECTOR_ELT(X, i));
    wvar = weightedVar(*n, x, y, h1);
    res[2*i] = wvar[0];
    res[2*i + 1] = wvar[1];
  }

  free(aux), free(wvar), free(y);

  UNPROTECT(1);
  return(RES);
}

double k3(int n, double* x){
  int i = 0;
  double mm, *resk3, aux;
  resk3 = (double *)malloc(n*sizeof(double));
  
  mm = vecSum(n, x)/n;
  for(i=0;i<n;i++){
    aux = x[i] - mm;
    resk3[i] = aux*aux*aux;
  }
  mm = vecSum(n, resk3)/n;
  
  free(resk3);
  return(mm);
}

double kappa3(int n, double* x){
  double reskappa3 = k3(n, x);
  double aux = vecSum(n, x)/n;
  reskappa3 = reskappa3/aux - 1;
  return(reskappa3);
}

double k(double d, double k3est){
  double ans = (k3est - 3*(d-1));
  ans = ans/((d-1)*(d-1));
  return(ans);
}

SEXP momentEstimates_wt_C(SEXP X, SEXP Y, SEXP N, SEXP P){
  int *n = INTEGER(N), *p = INTEGER(P), i = 0;
  double *res, *x, *moments, d, kappa, resk, a, *z = REAL(Y), *y, sumY, *aux, h1;
  SEXP RES;

  y = (double *)malloc((*n)*sizeof(double));

  // Check whether sum(Y)=1
  sumY = vecSum(*n, z);
  if(sumY != 1){
    // If not, make it happen
    for (i = 0; i<*n; i++)
      y[i] = z[i]/sumY;
  }
  else
    y[i] = z[i];
  
  aux = vecProd(*n, y, y);
  h1 = 1/(1 - vecSum(*n, aux));

  if(*p>1)
    PROTECT(RES = allocMatrix(REALSXP, 5, *p));
  else
    PROTECT(RES = allocVector(REALSXP, 5));
  res = REAL(RES);

  for(i=0;i<*p;i++){
    x = REAL(VECTOR_ELT(X,i));
    moments = weightedVar(*n, x, y, h1);
    d = moments[1]/moments[0];
    kappa = kappa3(*n, x);
    resk = k(d, kappa);
    a = (resk - 2)/(resk - 1);
    res[5*i] = moments[0];
    res[5*i + 1] = moments[1];
    res[5*i + 2] = d;
    res[5*i + 3] = resk;
    res[5*i + 4] = a;
  }

  free(y), free(moments), free(aux);
 
  UNPROTECT(1);
  return(RES);
}
