#include "tweeDEseq.h"

//UNBIASED VARIANCE COMPUTATION (x: obs.vector, n: number of obs.) 
double var(double* x, int n){
  int i;
  double aux = 0, aux2 = 0;

  for (i=0;i<n;i++){
    aux += x[i];
    aux2 += x[i]*x[i];
  }
  
  double m = (double)n;
  double res = (1/(m-1))*aux2 - (1/(m*m-m))*aux*aux;
  return res;
}

//MEAN COMPUTATION (x: obs.vector, n: number of obs.)         
double mean(double* x, int n){
  int i;
  double res = 0;
  
  for(i=0;i<n;i++)
    res += x[i];

  return res/n;
}

//STATISTIC COMPUTATION
double pval(double* x1, double* x2, int n1, int n2){
  double varx1 = var(x1,n1), varx2 = var(x2,n2), meanx1 = mean(x1,n1), meanx2 = mean(x2,n2);

  double m1 = (double)n1, m2 = (double)n2;
  double pooled = sqrt(varx1/m1 + varx2/m2);

  double tstat = (meanx1 - meanx2)/pooled;

  double res;
  if(tstat<0)
   res = tstat;
  else
   res = -tstat; 
}

//RANDOM PERMUTATION OF VECTOR x
int* perm(int* x, int n){
  int i, rnd;
  int *res = x;
  int aux;

  GetRNGstate();
  for(i=0;i<(n-1);i++){
    rnd = floor(runif(0,n-i));
    aux = res[i]; res[i] = res[i+rnd]; res[i+rnd] = aux;
  }
  PutRNGstate();  
  
  return res;
}

double ttest(double* x, int* class, SEXP M, int n1, int n2){
  int i, cont1=0, cont2=0, *m = INTEGER(M);
  double *x1, *x2;
  
  x1 = (double *)malloc(n1*sizeof(double));
  x2 = (double *)malloc(n2*sizeof(double));
  
  for (i=0;i<*m;i++){
    if(class[i] == 0){
      x1[cont1] = x[i];
      cont1++;
    }
    if(class[i] == 1){
      x2[cont2] = x[i];
      cont2++;
    }
  }
    
  double pvalue = pval(x1,x2,n1,n2);
  free(x1);
  free(x2);
  return pvalue;
}
  
  
//PERMUTATION TEST. X: DATA MATRIX, CLASS: VECTOR OF GROUPS(CASE-CONTROL), NPERM: NUMBER OF WISHED PERMUTATIONS, 
//N: NUMBER OF ROWS IN X, M: NUMBER OF COLUMNS IN X
SEXP permtest(SEXP X, SEXP CLASS, SEXP NPERM, SEXP N, SEXP M){
  int *n = INTEGER(N), *m = INTEGER(M), *nperm = INTEGER(NPERM), *class = INTEGER(CLASS);
  int i, j, *pclass, n1 = 0, n2 = 0;
  double *x;

  for(i = 0;i<*m;i++){
    if(class[i] == 0)
      n1++;
    else
      n2++;
  }

  double *pvalreals, **pvalperms;
  pvalreals = (double *)malloc((*n)*sizeof(double));
  pvalperms = (double **)malloc((*nperm)*sizeof(double *));
  for(i=0;i<*nperm;i++)
    pvalperms[i] = (double *)malloc((*n)*sizeof(double));


  for(i=0;i<*n;i++){
    if(*n==1)
      x = REAL(X);
    else
      x = REAL(VECTOR_ELT(X,i));
    pvalreals[i] = ttest(x,class,M,n1,n2);
  }
  
  for (i=0;i<*nperm;i++){
    pclass = perm(class,*m);
    for(j=0;j<*n;j++){
      if(*n==1)
	x = REAL(X);
      else
	x = REAL(VECTOR_ELT(X,j));
      pvalperms[i][j] = ttest(x,pclass,M,n1,n2);
    }
  }
  
  SEXP RES;
  double *res;
  PROTECT(RES = allocVector(REALSXP,*n));
  res = REAL(RES);

  int cont;
  for(i=0;i<*n;i++){
    cont=0;
    for(j=0;j<*nperm;j++){
      if(pvalreals[i]>pvalperms[j][i])
	cont++;
    }
    res[i] = cont/(double)(*nperm);
  }

  UNPROTECT(1);
  return RES;
}
  
