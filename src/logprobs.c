#include "tweeDEseq.h"

SEXP logprobs(SEXP N, SEXP ALPHA, SEXP DELTA, SEXP THETA, SEXP TOL){
  int *n = INTEGER(N), k, i, b = 0, ix = 0;
  double *alpha = REAL(ALPHA), *delta = REAL(DELTA), *theta = REAL(THETA), *res, aux, *tol = REAL(TOL);
  
  SEXP RES;
  
  PROTECT(RES = allocVector(REALSXP,*n+1));
  res = REAL(RES);
  
  //  double **C = logcfactor(*n,*alpha);
  
  double p0 = -(*delta)*(pow(*theta+1,*alpha) - pow(*theta,*alpha))/(*alpha);
  res[0] = p0;
  
  double ltol = log(*tol);
  
  double *newC = (double *)malloc(1*sizeof(double));
  double *oldC = (double *)malloc(1*sizeof(double));
  newC[0] = 0;
  oldC[0] = 0;
  
  for(k=1;k<=*n;k++){
    aux = 0;
    // ASSIGN NEW POSITION TO THE VECTOR newC
    if(k>1)
      newC = (double*)realloc(newC, k*sizeof(double));
    for(i=1;i<=k;i++){

      // COMPUTACION OF THE NEW C USING THE OLD ONE
      /* if(i == 1) */
      /* 	newC[0] = 0; */
      /* else{ */
      /* 	newC[0] = lgamma(k-*alpha)-lgamma(1-*alpha); */
      /* 	newC[k-1] = 0; */
      /* 	for(j=2; j<k; j++) */
      /* 	  newC[j-1] = logspace_add(oldC[j-2], oldC[j-1]+log((k-1)-j*(*alpha))); */
      /* } */

      if(i == k)
	newC[i-1] = 0;
      else if(i == 1)
	newC[i-1] = lgamma(k-*alpha)-lgamma(1-*alpha);
      else
	newC[i-1] = logspace_add(oldC[i-2], oldC[i-1]+log((k-1)-i*(*alpha)));

      if(i==1)
	aux = newC[i-1] + i*log(*delta) + (i*(*alpha)-k)*log(*theta + 1);
      else
	aux = logspace_add(aux, newC[i-1] + i*log(*delta) + (i*(*alpha)-k)*log(*theta + 1));  
    }
    // ASSIGN NEW POSITION TO THE VECTOR oldC
    oldC = (double*)realloc(oldC, (k+1)*sizeof(double));
    // COPY THE NEW C TO THE OLD ONE
    for(i=0; i<k; i++)
      oldC[i] = newC[i];
    res[k] = aux + p0 - lgamma(k+1);
    // IF PROBABILITIES ARE BELOW A CERTAIN TOLERANCE AND DECAYING STOP COMPUTING THEM AND ASSIGN DIRECTLY 0
    if ( (res[k] < ltol) && ( k > 1 ) ){
      if (res[k] < res[k-1]){
	b = 1;
	ix = k;
	break;
      }
    }
  }
  
  if (b == 1)
    for (i=ix+1; i<=*n; i++)
      res[i] = 0;
  
  free(newC), free(oldC);
  
  UNPROTECT(1);
  return RES;
}

      
    
/* SEXP doublelogprobs(SEXP N, SEXP ALPHA, SEXP DELTA1, SEXP DELTA2, SEXP THETA, SEXP TOL){ */
/*   int *n = INTEGER(N), k, i, b = 0, ix = 0; */
/*   double *alpha = REAL(ALPHA), *delta1 = REAL(DELTA1), *delta2 = REAL(DELTA2), *theta = REAL(THETA), *res, aux1, aux2, *tol = REAL(TOL); */
  
/*   SEXP RES; */
  
/*   PROTECT(RES = allocVector(REALSXP,2*(*n+1))); */
/*   res = REAL(RES); */
  
/*   //  double **C = logcfactor(*n,*alpha); */
  
/*   res[0] = -(*delta1)*(pow(*theta+1,*alpha) - pow(*theta,*alpha))/(*alpha); */
/*   res[*n+1] = -(*delta2)*(pow(*theta+1,*alpha) - pow(*theta,*alpha))/(*alpha); */
  
/*   double ltol = log(*tol); */
  
/*   double *newC = (double *)malloc(1*sizeof(double)); */
/*   double *oldC = (double *)malloc(1*sizeof(double)); */
/*   newC[0] = 0; */
/*   oldC[0] = 0; */
  
/*   for(k=1;k<=*n;k++){ */
/*     aux1 = 0; */
/*     aux2 = 0; */
/*     // ASSIGN NEW POSITION TO THE VECTOR newC */
/*     if(k>1) */
/*       newC = (double*)realloc(newC, k*sizeof(double)); */
/*     for(i=1;i<=k;i++){ */

/*       // COMPUTACION OF THE NEW C USING THE OLD ONE */
/*       /\* if(i == 1) *\/ */
/*       /\* 	newC[0] = 0; *\/ */
/*       /\* else{ *\/ */
/*       /\* 	newC[0] = lgamma(k-*alpha)-lgamma(1-*alpha); *\/ */
/*       /\* 	newC[k-1] = 0; *\/ */
/*       /\* 	for(j=2; j<k; j++) *\/ */
/*       /\* 	  newC[j-1] = logspace_add(oldC[j-2], oldC[j-1]+log((k-1)-j*(*alpha))); *\/ */
/*       /\* } *\/ */

/*       if(i == k) */
/* 	newC[i-1] = 0; */
/*       else if(i == 1) */
/* 	newC[i-1] = lgamma(k-*alpha)-lgamma(1-*alpha); */
/*       else */
/* 	newC[i-1] = logspace_add(oldC[i-2], oldC[i-1]+log((k-1)-i*(*alpha))); */

/*       if(i==1){ */
/* 	aux1 = newC[i-1] + i*log(*delta1) + (i*(*alpha)-k)*log(*theta + 1); */
/* 	aux2 = newC[i-1] + i*log(*delta2) + (i*(*alpha)-k)*log(*theta + 1); */
/*       } */
/*       else{ */
/* 	aux1 = logspace_add(aux1, newC[i-1] + i*log(*delta1) + (i*(*alpha)-k)*log(*theta + 1));   */
/* 	aux2 = logspace_add(aux2, newC[i-1] + i*log(*delta2) + (i*(*alpha)-k)*log(*theta + 1));   */
/*       } */
/*     } */
/*     // ASSIGN NEW POSITION TO THE VECTOR oldC */
/*     oldC = (double*)realloc(oldC, (k+1)*sizeof(double)); */
/*     // COPY THE NEW C TO THE OLD ONE */
/*     for(i=0; i<k; i++) */
/*       oldC[i] = newC[i]; */
/*     res[k] = aux1 + res[0] - lgamma(k+1); */
/*     res[k + *n + 1] = aux2 + res[*n + 1] - lgamma(k+1); */
/*     // IF PROBABILITIES ARE BELOW A CERTAIN TOLERANCE AND DECAYING STOP COMPUTING THEM AND ASSIGN DIRECTLY 0 */
/*     if ( (res[k] < ltol) && ( k > 1 ) && (res[k + *n + 1] < ltol) ){ */
/*       if ( (res[k] < res[k-1]) && (res[k + *n] < res[k + *n + 1]) ){ */
/* 	  b = 1; */
/* 	  ix = k; */
/* 	  break; */
/*       } */
/*     } */
/*   } */
    
/*   if (b == 1){ */
/*     for (i=ix+1; i<=*n; i++){ */
/*       res[i] = 0; */
/*       res[i + *n + 1] = 0; */
/*     } */
/*   } */

/*   free(newC), free(oldC); */
  
/*   UNPROTECT(1); */
/*   return RES; */
/* } */

      
