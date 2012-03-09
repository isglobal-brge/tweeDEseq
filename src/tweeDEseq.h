#include <R.h>
#include<Rinternals.h>
#include<Rmath.h>

double var(double* x, int n);
double mean(double* x, int n);
double pval(double* x1, double* x2, int n1, int n2);
int* perm(int* x, int n);
double ttest(double* x, int* class, SEXP M, int n1, int n2);
double **logcfactor(int n, double alpha);
double zhuprobs2(int n, SEXP A, double b, SEXP C, SEXP TOL);
double *vecProd(int n, double* x, double* y);
double vecSum(int n, double* x);
double weightedMean(int n, double* x, double* w);
double *weightedVar(int n, double* x, double* w, double h1);
double k3(int n, double* x);
double kappa3(int n, double* x);
double k(double d, double k3est);


SEXP permtest(SEXP X, SEXP CLASS, SEXP NPERM, SEXP N, SEXP M);
SEXP nprobs(SEXP N, SEXP ALPHA, SEXP DELTA, SEXP THETA);
SEXP zhuprobs(SEXP N, SEXP A, SEXP B, SEXP C, SEXP TOL);
SEXP loglikGlm(SEXP NOBS, SEXP NCOV, SEXP A, SEXP C, SEXP PAR, SEXP X, SEXP Y, 
	       SEXP TOL, SEXP OFFSET);
SEXP cov_wt_C(SEXP X, SEXP Y, SEXP N, SEXP P);
SEXP momentEstimates_wt_C(SEXP X, SEXP Y, SEXP N, SEXP P);
