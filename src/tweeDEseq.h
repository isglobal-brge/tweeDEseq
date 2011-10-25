#include <R.h>
#include<Rinternals.h>
#include<Rmath.h>

double var(double* x, int n);
double mean(double* x, int n);
double pval(double* x1, double* x2, int n1, int n2);
int* perm(int* x, int n);
double ttest(double* x, int* class, SEXP M, int n1, int n2);
double **logcfactor(int n, double alpha);

SEXP permtest(SEXP X, SEXP CLASS, SEXP NPERM, SEXP N, SEXP M);
SEXP nprobs(SEXP N, SEXP ALPHA, SEXP DELTA, SEXP THETA);
SEXP zhuprobs(SEXP N, SEXP A, SEXP B, SEXP C, SEXP TOL);
