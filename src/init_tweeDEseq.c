#include "tweeDEseq.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

R_CallMethodDef CallMethods[]  = {
  {"permtest", (DL_FUNC) &permtest, 5}, //permtest_in},        
  {"nprobs", (DL_FUNC) &nprobs, 4},// nprobs_in}, 
  {"zhuprobs", (DL_FUNC) &zhuprobs, 5},// zhuprobs_in}
  {NULL, NULL, 0}
};


void R_init_tweeDEseq(DllInfo *ddl){

  R_registerRoutines(ddl, NULL, CallMethods, NULL, NULL);
}
