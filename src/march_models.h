
#include <R.h>
#include <math.h>

#include "utilities.h"

// Bivariate DCC(1,1) Model Filter
void bidcc_filter(int *status, double *rho, double* eps, double *loglik, double *param, double *_y, int *T);

// MEWMA Model Filter
void mewma_filter(int *status, double *_s, double* _eps, double *loglik, double *param, double *_y, int *T, int *N);

