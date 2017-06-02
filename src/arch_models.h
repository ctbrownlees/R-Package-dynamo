
#include <R.h>
#include <math.h>

#include "utilities.h"

// GARCH(1,1) Model Filter
void garch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T);

// TARCH(1,1) Model Filter
void tarch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T);

// APARCH(1,1) Model Filter
void aparch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T);


