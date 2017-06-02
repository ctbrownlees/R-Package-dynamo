
#include <R.h>
#include <math.h>

#include "utilities.h"

typedef double (*hmm_meas_dens_t)(double y, double x, double *param);
typedef double (*hmm_trans_eqtn_t)(double x_t, double u_t, double *param);
typedef double (*hmm_init_t)(double x_1, double *param);

void particle_filter(double **x_pr, double **x_up, double *y, int T, double *param, double **u_pr, double **u_up, int P, hmm_meas_dens_t meas, hmm_trans_eqtn_t trans, hmm_init_t init);

void is_filter(double *loglik, double **is_pr, double **is_up, double **trans, double **meas, int T, int P);

double sv_meas_dens(double y, double x, double *param);

double sv_trans_eqtn(double x_t, double u_t, double *param);

double sv_init(double x_1, double *param);

void sv_filter(int *status, double *sigma2_pr, double *sigma2_up, double *loglik, double *param, double *y, double *u_pr, double *u_up, int *T, int *P);
