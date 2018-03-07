
#include <R.h>
#include <math.h>

#include "utilities.h"

typedef double (*hmm_meas_dens_t)(double y, double x, double *param);
typedef double (*hmm_trans_eqtn_t)(double x_t, double u_t, double *param);
typedef double (*hmm_init_eqtn_t)(double x_1, double *param);

typedef struct hmm {
  hmm_init_eqtn_t init;
  hmm_meas_dens_t meas;
  hmm_trans_eqtn_t tran;  
} hmm_t;

void sir(double *x_up, double *x_pr, double *m, double *u, int P);

void particle_filter(double **x_pr, double **x_up, double *y, int T, double *param, double **u_pr, double **u_up, double **m, double **tr, int P, hmm_t hmm);

void is_filter(double *loglik, double **is_pr, double **is_up, double **trans, double **meas, int T, int P);

double sv_meas_dens(double y, double x, double *param);

double sv_trans_eqtn(double x_t, double u_t, double *param);

double sv_init(double u_1, double *param);

void sv_filter(int *status, double *sigma2_pr, double *sigma2_up, double *loglik, double *param, double *y, double *u_pr, double *u_up, int *T, int *P);

double ll_meas_dens(double y, double x, double *param);

double ll_trans_eqtn(double x_t, double u_t, double *param);

double ll_init_eqtn(double u_1, double *param);

void ll_filter(int *status, double *_x_pr, double *_x_up, double *_meas, double *_tran, double *loglik, double *param, double *y, double *_u_pr, double *_u_up, int *T, int *P);
