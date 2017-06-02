
#include "particle_filtering.h"

void particle_filter(double **x_pr, double **x_up, double *y, int T, double *param, double **u_pr, double **u_up, int P, hmm_meas_dens_t meas, hmm_trans_eqtn_t trans, hmm_init_t init){

	int t;
	int p;

	for(t=1;t<T;++t){

		// prediction

		for(p=0;p<P;++p){
			x_pr[t][p] = trans(x_up[t-1][p],u_pr[t][p],param);
		}

		// filtering
	}

}

double sv_meas_dens(double y, double x, double *param){
	return 0;
}

double sv_trans_eqtn(double x_t, double u_t, double *param){
  double alpha0, alpha1, tau2;
  alpha0 = param[0];
  alpha1 = param[1];
  tau2   = param[2];
	return exp( alpha0 + alpha1*log(x_t) + sqrt(tau2)*u_t );
}

double sv_init(double x_1, double *param){
	return 0;
}

void sv_filter(int *status, double *_sigma2_pr, double *_sigma2_up, double *loglik, double *param, double *y, double *_u_pr, double *_u_up, int *T, int *P){
	
	double **sigma2_pr, **sigma2_up;
	double **u_up, **u_pr;

	//
	sigma2_pr = create_real_matrix(*T,*P);
	sigma2_up = create_real_matrix(*T,*P);
  u_pr      = create_and_copy_real_matrix(*T,*P,_u_pr);
  u_up      = create_and_copy_real_matrix(*T,*P,_u_up);

	// 
	particle_filter(sigma2_pr,sigma2_up,y,*T,param,u_pr,u_up,*P,sv_meas_dens,sv_trans_eqtn,sv_init);

	//  
	destroy_real_matrix(sigma2_pr,*T,*P);
	destroy_real_matrix(sigma2_up,*T,*P);
	destroy_real_matrix(u_pr,*T,*P);
  destroy_real_matrix(u_up,*T,*P);

}
