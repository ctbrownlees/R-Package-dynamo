
#include "particle_filtering.h"

void sir(double *x_up, double *x_pr, double *m, double *u, int P) {
    
    int i,j;
    double m_sum, m_l, m_u;
    
    m_sum = 0.0;    
    for (i = 0; i < P; i++) m_sum += m[i];
    
    // init
    j     = 0;
    m_u   = 0;
    
    // resample
    for (i = 0; i < P; i++) {
        m_l = m_u;
        m_u = m_u + m[i];
        
        while( (m_l/m_sum <= u[j]) && (u[j] < m_u/m_sum) ){
            x_up[j] = x_pr[i];
            if (j < P) j++;
            else break;
        }
    } 

}

void particle_filter(double **x_pr, double **x_up, double *y, int T, double *param, double **u_pr, double **u_up, double **m, double **tr, int P, hmm_t hmm){

	int t;
	int p;
  
  for(p=0;p<P;++p){
    x_pr[0][p] = hmm.init(u_pr[0][p],param);
  }

	for(t=0;t<T;++t){

		// filtering
    for( p=0; p<P; ++p) m[t][p] = hmm.meas(y[t],x_pr[t][p],param);
        
    sir(x_up[t],x_pr[t],m[t],u_up[t],P);
    
    // update
    if( t < T-1 ) for(p=0;p<P;++p){ 
        x_pr[t+1][p] = hmm.tran(x_up[t][p],u_pr[t][p],param); 
        tr[t][p] =  x_pr[t+1][p] + x_up[t][p];
    }
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

double sv_init(double u_1, double *param){
	return 0;
}

void sv_filter(int *status, double *_sigma2_pr, double *_sigma2_up, double *loglik, double *param, double *y, double *_u_pr, double *_u_up, int *T, int *P){
	
	double **sigma2_pr, **sigma2_up;
	double **u_up, **u_pr;
  double **m, **tr;

	//
	sigma2_pr = create_real_matrix(*T,*P);
	sigma2_up = create_real_matrix(*T,*P);
  u_pr      = create_and_copy_real_matrix(*T,*P,_u_pr);
  u_up      = create_and_copy_real_matrix(*T,*P,_u_up);
  m         = create_real_matrix(*T,*P);
  tr        = create_real_matrix(*T,*P);

	// 
	//particle_filter(sigma2_pr,sigma2_up,y,*T,param,u_pr,u_up,m,tr,*P,sv_meas_dens,sv_trans_eqtn,sv_init);

	//  
	destroy_real_matrix(sigma2_pr,*T,*P);
	destroy_real_matrix(sigma2_up,*T,*P);
	destroy_real_matrix(u_pr,*T,*P);
  destroy_real_matrix(u_up,*T,*P);
  destroy_real_matrix(m,*T,*P);
  destroy_real_matrix(tr,*T,*P);

}

double ll_init_eqtn(double u_1, double *param){
  double sigma2_x, sigma2_y;
  sigma2_x = param[0];
  sigma2_y = param[1];  

  return sqrt(sigma2_x) * u_1;
}

double ll_meas_dens(double y, double x, double *param){
  double sigma2_x, sigma2_y;
  sigma2_x = param[0];
  sigma2_y = param[1];  
  return sqrt( 1.0 / sqrt( 2 * PI * sigma2_y ) ) * exp( -0.5 * ( (y-x)*(y-x) / sigma2_y  ) );
}

double ll_trans_eqtn(double x_t, double u_t, double *param){
  double sigma2_x, sigma2_y;
  sigma2_x = param[0];
  sigma2_y = param[1];  
	return x_t + sqrt(sigma2_x) * u_t;
}

double ll_trans_dens(double x_aft, double x_bef, double *param){
  double sigma2_x, sigma2_y;
  sigma2_x = param[0];
  sigma2_y = param[1];  
  return sqrt( 1.0 / sqrt( 2 * PI * sigma2_x ) ) * exp( -0.5 * ( (x_aft-x_bef)*(x_aft-x_bef) / sigma2_x  ) );
}

void ll_filter(int *status, double *_x_pr, double *_x_up, double *_meas, double *_tran, double *loglik, double *param, double *y, double *_u_pr, double *_u_up, int *T, int *P){
  
	double **x_pr, **x_up;
	double **u_up, **u_pr;
  double **m, **tr;
  
  hmm_t ll;

  ll.init = ll_init_eqtn;  
  ll.meas = ll_meas_dens;
  ll.tran = ll_trans_eqtn;

	//
	x_pr = create_real_matrix(*T,*P);
	x_up = create_real_matrix(*T,*P);
  u_pr = create_and_copy_real_matrix(*T,*P,_u_pr);
  u_up = create_and_copy_real_matrix(*T,*P,_u_up);
  m    = create_real_matrix(*T,*P);
  tr   = create_real_matrix(*T,*P);

	// 
	particle_filter(x_pr,x_up,y,*T,param,u_pr,u_up,m,tr,*P,ll);

  //
  real_matrix_copy(x_pr,*T,*P,_x_pr);
  real_matrix_copy(x_up,*T,*P,_x_up);
  real_matrix_copy(m,*T,*P,_meas);
  real_matrix_copy(tr,*T,*P,_tran);
  
	//  
	destroy_real_matrix(x_pr,*T,*P);
	destroy_real_matrix(x_up,*T,*P);
	destroy_real_matrix(u_pr,*T,*P);
  destroy_real_matrix(u_up,*T,*P);
  destroy_real_matrix(m,*T,*P);
  destroy_real_matrix(tr,*T,*P);

}
