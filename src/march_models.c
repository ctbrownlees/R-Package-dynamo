
#include "march_models.h"

// Bivariate DCC(1,1) Model Filter
void bidcc_filter(int *status, double *rho, double* eps, double *loglik, double *param, double *_y, int *T){

  int t;
  double logden;
  double alpha,beta;
  double **Q, **y;
  double rho_bar;
  *loglik = 0;

  // sanity check
  if( !finite(param[0]) || !finite(param[1]) ){
	*loglik = -HUGE_VAL;
	return;
  }
  
  alpha = param[0];
  beta  = param[1];
  
  // check constraints
  if( alpha <= 1e-5 || beta < 0 || (alpha+beta)>1 ){
	*loglik = -HUGE_VAL;
	return;
  }

  // allocate
  Q = create_real_matrix(*T,3);
  y = create_and_copy_real_matrix(*T,2,_y);
    
  // init
  rho_bar = 0;
  for( t=0; t<*T; ++t ){ rho_bar += y[t][0]*y[t][1]; }
  rho_bar /= *T;
  
  Q[0][0] = 1;
  Q[0][1] = 1;
  Q[0][2] = rho_bar;
  rho[0]  = rho_bar;
  
  // loop
  *loglik = 0;
  for( t=1; t<*T; ++t ){
    Q[t][0] = (1-alpha-beta)         + alpha*y[t-1][0]*y[t-1][0] + beta*Q[t-1][0];
    Q[t][1] = (1-alpha-beta)         + alpha*y[t-1][1]*y[t-1][1] + beta*Q[t-1][1];
    Q[t][2] = rho_bar*(1-alpha-beta) + alpha*y[t-1][0]*y[t-1][1] + beta*Q[t-1][2];    
    rho[t]  = Q[t][2]/sqrt(Q[t][0]*Q[t][1]);
    
    logden  = -0.5*log(2*PI) - 0.5*log(1-rho[t]*rho[t]) - 0.5*(y[t][0]*y[t][0]+y[t][1]*y[t][1]-2*y[t][0]*y[t][1]*rho[t])/ (1.0-rho[t]*rho[t]);

    if( finite(logden) ){
    	*loglik += logden;    
    }
    else{
	Rprintf("problem at time %d\n",t);
    }

  }
  
  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }
  
  // cleanup
  destroy_real_matrix(Q,*T,3);
  destroy_real_matrix(y,*T,2);
}

// MEWMA Model Filter
void mewma_filter(int *status, double *_s, double* _eps, double *loglik, double *param, double *_y, int *T, int *N){

  int t,i,j;
  double logden;
	double lambda;
	double ***S, **y, **eps;
  double **Sig;
  double *work1, **work2;
  double rho_bar;
  *loglik = 0;
  
  lambda = param[0];
  
  // check constraints
	if( lambda <= 1e-5 || lambda>1 ){
		*loglik = -HUGE_VAL;
		return;
	}

  // allocate
  S     = create_real_array3d(*T,*N,*N);
  y     = create_and_copy_real_matrix(*T,*N,_y);
  eps   = create_and_copy_real_matrix(*T,*N,_eps);
  work1 = create_real_vector(*N);
  work2 = create_real_matrix(*N,*N);

  // init
  for( i=0; i<*N; ++i){
    for( j=0; j<=i; ++j ){
	    work2[i][j]=0;
      for( t=0; t<*T; ++t ){ work2[i][j] += y[t][i]*y[t][j]; }
	    work2[i][j] /= *T;
      work2[j][i] = work2[i][j];
    }
  }
  chol(S[0],work2,*N);
  *loglik = 0;

  // loop
  for( t=1; t<*T; ++t ){

    chol_up(S[t],S[t-1],y[t-1],*N,lambda,1.0-lambda,work1);    
    fwdinv(work2,S[t],*N);
    matvec(eps[t],work2,y[t],*N);
    
    logden = -0.5*(*N)*log(2*PI);
    for(i=0;i<*N;++i) logden += -log( S[t][i][i] )-0.5*eps[t][i]*eps[t][i];
    
    *loglik += logden;
  }
  
  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }
  
  // copy results
  real_array3d_copy(S,*T,*N,*N,_s);
  real_matrix_copy(eps,*T,*N,_eps);
  
  // cleanup
  destroy_real_array3d(S,*T,*N,*N);
  destroy_real_matrix(y,*T,*N);
  destroy_real_matrix(eps,*T,*N);
  destroy_real_vector(work1,*N);
  destroy_real_matrix(work2,*N,*N);

}
