
#include "arch_models.h"

void garch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

	double logden;
	double omega, alpha, beta;
	int t;

	omega = param[0];
	alpha = param[1];
	beta  = param[2];

	// check constraints
	if( alpha <= 1e-6 || beta < 0 || omega<=0 || (alpha+beta)>1 ){
		*loglik = -HUGE_VAL;
		return;
	}

	// init
	sigma2[0]=0;
	for( t=0; t<10; ++t ){ sigma2[0] += y[t]*y[t]; }
	sigma2[0] /= 10;
	eps[0] = y[0]/sqrt( sigma2[0] );
	*loglik = 0;
  
	// loop
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];	
		eps[t]    = y[t]/sqrt( sigma2[t] );

		logden    = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];
		*loglik   += logden;
	}
  
  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

}

// TARCH(1,1) Model Filter
void tarch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

	double logden;
	double omega, alpha, gamma, beta;
	int t;

	omega = param[0];
	alpha = param[1];
	gamma = param[2];
	beta  = param[3];

	// check constraints
	if( alpha <= 0 || beta < 0 || omega<0 || (alpha+beta)>1 ){
		*loglik = -HUGE_VAL;
		return;
	}

	// init
	sigma2[0]=0;
	for( t=0; t<10; ++t ){ sigma2[0] += y[t]*y[t]; }
	sigma2[0] /= 10;
	eps[0] = y[0]/sqrt( sigma2[0] );
  
	// loop 
	*loglik = 0;
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + gamma * y[t-1]*y[t-1]*((double)(y[t-1]<0)) + beta * sigma2[t-1];  
		eps[t]    = y[t]/sqrt( sigma2[t] );

		logden    = -0.5 *log(2*PI) -0.5*log(sigma2[t]) -0.5*(y[t]*y[t])/sigma2[t];
		*loglik   += logden;
	}
  
  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }
}

// APARCH(1,1) Model Filter
void aparch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){
  
  double logden;
  double omega, alpha, gamma, beta, mu, delta;
  int t;
  
  mu    = param[0];
  omega = param[1];
  alpha = param[2];
  gamma = param[3];
  beta  = param[4];
  delta = param[5];
  
  
  // check constraints
  if( alpha <= 0.0 || beta <= 0.0 || omega <= 0.0 || delta <= 0.0 ){ // ) || alpha + beta > 1){
    *loglik = -HUGE_VAL;
    return;
  }
  
  // init
  sigma2[0]=0;
  for( t=0; t<10; ++t ){ sigma2[0] += (y[t]-mu)*(y[t]-mu); }
  sigma2[0] /= 10;
  sigma2[0] = pow(sigma2[0], delta/2 );
  eps[0] = (y[0]-mu)/ pow( sigma2[0], 1/delta );

  // loop
  *loglik = 0;
  logden    = -0.5 *log(2*M_PI) -0.5*(1/delta)*log(sigma2[0]) -0.5*((y[0]-mu)*(y[0]-mu))/pow(sigma2[0], 1/delta);
  *loglik   += logden;
  for( t=1 ; t<*T ; ++t ){
    sigma2[t] = omega + alpha * pow( (fabs(y[t-1] - mu) - gamma * (y[t-1] - mu)) , delta) + beta * sigma2[t-1];
    eps[t]    = (y[t]-mu)/ pow( sigma2[t], 1/delta );
    logden    = -0.5 *log(2*M_PI) - (1/delta)*log(sigma2[t]) -0.5*((y[t]-mu)*(y[t]-mu))/pow(sigma2[t], 2/delta);
    *loglik   += logden;
  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

}

