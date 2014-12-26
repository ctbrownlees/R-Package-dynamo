
#include <R.h>
#include <math.h>

// GARCH(1,1) Model Filter
void garch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

	double logden;
	double omega, alpha, beta;
  int t;

	omega = param[0];
	alpha = param[1];
	beta  = param[2];
  
  if( alpha <= 0 || beta < 0 || omega<0 || (alpha+beta)>1 ){
    *loglik = -10000000;
    return;
  }

  *loglik = 0;

	// init sigma2
	sigma2[0]=0;
	for( t=0; t<10; ++t ){
		sigma2[0] += y[t]*y[t];
	}
	sigma2[0] /= 10;
  
  //
	eps[0] = y[0]/sqrt( sigma2[0] );
  
	// garch loop 
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];	
		eps[t]    = y[t]/sqrt( sigma2[t] );
		logden    = log( sigma2[t] ) -0.5 * (eps[t]*eps[t])/sigma2[t];
    
    *loglik   += logden;
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
  
  if( alpha <= 0 || beta < 0 || omega<0 || (alpha+0.5*gamma+beta)>1 ){
    *loglik = -10000000;
    return;
  }

  *loglik = 0;

	// init sigma2
	sigma2[0]=0;
	for( t=0; t<10; ++t ){
		sigma2[0] += y[t]*y[t];
	}
	sigma2[0] /= 10;
  
  //
	eps[0] = y[0]/sqrt( sigma2[0] );
  
	// garch loop 
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + gamma * y[t-1]*y[t-1]*((double)(y[t-1]<0)) + beta * sigma2[t-1];	
		eps[t]    = y[t]/sqrt( sigma2[t] );
		logden    = log( sigma2[t] ) -0.5 * (eps[t]*eps[t])/sigma2[t];
    
    *loglik   += logden;
	}

}