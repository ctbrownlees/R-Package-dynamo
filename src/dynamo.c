
#include <R.h>
#include <math.h>

// Gaussian GARCH(1,1) Filter
void garch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

	double logden;
	double omega, alpha, beta;
	int t;

	omega = param[0];
	alpha = param[1];
	beta  = param[2];

	// check constraints
	if( alpha <= 0 || beta < 0 || omega<0 || (alpha+beta)>1 ){
		*loglik = -HUGE_VAL;
		return;
	}

	// init
	sigma2[0]=0;
	for( t=0; t<10; ++t ){
		sigma2[0] += y[t]*y[t];
	}
	sigma2[0] /= 10;
  eps[0] = y[0]/sqrt( sigma2[0] );
  
	// loop 
	*loglik = 0;
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];	
		eps[t]    = y[t]/sqrt( sigma2[t] );

		logden    = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];
		*loglik   += logden;
	}
}

// Gaussian GARCH(1,1) numerical OPG
void garch_opg_num(int *status, double *OPG, double *param, double *y, int *T, double *eps){

  double logden_plus, logden_minus;
  double param_plus[4], param_minus[4];
  double dlogden[4];
	double omega, alpha, beta;
	int t,i,j,p;

  /*
  // init
	sigma2[0]=0;
	for( t=0; t<10; ++t ){
		sigma2[0] += y[t]*y[t];
	}
	sigma2[0] /= 10;
  
	// loop 
	for( t=1; t<*T; ++t ){
    for( p=0; p<4; ++p ){
      
      memcpy(param_plus ,param,4*sizeof(double));
      param_minus[p] += eps[p];
      omega = param_plus[0]; alpha = param_plus[1]; beta  = param_plus[2];    
		  sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];	
		  logden_plus = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];

      memcpy(param_minus,param,4*sizeof(double));
      param_plus[p] += eps[p];
      omega = param_minus[0]; alpha = param_minus[1]; beta  = param_minus[2];    
  	  sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];	
		  logden_minus = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];
	
      dlogden[p] = (logden_plus-logden_minus)/(2*eps[p]);

    }
    
    for( i=0; i<4; ++i ){
      for( j=0; j<=i; ++j){
        dlogden[i]*dlogden[j];
      {
    }
	}
  */

}

// Gaussian GARCH(1,1) numerical Hessian
void garch_hessian_num(int *status, double *OPG, double *param, double *y, int *T, double *eps){

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
	for( t=0; t<10; ++t ){
		sigma2[0] += y[t]*y[t];
	}
	sigma2[0] /= 10;
  eps[0] = y[0]/sqrt( sigma2[0] );
  
	// loop 
	*loglik = 0;
	for( t=1 ; t<*T ; ++t ){
		sigma2[t] = omega + alpha * y[t-1]*y[t-1] + gamma * y[t-1]*y[t-1]*((double)(y[t-1]<0)) + beta * sigma2[t-1];  
		eps[t]    = y[t]/sqrt( sigma2[t] );

		logden    = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];
		*loglik   += logden;
	}
}

// MEWMA Model Filter
/*
void mewma_filter(int *status, double *Sigma, double* eps, double *loglik, double *param, double *y, int *T, int *N){

	double logden;
	double lambda;
	int t,i,j;
	double ***S, **E, **Y;
 
  	*loglik = 0;

	S = Calloc(T,double**);
	for( t=0; t<T; t++) { 
		S[t] = Calloc(N,double *); 
		for( i=0; i<N; ++i ) {
			S[t][i] = Calloc(N,double); 
			for( i=0; i<N; ++i ){
				S[t][i][i] = 0;
			}
		}
	}

	Y = Calloc(T,double*);
	E = Calloc(T,double*);
	for( i=0; i<T; i++) { 
		Y[i] = Calloc(N,double); 
		E[i] = Calloc(N,double); 
		for( j=0; j<N; j++ ) {
			Y[i][j] = y[ M * j + i ];
			E[i][j] = 0;
		}
	}

	// clean up	
	for(t = 0; t < T; t++){ Free(Y[i]); Free(E[i])}
	Free(Y);
	Free(E);
}

*/
