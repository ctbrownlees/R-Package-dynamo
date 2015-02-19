
#include <R.h>
#include <math.h>

// utilities
double **create_real_matrix(int rows, int cols){
  int i;
  double **mat = (double **) Calloc(rows,double*);
  for(i = 0; i < rows; i++) mat[i] = (double *) Calloc(cols,double);
  return mat;
}

double **create_and_copy_real_matrix(int rows, int cols, double *data){
  int i,j;
  double **mat = (double **) Calloc(rows,double*);
  for(i = 0; i < rows; i++){
    mat[i] = (double *) Calloc(cols,double);
    for(j=0;j<cols;++j) mat[i][j] = data[j * rows + i];
  }
  return mat;
}

void destroy_real_matrix(double **matrix, int rows, int cols){
  int i;
  for(i = 0; i < rows; i++){ Free(matrix[i]); }
	Free(matrix);
}

// Gaussian GARCH(1,1) Filter
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
  
	// loop 
	*loglik = 0;
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
	for( t=0; t<10; ++t ){ sigma2[0] += y[t]*y[t]; }
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

// BIDCC Model Filter
void bidcc_filter(int *status, double *rho, double* eps, double *loglik, double *param, double *_y, int *T){

  int t;
	double logden;
	double alpha,beta;
	double **Q, **y;
  double rho_bar;
  *loglik = 0;
  
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
  	*loglik += logden;    
  }
  
  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }
  
  // cleanup
  destroy_real_matrix(Q,*T,3);
  destroy_real_matrix(y,*T,2);
}
