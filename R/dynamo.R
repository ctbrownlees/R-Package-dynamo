
.packageName <- "dynamo"

.First.lib <- function(lib, pkg)
{
     library.dynam("dynamo", pkg, lib)
}

# dynamo interface 
`dm` <- function(formula,data=parent.frame(),opts=NULL){

	call <- match.call()

	job <- dm.spec(formula,data)

	est <- switch( job$mdl ,
	  garch=garch.fit(job$y,opts),
	  tarch=tarch.fit(job$y,opts),
	  mewma=mewma.fit(job$y,opts)
	)

	obj <- list( call=call )
	class(obj) <- 'dm'

	# model info
	obj$model  <- job$mdl

	# coef
	obj$coef          <- est$param
	obj$vcv           <- est$vcv
	se.coef           <- sqrt(diag(obj$vcv))
	tval              <- obj$coef/se.coef
	matcoef           <- cbind( obj$coef , se.coef, tval, 2*(1-pnorm(abs(tval))))
	dimnames(matcoef) <- list( est$param.names , c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
	obj$matcoef       <- matcoef

	return(obj)
}

`predict.dm` <- function(x) {

}

`print.dm` <- function (x, digits = max(3, getOption("digits") - 3), ...) {
	cat("\nCall:\n", deparse(x$call) , sep = "")
	cat("\n\nCoefficient(s):\n")
	printCoefmat(x$matcoef, digits = digits, signif.stars = TRUE)
	invisible(x)
}

`plot.dm` <- function( x , ... ){
}

`summary.dm` <- function( x , ... ){
	invisible(x)
}

print.summary.dm <- function( x , ... ){
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

	#cat("\nCoefficient(s):\n")
	#printCoefmat(x$matcoef, digits = 6, signif.stars = TRUE)
}

# Utilities
dm.spec <- function( frml, data){

    # variable initialization
    frml.char <- as.character(frml) 
    y   <- data[[ frml.char[2] ]]
    mdl <- frml.char[3]

    job <- list( mdl=mdl , y=y )
}

vcv.mle <- function( x , obj , epsilon ){
	npar    <- length(x)
	Hessian <- matrix(0,npar,npar)
	for (i in 1:npar) {
		for (j in 1:npar) {
			x1 = x2 = x3 = x4 = x
			x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
			x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
			x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
			x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
			Hessian[i, j] = (obj(x1)-obj(x2)-obj(x3)+obj(x4))/ (4*epsilon[i]*epsilon[j])
		}
	}
	vcv <- solve(Hessian)
}

# Gaussian GARCH(1,1) functions
garch.filter <- function( y , param ){

	T      <- length(y)

	filter <- .C('garch_filter', 
		status = as.integer(0), 
		sigma2 = as.double(rep(0,T)) , 
		eps    = as.double(rep(0,T)) , 
		loglik = as.double(0) , 
		as.double(param) , 
		as.double(y) , 
		as.integer(T) , 
		PACKAGE="dynamo" )

        filter = list( loglik=filter$loglik , sigma2=filter$sigma2 , eps=filter$eps )
  
	return(filter)
}

garch.fit <- function(y,opts){

	# INPUT 
	if( is.null(opts$param) ){  
		param.init <- c( var(y)*(0.05) , 0.05 , 0.90 ) 
	}
	else {
		param.init <- opts$param.init 
	}
	if( is.null(opts$fit) ){ 
		fit <- TRUE
	}
	else { 
		fit <- as.logical( opts$fit ) 
	}

	# MAIN
	obj  <- function(x){ return( -garch.filter(y,x)$loglik ) }  

	if( fit==TRUE ){ 
		res <- nlminb( param.init, obj, lower=c(0,0,0), upper=c(1,1,1) )
		param.est <- res$par
	}
	else {
		param.est <- param.init 
	}

	filter <- garch.filter(y,param.est)
	vcv    <- vcv.mle( param.est , obj , 0.0001 * param.est )

	list(   param=param.est   , param.names=c('intercept','ARCH','GARCH') , 
		fit=filter$sigma2 , fit.names='sigma2', 
		resid=filter$eps  , res.names='res',
		vcv=vcv ,
                loglik=filter$loglik )
}

garch.predict <- function( x ){

}

# Gaussian TARCH(1,1) functions
tarch.filter <- function( y , param ){
  
  T      <- length(y)
  
  result <- .C( 'tarch_filter', 
                status = as.integer(0), 
                sigma2 = as.double(rep(0,T)) , 
                eps    = as.double(rep(0,T)) , 
                loglik = as.double(0) , 
                as.double(param) , 
                as.double(y) , 
                as.integer(T) , 
                PACKAGE="dynamo" )
  
  filter = list( sigma2=result$sigma2 , loglik=result$loglik )
  
  return(filter)
}

tarch.fit <- function(y){
  
  obj  <- function(x){ return( -tarch.filter(y,x)$loglik ) }  
  der  <- function(x){ return( nl.grad(x,obj) ) }
  
  # initial values
  x0 <- c( var(y)*(0.05) , 0.05 , 0.0 , 0.90 )
  
  opts <- list("algorithm"="NLOPT_LD_LBFGS",
               "xtol_rel"=1.0e-8)
  
  #
  res <- nloptr( x0=x0, 
                 eval_f     =obj,
                 eval_grad_f=der, 
                 opts=opts)
  
  print( res )
  
  list( param=res$solution )
}

# Multivariate EWMA
mewma.filter <- function( y , param ){

	T      <- nrow(y)
	N      <- ncol(y)

	result <- .C( 'mewma_filter', 
		status = as.integer(0), 
		Sigma  = as.double(rep(0,T*N)) , 
		eps    = as.double(rep(0,T*N)) , 
		loglik = as.double(0) , 
		as.double(param) , 
		as.double(y) , 
		as.integer(T) , 
		as.integer(N) , 
		PACKAGE="dynamo" )

  	filter = list( sigma2=result$sigma2 , loglik=result$loglik )
  
	return(filter)
}

mewma.fit <- function(y,param,fit){
  
  obj  <- function(x){ return( -mewma.filter(y,x)$loglik ) }  
  der  <- function(x){ return( nl.grad(x,obj) ) }
  
  # initial values
  if( fit==TRUE ){
	  x0 <- c( var(y)*(0.05) , 0.05 , 0.90 )
  
	  opts <- list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8)
  
	  res <- nloptr( x0=x0, 
                 eval_f     =obj,
                 eval_grad_f=der, 
                 opts=opts)
  
	  print( res )
  }
  else {
	out <- mewma.filter(y,param)
  }  

  #C <- solve( hessian( obj , res$solution ) )
  
  list( param=res$solution , param.names <-c('intercept','ARCH','GARCH') , vcv=C)
}

