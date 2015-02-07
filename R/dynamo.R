
.packageName <- "dynamo"

.First.lib <- function(lib, pkg)
{
     library.dynam("dynamo", pkg, lib)
}

`dm.spec` <- function( frml, data){

    # variable initialization
    frml.char <- as.character(frml) 
    y   <- data[[ frml.char[2] ]]
    mdl <- frml.char[3]

    job <- list( mdl=mdl , y=y )
}

`dm` <-
function(formula,data=parent.frame(),opts=NULL){

  call <- match.call()
  
  job <- dm.spec(formula,data)
  
  est <- switch( job$mdl ,
          garch=garch.fit(job$y,opts),
          tarch=tarch.fit(job$y,opts),
          mewma=mewma.fit(job$y,opts)
  )
  
  obj <- list( call=call )
  class(obj) <- 'dm'

  obj$model  <- job$mdl
  obj$coef   <- est$param
  obj$vcv    <- est$vcv

  obj
}

print.dm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n")
    invisible(x)
}

plot.dm <- function( x , ... ){
}

summary.dm <- function( x , ... ){
}

print.summary.dm <- function( x , ... ){
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
}

garch.filter <- function( y , param ){

	T      <- length(y)

	result <- .C( 'garch_filter', 
		status = as.integer(0), 
		sigma2 = as.double(rep(0,T)) , 
		eps    = as.double(rep(0,T)) , 
		loglik = as.double(0) , 
		as.double(param) , 
		as.double(y) , 
		as.integer(T) , 
		PACKAGE="dynamo" )

        filter = list( loglik=result$loglik , sigma2=result$sigma2 , eps=results$eps )
  
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
	if( fit==TRUE ){ 

		obj  <- function(x){ return( -garch.filter(y,x)$loglik ) }  

		res <- nlminb( param.init, obj, lower=c(0,0,0), upper=c(1,1,1) )
		param.est <- res$par

		print( param.est )
	}
	else {
		param.est <- param.init 
	}

	filter <- garch.filter(y,param.est)

	#C <- solve( hessian( obj , res$solution ) )

	list( param=param.est , param.names=c('intercept','ARCH','GARCH') , 
		fitted=filter$sigma2 , fitted.names='sigma2', 
		residuals=1 )
}

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

