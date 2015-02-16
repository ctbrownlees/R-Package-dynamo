
.packageName <- "dynamo"

.First.lib <- function(lib, pkg){
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
	dimnames(matcoef) <- list( dimnames(est$param)[[1]] , c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
	obj$matcoef       <- matcoef
  obj$objfun        <- est$obj
  
  # fitted values and residuals
  obj$y     <- job$y
  obj$fit   <- est$fit
  obj$resid <- est$resid
  
  for( i in 1:length(names(est$fit)) ){
    obj[[ names(est$fit)[i] ]] <- est$fit[,i]
  }

	return(obj)
}

`predict.dm` <- function(x,n.ahead=NULL,y.out=NULL) {

  pred <- switch( x$model ,
          garch=garch.predict(x,n.ahead,y.out)
  )

  obj        <- x 
  #obj$type   <- 'static'
  #obj$y.out  <- y.out
  obj$pred   <- pred$pred
  class(obj) <- "predict.dm"
  
  return(obj)  
}

`print.dm` <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  cat("\nCoefficients:\n")
  print( round( x$coef , 3 ) )
	invisible(x)
}

`plot.dm` <- function( x , type=NULL , ... ){
  switch( x$model ,
          garch=arch.plot(x,type),
          tarch=arch.plot(x,type)
  )
  invisible(x) 
}

`summary.dm` <- function( x , ... ){

  obj        <- x
  obj$aic    <- -2*x$objfun/length(x$y) + 2/length(x$y) * length(x$coef)
  obj$bic    <- -2*x$objfun/length(x$y) + log(length(x$y))/length(x$y) * length(x$coef)
  class(obj) <- "summary.dm"
  
  return(obj)
}

`print.summary.dm` <- function( x , digits = max(3, getOption("digits") - 3), ... ){
	cat("\nCall:\n", deparse(x$call), "\n", sep = "")
	cat("\nCoefficients:\n")
	printCoefmat(x$matcoef, digits = digits, signif.stars = TRUE)
	cat("\nLogLik: ",x$objfun,", AIC: ",x$aic,", BIC: ",x$bic,sep='')
  invisible(x)
}

# utilities
`dm.spec` <- function( frml, data ){

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

# TESTS
archlm.test <- function(){ }

dm.test <- function(){ }