
# Gaussian GARCH(1,1) functions
garch.filter <- function( y , param ){

  T      <- length(y)
  
  if( any(!is.finite(param)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) , eps=rep(NA,T) )    
    return( filter ) 
  }
  
  filter <- .C('garch_filter', 
               status = as.integer(0), 
               sigma2 = as.double(rep(0,T)), 
               eps    = as.double(rep(0,T)), 
               loglik = as.double(0), 
               as.double(param),
               as.double(y),
               as.integer(T),
               PACKAGE="dynamo")
  
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
    res <- nlminb( param.init, obj, lower=c(0,1e-5,0), upper=c(1,1,1) )
    param.est <- res$par
  }
  else {
    param.est <- param.init 
  }
  
  filter <- garch.filter(y,param.est)
  vcv    <- vcv.mle( param.est , obj , 0.0001 * abs(param.est) )
  sigma2 <- data.frame( sigma2=filter$sigma2 )
  eps    <- data.frame( eps=filter$eps )
  
  param.est           <- as.array(param.est)
  dimnames(param.est) <- list( c('(Intercept)','ARCH','GARCH') )
  
  list( param=param.est , 
        fit=sigma2, 
        resid=eps,
        vcv=vcv ,
        obj=filter$loglik )
}

garch.predict <- function( x , n.ahead=NULL, y.out=NULL ){
  
  if( !is.null(n.ahead) ){ type <- 'dynamic' }
  if( !is.null(y.out)   ){ type <- 'static' }
  
  if( type=='static' ){
    y      <- c( x$y , y.out , 0 )
    filter <- garch.filter( y , x$coef )
    pred   <- filter$sigma2[length(x$y):length(y)]
  }
  if( type=='dynamic' ){
    #y      <- c( x$y ,  )
    #filter <- garch.filter( y , x$coef )
    #pred   <- filter$sigma2[length(x$y):length(y)]
  }
  
  list( type=type , pred=pred )
}

arch.plot <- function( x , type ){
  
  if( is.null(type) ){ type <- 'vol' }
  
  if( type=='vol' ){

    plot( sqrt(x$sigma2) , t='l' , col='red2' , lwd=2 , tck=0.02 , xaxs='i' , xlab='time' , ylab='conditional volatility' ) 
    grid()

    return()
  } else if( type=='series' ){

    plot( x$y , t='p' , col='orange2' , lwd=2 , tck=0.02 , xaxs='i' , xlab='time' , ylab='') 
    lines( 2*sqrt(x$sigma2) , col='red2') 
    lines( -2*sqrt(x$sigma2) , col='red2' ) 
    grid()

    return()
  }
  
}

# Gaussian TARCH(1,1) functions
tarch.filter <- function( y , param ){
  
  T      <- length(y)

  if( any(!is.finite(param)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'tarch_filter', 
                status = as.integer(0), 
                sigma2 = as.double(rep(0,T)) , 
                eps    = as.double(rep(0,T)) , 
                loglik = as.double(0) , 
                as.double(param) , 
                as.double(y) , 
                as.integer(T) , 
                PACKAGE="dynamo" )
  
  filter = list( loglik=result$loglik , sigma2=result$sigma2 )
  
  return(filter)
}

tarch.fit <- function(y){
  
  obj  <- function(x){ return( -tarch.filter(y,x)$loglik ) }  
  
  list( a=0 )
}

# Gaussian APARCH filter
aparch.filter <- function( x , params ){
  
  T   <- length(x)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'aparch_filter', 
                status = as.integer(0), 
                sigma2 = as.double(rep(0,T)) , 
                eps    = as.double(rep(0,T)) , 
                loglik = as.double(0) , 
                as.double(params) , 
                as.double(x) , 
                as.integer(T)
                )
  
  return(list( loglik=result$loglik , sigma2=result$sigma2 ))
}

aparch.fit <- function(x){
  
  # Initialise parameters and set bounds
  Tx <<- x
  Meanx = mean(Tx); Varx = var(Tx); S = 1e-3
  
  params_init = c(mu = Meanx, omega = 0.1*Varx, alpha = 0.1, gam1= 0.02, beta = 0.81,delta=2)
  lowerBounds = c(mu = -10*abs(Meanx), omega = S, alpha = S, gam1= -(1-S), beta = S,delta=0.1)
  upperBounds = c(mu = 10*abs(Meanx), omega = 10*Varx, alpha = 1-S, gam1 = (1-S), beta = 1-S,delta=4)
  
  # Optimise -log-likelihood and calculate Hessian matrix
  
  fit = nlminb(start = params_init, objective = llh,
               lower = lowerBounds, upper = upperBounds) # , control = list(trace=3))

  hess <- Hessian(fit$par) 
  
  cat("Log likelihood at MLEs: ","\n")
  print(-llh(fit$par))
  
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(abs(diag(solve(hess))))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  # compute output
  est=fit$par
  mu = est[1]; omega = est[2]; alpha = est[3]; gam1=est[4]; beta = est[5]; delta = est[6]
  z= Tx-mu
  sigma.t = aparch.filter(x,est)$sigma2
  
  return(list(summary = matcoef, residuals = z, volatility = sigma.t, par=est, n.loglik = -fit$obj))

}

# Gaussian sv(1) functions
sv.filter <- function( y , u , z , param ){
  
  T      <- length(y)
  P      <- ncol(u)
  
  if( any(!is.finite(param)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) , eps=rep(NA,T) )    
    return( filter ) 
  }
  
  filter <- .C('sv_filter', 
               status    = as.integer(0), 
               sigma2.pr = as.double(rep(0,T*P)), 
               sigma2.up = as.double(rep(0,T*P)), 
               loglik    = as.double(0), 
               as.double(param),
               as.double(y),
               as.double(z),
               as.double(u),
               as.integer(T),
               as.integer(P),
               PACKAGE="dynamo")
  
  filter = list( loglik=filter$loglik , sigma2=filter$sigma2 , eps=filter$eps )
  
  return(filter)
}
