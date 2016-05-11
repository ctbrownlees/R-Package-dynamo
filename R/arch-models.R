
# Gaussian GARCH(1,1) functions
garch.filter <- function( y , param ){
  
  T      <- length(y)
  
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

# garch.opg <- function( y , param ){
#   
#   T      <- length(y)
#   
#   filter <- .C('garch_opg_num', 
#                status = as.integer(0), 
#                OPG     = as.double(rep(0,4*4)), 
#                as.double(param), 
#                as.double(y), 
#                as.integer(T), 
#                PACKAGE="dynamo" )
#   
#   return(OPG)
# }

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
  
  list( a=0 )
}
