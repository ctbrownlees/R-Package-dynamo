
# Bivariate DCC
bidcc.filter <- function(y,param){
  
  T      <- nrow(y)
  
  result <- .C( 'bidcc_filter', 
                status = as.integer(0), 
                rho    = as.double(rep(0,T)), 
                eps    = as.double(rep(0,T*2)),
                loglik = as.double(0),
                as.double(param),
                as.double(y),
                as.integer(T), 
                PACKAGE="dynamo")

  filter = list( rho=result$rho , loglik=result$loglik )
  
  return(filter)
}

bidcc.fit <- function(y,opts){
  
  # INPUT 
  if( is.null(opts$param) ){  
    param.init <- c( 0.05 , 0.90 ) 
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
  obj  <- function(x){ return( -bidcc.filter(y,x)$loglik ) }  
  
  if( fit==TRUE ){ 
    res <- nlminb( param.init, obj, lower=c(1e-5,0), upper=c(1,1) )
    param.est <- res$par
  }
  else {
    param.est <- param.init 
  }
  
  filter <- bidcc.filter(y,param.est)
  vcv    <- vcv.mle( param.est , obj , 0.0001 * param.est )
  #vcv    <- matrix(c(1,0,0,1),2,2)
  rho    <- data.frame( rho=filter$rho )
  eps    <- data.frame( eps=filter$eps )
  
  param.est           <- as.array(param.est)
  dimnames(param.est) <- list( c('alpha','beta') )
  
  print( param.est) 
  
  list( param=param.est , 
        fit=rho, 
        resid=eps,
        vcv=vcv ,
        obj=filter$loglik )
}
