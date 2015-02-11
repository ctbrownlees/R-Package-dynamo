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
