# This is equation (2) on page 340
# x is the point to evaluate
# z is the vector of envelop intersection points
# sampleFunc is h(x)
# derivSampleFunc is h'(x)


#I'm not sure what the output should be at this moment
envelope <- function(x,z,sampleFunc){
  
  x <- as.numeric(sort(x))
  derivX <- numericDeriv(quote(sampleFunc),'x')
  dX <- diag(attr(derivX,'gradient'))
  dimZ <- length(z)
  z <- sort(z)
  u <- vector('list',dimZ)
  
  for (i in 0:(dimZ-1)){
    u[[i+1]] <- function(x)sampleFunc(x[i+1])+(quote(x)-x[i+1])*dX[i+1]
  }
  return(u)
}


##need to modify later
updateStep <- function(z,xstar){
  z <- c(z,xstar)
  z <- sort(z)
  #call u_k+1(x)
  #call s_k+1(x)
  #call l_k+1(x)
  
}