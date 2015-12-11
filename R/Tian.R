# This is equation (2) on page 340
# x is the point to evaluate
# z is the vector of envelop intersection points
# sampleFunc is h(x)
# derivSampleFunc is h'(x)

#take the input vector of z(sorted), one value of xstar, vector of h(x) and vector of h'(x)
#and output the value of u_k(xstar)
envelope <- function(z,xstar,sampleFunc,derivX){
  #calculate the vector of absolute difference between z and xstar
  diffs <- abs(xstar-z)
  #calcuate the index of the min of the differences
  index <- which.min(diffs)
  
  #for the z_index which has the min absolute difference with xstar,
  #either it's the upper bound or lower bound of the interval in which
  #xstar falls
  if(xstar>=z[index-1] && xstar<=z[index]){
    ustar <- sampleFunc[index]+(xstar-z[index])*derivX[index]
  }else if (xstar>z[index] && xstar<=z[index+1]){
    ustar <- sampleFunc[index+1]+(xstar-z[index+1])*derivX[index+1]
  }
  
  #return ustar
  return(ustar)
}

#to test the function
z <- rnorm(10)
z <- sort(z)
minZ <- min(z)
maxZ <- max(z)
xstar <- runif(1,minZ,maxZ)
sampleFunc <- rnorm(10)
derivX <- rnorm(10)

envelope(z,xstar,sampleFunc,derivX)

##
#take the vector of z, current value of xstar, and the result of rejection-or-not
updateStep <- function(z,xstar,result,x,hx,dx){
  #if the result is true, add xstar to T_x, and update z
  if (result=TRUE){
    #add xstar to T_x and sort
    x <- c(x,xstar)
    x<- sort(x)
    
    #update z
    z <- envelopeIntersectPoints(x,distFun,derivFun = NULL)
    
    ##update hx and dx,maybe we can use envelopeIntersectPoints fucntion
    ##to also output hx and dx?
    hx <- hx(x)
    dx <- dx(x)
  }
  
  
}



