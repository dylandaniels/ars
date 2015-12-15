#***********document*********************************************************
#we import h(x),abscissae,and x value And returns l(x)(which is squeezing here),after that if we use exp(l(x)),we get T(x)
#*****************************************************************************
squeezing=function(hx,abscissae,x){
  k=length(abscissae)
  if(x<abscissae[1] || x>abscissae[k])
  {return (-Inf)}
  for(i in 1:(k-1))
  {
    if(x>=abscissae[i] && x<=abscissae[i+1])
    {
      return (((abscissae[i+1]-x)*hx[i]+(x-abscissae[i])*hx[i+1])/(abscissae[i+1]-abscissae[i]))
    }
  }
}




