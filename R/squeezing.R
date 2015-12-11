#***********document*********************************************************
#we import h(x),abscissae,and x value And returns l(x)(which is squeezing here),after that if we use exp(l(x)),we get T(x)
#*****************************************************************************
squeezing=function(h,abscissae,x){
  k=length(abscissae)
  a=0
  for(i in 1:(k-1))
  {
    if (abscissae[i]==abscissae[i+1])
    {
      print('abscissae is wrong.')
      a=1
      break
    }
  }
  if (a==0)
  {
  if(x<abscissae[1] || x>abscissae[k])
  {return (-Inf)}

  {for(i in 1:(k-1))
    {
      if(x>=abscissae[i] && x<=abscissae[i+1])
    {
      return ((abscissae[i+1]-x)*h(abscissae[i])+(x-abscissae[i])*h(abscissae[i+1])/(abscissae[i+1]-abscissae[i]))
    }
    }
    }
  }
}




