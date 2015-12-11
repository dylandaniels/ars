Mainline=function(g,x,leftbound=-Inf,rightbound=Inf)
{ y=x
  h=function(y){
    return (g(y))
  }
  h_der=function(h,x){
    y=x
    dh_x=diag(attributes(numericDeriv(quote(h(y)),'y'))$gradien)
return(dh_x)
  }

abscissae=ininitialize()
z=envelopeIntersectionPoints(x, z, abscissae, h)
u=envelope (x, z, abscissae, h)
s=normalizedEnvelope(x, z, abscissae, h)
l=squeezing(h,abscissae,x)
sample=sampleFromEnvelope(abscissae, z, u, h,x)
acceptReject <- function (xStar, l, u, h)
updatestep()
}

# Need to update later.
