Mainline <- function(g,abscissae,leftbound=-Inf,rightbound=Inf) { 
  h <- function(y) {
    return (log(g(y)))
  }
  h_der=function(x){
    y=x
    dh_x=diag(attributes(numericDeriv(quote(h(y)),'y'))$gradien)
return(dh_x)
  }

#abscissae=ininitialize() # or get abscissae from user input
hx = h(abscissae)
dhx = h_der(abscissae)
z=envelopeIntersectionPoints(abscissae, hx, dhx)
u=envelope (z, x, xstar, hx, dhx)
s=normalizedEnvelope(x, z, abscissae, h)
l=squeezing(h,abscissae,x)
sample=sampleFromEnvelope(abscissae, z, u, h,x)
acceptReject <- function (xStar, l, u, h)
updatestep()
}

# Need to update later.
