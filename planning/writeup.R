devtools::load_all('.')
library(optimx)
ddgamma=Vectorize(function(x,alpha=4,beta=5)
{
  return(dgamma(x,shape=alpha,rate=beta))
})
ddbeta=Vectorize(function(x,alpha=2,beta=2)
{
  return(dbeta(x,shape1=alpha,shape2=beta))
})
ddchisqu=Vectorize(function(x,k=6)
{
  return(1/(2^(k/2)*gamma(k/2))*x^(k/2-1)*exp(-x/2))
})

ppgamma=Vectorize(function(x,alpha=4,beta=5)
{
  return(pgamma(x,shape=alpha,rate=beta))
}
ppbeta=Vectorize(function(x,alpha=2,beta=2)
{
  return(pbeta(x,shape1=alpha,shape2=beta))
})

ddchisq=Vectorize(function(x,df=6)
{
  return(dchisq(x,df))
}
)
ppchisq=Vectorize(function(x,df=6)
{
  return(pchisq(x,df))
}
)




set.seed(0)
samples=ars(n=1000, g=dnorm, dg=NULL, initialPoints=NULL, leftbound=-Inf, rightbound=Inf)
png(filename="the_distribution_of_samples_vs_normal_distribution.png")
hist(samples,probability=TRUE,main=NULL)
plot(dnorm,xlim=c(-3,3),add=TRUE,col="blue",main='the distribution of samples vs normal distribution')
title(main='the distribution of samples vs normal(0,1) distribution')
dev.off()
ks.test(samples,pnorm)$p.value
#From the figure, we can see that the samples approximately follow normal(0,1) distribution.
#And because the p value is 0.1382014,we do not reject that he samples follow normal(0,1) distribution,
#so it means that method is OK.

set.seed(0)
samples=ars(1000, g=dunif,leftbound=0, rightbound=1)
png(filename="the_distribution_of_samples_vs_uniform(0,1)_distribution.png")
hist(samples,probability=TRUE,main=NULL)
plot(dunif,xlim=c(0,1),add=TRUE,col="blue",main='the distribution of samples vs uniform distribution')
title(main='the distribution of samples vs uniform(0,1) distribution')
dev.off()
ks.test(samples,punif)$p.value
#From the figure, we can see that the samples approximately follow uniform(0,1) distribution.
#And because the p value is 0.2363632,we do not reject that he samples follow uniform(0,1) distribution,
#so it means that method is OK.

set.seed(0)
samples=ars(1000, g=ddgamma,leftbound=0, rightbound=Inf)
png(filename="the_distribution_of_samples_vs_gamma(4,5)_distribution.png")
hist(samples,probability=TRUE,main=NULL)
plot(ddgamma,xlim=c(0,10),add=TRUE,col="blue",main='the distribution of samples vs gamma(4,5) distribution')
title(main='the distribution of samples vs gamma(4,5) distribution')
dev.off()
ks.test(samples,ppgamma)$p.value
#From the figure, we can see that the samples approximately follow gammma(4,5) distribution.
#And because the p value is 0.6165001,we do not reject that he samples follow gamma(4,5) distribution,
#so it means that method is OK.


set.seed(0)
samples=ars(1000, g=ddbeta,leftbound=0, rightbound=1)
png(filename="the_distribution_of_samples_vs_beta(2,2)_distribution.png")
hist(samples,probability=TRUE,main=NULL)
plot(ddbeta,xlim=c(0,10),add=TRUE,col="blue",main='the distribution of samples vs beta(2,2) distribution')
title(main='the distribution of samples vs beta(2,2) distribution')
dev.off()
ks.test(samples,ppbeta)$p.value
#From the figure, we can see that the samples approximately follow beta(2,2) distribution.
#And because the p value is  0.1823897,we do not reject that he samples follow beta(2,2) distribution,
#so it means that method is OK.


set.seed(0)
samples=ars(1000, g=ddchisq,leftbound=0, rightbound=Inf)
png(filename="the_distribution_of_samples_vs_chisq(6)_distribution.png")
hist(samples,probability=TRUE,main=NULL)
plot(ddchisqu,xlim=c(0,20),add=TRUE,col="blue",main='the distribution of samples vs chisq(6) distribution')
title(main='the distribution of samples vs chisq(6) distribution')
dev.off()
ks.test(samples,ppchisq)$p.value
#From the figure, we can see that the samples approximately follow chisq(6) distribution.
#And because the p value is  0.149094,we do not reject that he samples follow chisq(6) distribution,
#so it means that method is OK.
