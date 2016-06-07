#################################
#
# co-ordination hypothesis (Chen 1993, Maire 2012) to determine photosynthetic parameters
# As implemented in SDGVM
#
# Anthony Walker (walkerap@ornl.gov)
# 21 Jun 2016
#
#################################

library(lattice)

SDGVM_current <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V*c*(4*c + 8*gs) * (1 + (a*Q/(exp(1)*V^0.89))^2)^0.5 - a*Q*c*(c+k)
}

SDGVM_current1 <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V*(4*c + 8*gs) * (1 + (a*Q/(exp(1)*V^0.89))^2)^0.5 - a*Q*(c+k)
}


recent <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V - a*Q / ( (1 + (a*Q/(exp(1)*V^0.89))^2)^0.5 ) * (c+k)/(4*(c+2*gs))
}

recent1 <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V - a*Q / ( (1 + (a*Q/(20+1.63*V))^2)^0.5 ) * (c+k)/(4*(c+2*gs))
}

alt <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V - a*Q / ( (1 + (a*Q/(a*2000))^2)^0.5 ) * (c+k)/(4*(c+2*gs))
}

alt1 <- function(V,a=alpha,Q=q,c=ci,k=Km,gs=gstar,Sv=1,Sj=1) {
  # to be solved with uniroot
  V - a*Q / ( (1 + 1)^0.5 ) * (c+k)/(4*(c+2*gs))
}

recent(0.1)

help(uniroot)

uniroot(recent,interval=c(10,400))
uniroot(SDGVM_current,interval=c(10,100))

alpha <- 0.7
q     <- 100
ci    <- 20
Km    <- 40
gstar <- 4

vcmx <- 5:40
vcmx <- 0.1:10
xyplot(recent(vcmx)~vcmx,abline=0,type='p')
xyplot(alt1(vcmx)~vcmx,abline=0,type='p')

xyplot(recent(vcmx)~vcmx,abline=0,type='l',groups=rep(10:1000,each=length(vcmx)))
xyplot(recent1(vcmx)~vcmx,abline=0,type='l',groups=rep(10:1000,each=length(vcmx)))

xyplot(recent(vcmx)+SDGVM_current(vcmx)~vcmx,abline=0,type='l')

vcmx <- 10:20
xyplot(SDGVM_current(vcmx)+SDGVM_current1(vcmx)~vcmx,abline=0,type='l')

f1 <- function(a,b,c) (a^2*b^2/c^2)
f2 <- function(a,b,c) (a*b/c)^2

a <- rnorm(10)
b <- rnorm(10)
c <- rnorm(10)
xyplot(f1(a,b,c)~f2(a,b,c),abline=list(a=0,b=1))