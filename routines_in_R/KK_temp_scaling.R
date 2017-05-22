#################################
#
# Kattge & Knorr 2007 temperature scaling and acclimation of photosynthetic parameters
# As implemented in SDGVM
#
# Anthony Walker (walkerap@ornl.gov)
# 21 Jan 2016
#
#################################

# functions
# use fixed delta S based on optimum temperature and other parameters
f_topt_dS <- function(topt,pars,R=8.31446) {
  # no statement to account for Jmax difference in topt
  # - Eq 2 K&K 2007
  
  pars$Hd/(topt+273.15) + ( R*log(pars$Ha/(pars$Hd-pars$Ha)) )
}

# employ Kattge&Knorr scaling of delta S based on mean growth temp 
f_tacc_dS <- function(tgrowth,jv='v') {
  # - Table 3 K&K 2007
  if(jv=='v') 668.39 - 1.07*tgrowth else 659.70 - 0.75*tgrowth  
}

# see Medlyn etal 2002 or Kattge&Knorr 2007  
f_t_scalar <- function(Tsk,Trk,pars,R=8.31446) {
  # This function is written to give the temperature scalar on a parameter at a reference temperature (Trk) 
  # to an environment temperature (Tsk) 
  
  exp(pars$Ha*(Tsk-Trk) / (R*Tsk*Trk)) * ( 
    (1 + exp((Trk*pars$dS-pars$Hd) / (Trk*R)) ) 
    /(1 + exp((Tsk*pars$dS-pars$Hd) / (Tsk*R)) ) )    
}

f_JV_ratio <- function(tgrowth,Tr) {
  # - Table 3 K&K 2007
  # - this function returns the ratio of the JV ratio at the acclimated temp to the JV ratio at the reference temp
  # - the reference temperture (usuallly 25oC) is assumed to be both the measurement temp and the acclimated temp of the reference J and V values   
  
  ((2.59-0.035*tgrowth)/(2.59-0.035*Tr))
}

# see K&K 2007
f_kkt_scalar <- function(Ts,Tr=25,pars,topt=25,tacc=F,tgrowth=25) {
  
  # calculate temps in kelvin
  Tsk <- Ts + 273.15
  Trk <- Tr + 273.15
  
  # calculate deltaS
  # K&K 2007 Table 3 - only deltaS and the Jmax:Vcmax ratio significantly acclimate to growth temperature
  deltaS <- if(tacc) f_tacc_dS(tgrowth,pars$jv) else f_topt_dS(topt,pars)
  
  # add to pars
  pars <- c(pars,dS=deltaS)
#   for(i in 2:4) pars[[i]] <- pars[[i]]*1e-3
#   print(pars)
  
  # calculate temperature scalar
  t_scalar <- f_t_scalar(Tsk,Trk,pars) 
  
  # adjust jmax:vcmax ratio based on growth temp
  if(pars$jv=='j'&tacc) t_scalar * f_JV_ratio(tgrowth,Tr)  else t_scalar
  
}


t <- 25
t <- 1:50
tk <- t + 273.15
# below are the kinetic parameters from Farquhar etal 1980
kc  = exp(35.8 - 80.5/(0.00831*tk))
ko  = exp(9.6 - 14.51/(0.00831*tk))*1000.0
tau = exp(-3.949 + 28.99/(0.00831*tk))

# the above are now deprecated for the in vivo parameters from Bernacchi etal 2001
kc  = 40.49 * exp((79430.0 / (8.31 * 298.15)) * (1.0 - (298.15) / (273.15 + t)))
ko  = 27840. * exp((36380.0 / (8.31 * 298.15)) * (1.0 - (298.15) / (273.15 + t)))
c_p = 4.275 * exp((37830.0 / (8.31 * 298.15)) * (1.0 - (298.15) / (273.15 + t)))
tau = 0.5*21000/c_p






# non temperature dependent parameters
# - Table 3 K&K 2007
p_vcmax <- list(jv='v',Ha=71513,Hd=200000) 
p_jmax  <- list(jv='j',Ha=49884,Hd=200000) 
         


# plots
library(lattice)
ts <- 1:50
p1 <- 
xyplot(f_kkt_scalar(ts,pars=p_vcmax)+
         f_kkt_scalar(ts,pars=p_vcmax,tacc=T)+
         f_kkt_scalar(ts,pars=p_vcmax,tacc=T,tgrowth=30)+
         f_kkt_scalar(ts,pars=p_jmax)+
         f_kkt_scalar(ts,pars=p_jmax,tacc=T)+
         f_kkt_scalar(ts,pars=p_jmax,tacc=T,tgrowth=30)~
         ts,
       type='l',lty=1:3,col=rep(c('blue','red'),each=3),ylim=c(0,3),
       xlab=expression('Environment Temperature  ['*degree*C*']'),
       ylab=expression('Scalar [from 25'*degree*C*']'),
       abline=list(h=0:1,v=c(25)),
       key=list(corner=c(0,1),x=0,y=0.95,border=T,
                lines=list(lty=1:3,col=rep(c('blue','red'),each=3)),
                text=list(c(expression(V[cmax]*' Topt 25'*degree*c),
                          expression(V[cmax]*' Tacc 25'*degree*c),
                          expression(V[cmax]*' Tacc 30'*degree*c),
                          expression(J[max]*'   Topt 25'*degree*c),
                          expression(J[max]*'   Tacc 25'*degree*c),
                          expression(J[max]*'   Tacc 30'*degree*c)
                          ))) )

p2 <-
xyplot(f_kkt_scalar(ts,pars=p_vcmax)+
         f_kkt_scalar(ts,pars=p_vcmax,tacc=T,tgrowth=5)+
         f_kkt_scalar(ts,pars=p_vcmax,tacc=T,tgrowth=10)+
         f_kkt_scalar(ts,pars=p_jmax)+
         f_kkt_scalar(ts,pars=p_jmax,tacc=T,tgrowth=5)+
         f_kkt_scalar(ts,pars=p_jmax,tacc=T,tgrowth=10)~
         ts,
       type='l',lty=1:3,col=rep(c('blue','red'),each=3),ylim=c(0,2),
       xlab=expression('Environment Temperature  ['*degree*C*']'),
       ylab=expression('Scalar [from 25'*degree*C*']'),
       abline=list(h=0:1,v=c(25)),
       key=list(corner=c(0,1),x=0,y=0.95,border=T,
                lines=list(lty=1:3,col=rep(c('blue','red'),each=3)),
                text=list(c(expression(V[cmax]*' Topt 25'*degree*c),
                          expression(V[cmax]*' Tacc 5'*degree*c),
                          expression(V[cmax]*' Tacc 10'*degree*c),
                          expression(J[max]*'   Topt 25'*degree*c),
                          expression(J[max]*'   Tacc 5'*degree*c),
                          expression(J[max]*'   Tacc 10'*degree*c)
                          ))) )

pdf('K&K 2007 Temp acclimation.pdf',width=8,height=8)
print(p1,split=c(1,1,1,2),more=T)
print(p2,split=c(1,2,1,2),more=F)
dev.off()



