###############################
#
# Respiration routines in SDGVM
#
# AWalker May 2019
#
###############################


library(lattice)

t <- seq(-10,40,0.1)

sapresp <- function(t) {
  si <- 0.4 * exp(0.02*t) * 0.35
  si * 3600 * 24 * 1e-6
}

xyplot(sapresp(t)~t)

xyplot(I(0.17*sapresp(t))~t)


