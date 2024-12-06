# R script to generate all figures and analyses for 
# "Overviewing the exponential model of demand and introducing a simplification that solves issues of span, scale, and zeros" 

# Libararies ----
library(tidyverse)
library(dplyr)
library(tidyr)
library(beezdemand)
library(data.table)
library(nls.multstart)
library(nlme)
library(emmeans)
library(psych)
library(weights)
library(VGAM)

# Handy Functions ----
tickFunc <- function(from = .0001, to = .001, jumps = 6){
    tickHolder <- vector("numeric")
    for(a in 1:jumps){
        ticks <- seq(from*(10^a), to*(10^a), by = to*(10^a)/10)
        tickHolder<- append(tickHolder, ticks)
    }
    return(tickHolder)
}

EXPDToSND <- function(k, base = 10){
    - 1/(log(1 - (1/(k*log(base)))))
}
# HS Normalized Decay ----
svg("HSNormalizedDecay.svg", 7.35, 7.35)
par(mfrow = c(2,2), pty = "s", mar = c(4, 3, 3, .5))

q0 <- 25
k <- 3
alpha <- .04
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy",
      main = expression("Non-Normalized\n     Consumption"~":"~italic(Q)[0] %*%~ 10 ^ {italic(k)  (italic(e) ^{-ring(alpha)*italic(C)}-1)}),
      xlab = NA, ylab = NA,
      from = .0001, to = 10000, cex.main = 1.35, axes = FALSE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(from = .0002, .001), tcl = -.25, labels = NA)
axis(2, at = tickFunc(from = .0002, .001), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

q0 <- 500
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*x)- 1)), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)

q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*x)- 1)), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lwd = 2)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500, italic(P)[max]==4.3)),
                                expression(list(italic(Q)[0]==25, italic(P)[max]==4.3)),
                                expression(list(italic(Q)[0]==1, italic(P)[max]==4.3))
), bty = "n",  lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)


q0 <- 25
k <- 3
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression("   Normalized\nConsumption"~":"~italic(Q)[0] %*%~ 10 ^ {italic(k)  (italic(e) ^{-alpha*italic(Q)[0]*italic(C)}-1)}),
      ylab = NA, xlab = NA,
      from = .0001, to = 10000, cex.main = 1.35, axes = FALSE, lty = 5, lwd = 1.25, line = 1.25)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)

q0 <- 500
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)


q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lwd = 2)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500, italic(P)[max]==.01)),
                                expression(list(italic(Q)[0]==25, italic(P)[max]==.17)),
                                expression(list(italic(Q)[0]==1, italic(P)[max]==4.3))
), bty =  "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)


q0 <- 25
k <- 3
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy",
      main = expression("Non-Normalized\n       Expenditure"~":"~italic(C) %*%~ italic(Q)[0] %*%~ 10 ^ {italic(k)  (italic(e) ^{-ring(alpha)*italic(C)}-1)}),
      ylab = NA, xlab = NA, from = .0001, to = 10000, 
      cex.main = 1.35, axes = FALSE, frame = TRUE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)

points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), pch = 4)
Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1)))
abline(h = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), lty = 3)

q0 <- 500
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*x)- 1))), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), pch = 4)
Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1)))
abline(h = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), lty = 3)

q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*x)- 1))), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lwd = 2)
points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), pch = 4)
Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1)))
abline(h = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), lty = 3)

legend("bottomright", legend = c(expression(list(italic(Q)[0]==500, italic(O)[max]==721.73)),
                                expression(list(italic(Q)[0]==25, italic(O)[max]==36.09)),
                                expression(list(italic(Q)[0]==1, italic(O)[max]==1.44))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


q0 <- 25
k <- 3
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression(" Normalized\nExpenditure"~":"~italic(C) %*%~ italic(Q)[0] %*%~ 10 ^ {italic(k)  (italic(e) ^{-alpha*italic(Q)[0]*italic(C)}-1)}),
      ylab = NA, xlab = NA,
      from = .0001, to = 10000, cex.main = 1.35, axes = FALSE, frame = TRUE, lty = 5, lwd = 1.5, line = 1.25)
points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 4)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)

q0 <- 500
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 4)

q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(q0*alpha)
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lwd = 2)
points(x = Pmax, y = Pmax * (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 4)
abline(h = Pmax * (q0 * 10^(k*(exp(-alpha*Pmax)- 1))), lty = 3)
legend("bottomright", legend = c(expression(list(italic(Q)[0]==500, italic(O)[max]==1.44)),
                                expression(list(italic(Q)[0]==25, italic(O)[max]==1.44)),
                                expression(list(italic(Q)[0]==1, italic(O)[max]==1.44))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

dev.off()

# EXPD standardized ----
svg("EXPDStandardized.svg", 7.35, 7.35)
par(mfrow = c(2,2), pty = "s", mar = c(4, 3, 3, .5))
alpha <- .01
k <- 4
x <- 10^seq(-3, 3, by = .1)
q01 <- 100
y1 <- q01*10^(k*(exp(-alpha*q01*x)-1))
q02 <- 10
y2 <- q02*10^(k*(exp(-alpha*q02*x)-1))
q03 <- 1
y3 <- q03*10^(k*(exp(-alpha*q03*x)-1))
pmax1 <- -VGAM::lambertW(((-1/log(10^(c(k))))))/(alpha*q01)
pmax2 <- -VGAM::lambertW(((-1/log(10^(c(k))))))/(alpha*q02)
pmax3 <- -VGAM::lambertW(((-1/log(10^(c(k))))))/(alpha*q03)

plot(x, y1, log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), pch = 21, 
     main = expression("Untransformed Consumption"), 
     bg = "orange", cex = 2, xlab = NA, ylab = NA, axes = FALSE, cex.main = 1.5)
points(x, y2, pch = 21, bg = "yellow", cex = 1.33)
points(x, y3, pch = 21, bg = "white", cex = 1)
title(ylab = "Consumption", xlab = "Cost per Unit", line = 2.5)
legend("topright", c(expression(list(italic(Q)[0]==100, italic(P)[max]==0.12)),
                     expression(list(italic(Q)[0]==10, italic(P)[max]==1.23)),
                     expression(list(italic(Q)[0]==1, italic(P)[max]==12.28))), pch = 21,
       pt.bg = c("orange", "yellow", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
lines(x = exp(log(pmax3) +c(-log(100000), log(100000))), 
      y = exp(log(q03 * 10^(k*(exp(-alpha*q03*pmax3)- 1)))+ c(log(100000), -log(100000))), lty = 2)



plot(x*q01, y1/q01, log = "xy", ylim = c(.00001, 10), xlim = c(.001, 1000), pch = 21,
     main = expression("Standardized Consumption"), bg = "orange", cex = 2, cex.main = 1.5,
     xlab = NA, ylab = NA, axes = FALSE)
title(ylab = "Standardized Consumption", xlab = "Standardized Cost", line = 2.5)
points(x*q02, y2/q02, pch = 21, bg = "yellow", cex = 1.43)
points(x*q03, y3/q03, pch = 21, bg = "white", cex = .75)
legend("bottomleft", c(expression(list(italic(Q)[0]==100, italic(P)[max]==12.28)),
                       expression(list(italic(Q)[0]==10, italic(P)[max]==12.28)),
                       expression(list(italic(Q)[0]==1, italic(P)[max]==12.28))), pch = 21,
       pt.bg = c("orange", "yellow", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(.000002, .00001, jumps = 5), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:0), labels = c(as.character(c(.001, .01, .1, 1))),
     cex.axis = .9)
axis(2, at = 10^(-5:-4), 
     cex.axis = .9)
axis(2, at = 10^(-4), 
     cex.axis = .9)
lines(x = exp(log(pmax3) +c(-log(100000), log(100000))), 
      y = exp(log(q03 * 10^(k*(exp(-alpha*q03*pmax3)- 1)))+ c(log(100000), -log(100000))), lty = 2)


plot(x, x*y1, log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), pch = 23, 
     main = expression("Untransformed Expenditure"), cex.main = 1.5,
     bg = "orange", cex = 1.66, xlab = NA, ylab = NA,
     frame = TRUE, axes = FALSE)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x, x*y2, pch = 23, bg = "yellow", cex = 1.33)
points(x, x*y3, pch = 23, bg = "white", cex = 1)
abline(h = pmax3*q03*10^(k*(exp(-alpha*pmax3)-1)), lty = 3)
legend("topleft", c(expression(list(italic(Q)[0]==100, italic(O)[max]==4.24)),
                    expression(list(italic(Q)[0]==10, italic(O)[max]==4.24)),
                    expression(list(italic(Q)[0]==1, italic(O)[max]==4.24))), pch = 23,
       pt.bg = c("orange", "yellow", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


plot(x*q01, (x*q01)*(y1/q01), log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), 
     pch = 23, main = expression("Standardized Expenditure"), cex.main = 1.5,
     bg = "orange", cex = 2,
     xlab = NA,
     ylab = NA,
     frame = TRUE, axes = FALSE)
title(ylab = "Expenditure*", xlab = "Standardized Cost", line = 2.5)
points(x*q02, (x*q02)*(y2/q02), pch = 23, bg = "yellow", cex = 1.43)
points(x*q03, (x*q03)*(y3/q03), pch = 23, bg = "white", cex = .75)
legend("topleft", c(expression(list(italic(Q)[0]==100, italic(O)[max]==4.24)),
                    expression(list(italic(Q)[0]==10, italic(O)[max]==4.24)),
                    expression(list(italic(Q)[0]==1, italic(O)[max]==4.24))), pch = 23,
       pt.bg = c("orange", "yellow", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

abline(h = pmax3*q03*10^(k*(exp(-alpha*pmax3)-1)), lty = 3)
dev.off()

# K Differences ----
svg("KDifferences.svg", 7.35, 7.35)
par(mfrow = c(2,2), pty = "s", mar = c(4, 3, 3, .5))

q0 <- 500
k <- 3
alpha <- .006
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", 
      log = "xy", main = expression("Varied "*italic(Q)[0]*", Fixed "*italic(k)*", Fixed "*alpha), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, 
      from = .0001, to = 10000, lty = 4, lwd = 1.5)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)

points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

q0 <- 25
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",add = TRUE, 
      from = .0001, to = 10000, lty = 5, lwd = 1.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)

q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2", add = TRUE, lwd = 2, from = .0001, to = 10000)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500,italic(k)==3,italic(P)[max]==.06)),
                                 expression(list(italic(Q)[0]==25,italic(k)==3,italic(P)[max]==1.15)),
                                 expression(list(italic(Q)[0]==1,italic(k)==3,italic(P)[max]==28.65))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)

q0 <- 500
k <- 3
alpha <- .006
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), 
      col = "orange2", log = "xy", 
      main = expression("Varied "*italic(Q)[0]*", Varied "*italic(k)*", Fixed "*alpha), cex.main = 1.5,
      ylab = NA, xlab = NA,
      from = .0001, to = 10000, axes = FALSE, lty = 4, lwd = 1.5)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)



q0 <- 25
k <- 1.8
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lty = 5, lwd = 1.5)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)


q0 <- 1
k <- 1.2
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax
curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",
      add = TRUE, lwd = 2, from = .0001, to = 10000)
points(x = Pmax, y = q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 3)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500,italic(k)==3,italic(P)[max]==.06)),
                                expression(list(italic(Q)[0]==25,italic(k)==1.8,italic(P)[max]==2.26)),
                                expression(list(italic(Q)[0]==1,italic(k)==1.2,italic(P)[max]==138.31))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)

q0 <- 500
k <- 3
alpha <- .006
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", 
      log = "xy", main = expression("Varied "*italic(Q)[0]*", Fixed "*italic(k)*", Fixed "*alpha), cex.main = 1.5,
      ylab = NA, xlab = NA, 
      from = .0001, to = 10000, axes = FALSE, frame = TRUE, lty = 4, lwd = 1.5)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


q0 <- 25
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lty = 5, lwd = 1.5)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))

q0 <- 1
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2", 
      add = TRUE, lwd = 2, from = .0001, to = 10000)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))
abline(h = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), lty = 3)
legend("bottomright", legend = c(expression(list(italic(Q)[0]==500,italic(k)==3,italic(O)[max]==9.62)),
                                expression(list(italic(Q)[0]==25,italic(k)==3,italic(O)[max]==9.62)),
                                expression(list(italic(Q)[0]==1,italic(k)==3,italic(O)[max]==9.62))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)




q0 <- 500
k <- 3
alpha <- .006
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(.001, 1000), xlim = c(.001, 1000), 
      col = "orange2", log = "xy", main = expression("Varied "*italic(Q)[0]*", Varied "*italic(k)*", Fixed "*alpha), cex.main = 1.5,
      ylab = NA, xlab = NA, 
      from = .0001, to = 10000, axes = FALSE, frame = TRUE, lty = 4, lwd = 1.5)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))
abline(h = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), lty = 3)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)



q0 <- 25
k <- 1.8
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lty = 5, lwd = 1.5)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))
abline(h = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), lty = 3)


q0 <- 1
k <- 1.2
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange2",
      add = TRUE, lwd = 2, from = .0001, to = 10000)
points(x = Pmax, y = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), pch = 4)
Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))
abline(h = Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)), lty = 3)

legend("bottomright", legend = c(expression(list(italic(Q)[0]==500,italic(k)==3,italic(O)[max]==9.62)),
                                 expression(list(italic(Q)[0]==25,italic(k)==1.8,italic(O)[max]==17.16)),
                                 expression(list(italic(Q)[0]==1,italic(k)==1.2,italic(O)[max]==29.12))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "orange2", cex = .9)

dev.off()

# K Corrections ----
svg("KCorrections.svg", 9, 6.5)
par(mfrow = c(2,3), pty = "s", mar = c(4, 4, 4, 1))
q0 <- 25
alpha <- .01
curve((q0 * 10^(1.4 * (exp(-alpha * q0 * x) - 1))),
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main = expression(italic(Q)[0]==25*","~alpha==.01), cex.main = 1.5,
      axes = FALSE,
      ylab = NA, xlab = NA, lty = 2, 
      col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
curve((q0 * 10^(1.2 * (exp(-alpha * q0 * x) - 1))),
     add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(2 * (exp(-alpha * q0 * x) - 1))), 
     add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(3 * (exp(-alpha * q0 * x) - 1))), 
      add = TRUE, lty = 5, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(4 * (exp(-alpha * q0 * x) - 1))), 
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0)
Pmax
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 3)
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(P)[max]==3.32)), 
                  expression(list(italic(k)==1.4, italic(P)[max]==2.09)), 
                  expression(list(italic(k)==2, italic(P)[max]==1.16)),
                  expression(list(italic(k)==3, italic(P)[max]==0.68)),
                  expression(list(italic(k)==4, italic(P)[max]==0.49))),
        box.col = "white", bg = "white", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


alpha <- .01
curve((q0 * 10^(1.4 * (exp(-alpha * q0 * x/1.4) - 1))), 
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main = expression(italic(Q)[0]==25*","~alpha==.01/italic(k)), cex.main = 1.5,
      axes = FALSE, lwd = 1.25,
      ylab = NA, xlab = NA, lty = 2, col = "orange2", from = .0001, to = 10000) 
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
curve((q0 * 10^(1.2 * (exp(-alpha * q0 * x/1.2) - 1))), 
      add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(2 * (exp(-alpha * q0 * x/2) - 1))),
      add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(3 * (exp(-alpha * q0 * x/3) - 1))), 
      add = TRUE, lty = 5, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(4 * (exp(-alpha * q0 * x/4) - 1))),
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0/k)
Pmax
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax/k)- 1))), pch = 3)
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(P)[max]==3.98)), 
                  expression(list(italic(k)==1.4, italic(P)[max]==2.93)), 
                  expression(list(italic(k)==2, italic(P)[max]==2.32)),
                  expression(list(italic(k)==3, italic(P)[max]==2.06)),
                  expression(list(italic(k)==4, italic(P)[max]==1.96))),
       box.col = "white", bg = "white", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
alpha <- .01
curve((q0 * 10^(1.4 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(1.4)))) - 1))), 
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main = expression(italic(Q)[0]==25*","~alpha==.01%*%~-italic(W)(-1/ln(10^italic(k)))), cex.main = 1.5,
      axes = FALSE, lwd = 1.25,
      ylab = NA, xlab = NA, lty = 2, col = "orange2", from = .0001, to = 10000) 
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
curve((q0 * 10^(1.2 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(1.2)))) - 1))), 
      add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(2 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(2)))) - 1))),
      add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(3 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(3)))) - 1))), 
      add = TRUE, lty = 2, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((q0 * 10^(4 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(4)))) - 1))),
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0*-VGAM::lambertW(-1/log(10^(k))))
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax*-VGAM::lambertW(-1/log(10^(k))))- 1))), pch = 3)
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(P)[max]==4)), 
                  expression(list(italic(k)==1.4, italic(P)[max]==4)), 
                  expression(list(italic(k)==2, italic(P)[max]==4)),
                  expression(list(italic(k)==3, italic(P)[max]==4)),
                  expression(list(italic(k)==4, italic(P)[max]==4))),
       bty = "n", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
alpha <- .01
curve((x * q0 * 10^(1.4 * (exp(-alpha * q0 * x) - 1))), 
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main =  expression(italic(Q)[0]==25*","~alpha==.01), cex.main = 1.5,
      axes = FALSE, frame = TRUE, lwd = 1.25,
      ylab = NA, xlab = NA, lty = 2, col = "orange2", from = .0001, to = 10000) 
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
curve((x * q0 * 10^(1.2 * (exp(-alpha * q0 * x) - 1))),
      add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(2 * (exp(-alpha * q0 * x) - 1))), 
      add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(3 * (exp(-alpha * q0 * x) - 1))),
      add = TRUE, lty = 2, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(4 * (exp(-alpha * q0 * x) - 1))),
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0)
points(x = Pmax, y = (Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 4)
(Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(O)[max]==17.47)), 
                  expression(list(italic(k)==1.4, italic(O)[max]==14.07)), 
                  expression(list(italic(k)==2, italic(O)[max]==9.1)),
                  expression(list(italic(k)==3, italic(O)[max]==5.77)),
                  expression(list(italic(k)==4, italic(O)[max]==4.24))),
       box.col = "white", bg = "white", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
alpha <- .01
curve((x * q0 * 10^(1.4 * (exp(-alpha * q0 * x/1.4) - 1))),
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main = expression(italic(Q)[0]==25*","~alpha==.01/italic(k)), cex.main = 1.5,
      axes = FALSE, frame = TRUE, lwd = 1.25,
      ylab = NA, xlab = NA, lty = 2, col = "orange2", from = .0001, to = 10000) 
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
curve((x * q0 * 10^(1.2 * (exp(-alpha * q0 * x/1.2) - 1))), 
      add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(2 * (exp(-alpha * q0 * x/2) - 1))),
      add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(3 * (exp(-alpha * q0 * x/3) - 1))), 
      add = TRUE, lty = 2, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(4 * (exp(-alpha * q0 * x/4) - 1))),
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0/k)
points(x = Pmax, y = (Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax/k)- 1))), pch = 4)
(Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax/k)- 1)))
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(O)[max]==20.97)), 
                  expression(list(italic(k)==1.4, italic(O)[max]==19.7)), 
                  expression(list(italic(k)==2, italic(O)[max]==18.2)),
                  expression(list(italic(k)==3, italic(O)[max]==17.32)),
                  expression(list(italic(k)==4, italic(O)[max]==16.94))),
       bty = "n", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
alpha <- .01
curve((x * q0 * 10^(1.4 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(1.4)))) - 1))),
      ylim = c(.1, 100), xlim = c(.1, 100), log = "xy",
      main = expression(italic(Q)[0]==25*","~alpha==.01%*%~-italic(W)(-1/ln(10^italic(k)))), cex.main = 1.5,
      axes = FALSE, frame = TRUE, lwd = 1.25,
      ylab = NA, xlab = NA, lty = 2, col = "orange2", from = .0001, to = 10000) 
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
curve((x * q0 * 10^(1.2 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(1.2)))) - 1))),
      add = TRUE, lty = 3, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(2 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(2)))) - 1))), 
      add = TRUE, lty = 4, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(3 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(3)))) - 1))), 
      add = TRUE, lty = 2, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
curve((x * q0 * 10^(4 * (exp(-alpha * q0 * x*-VGAM::lambertW(-1/log(10^(4)))) - 1))), 
      add = TRUE, lty = 6, col = "orange2", from = .0001, to = 10000, lwd = 1.25) 
k <- c(1.2, 1.4, 2:4)
Pmax <- -VGAM::lambertW((-1/log(10^(k))))/(alpha*q0*-VGAM::lambertW(-1/log(10^(k))))
points(x = Pmax, y = (Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax*-VGAM::lambertW(-1/log(10^(k))))- 1))), pch = 4)
(Pmax * q0 * 10^(k*(exp(-alpha*q0*Pmax*-VGAM::lambertW(-1/log(10^(k))))- 1)))
legend("bottomleft", 
       legend = c(expression(list(italic(k)==1.2, italic(O)[max]==21.05)), 
                  expression(list(italic(k)==1.4, italic(O)[max]==26.87)), 
                  expression(list(italic(k)==2, italic(O)[max]==31.34)),
                  expression(list(italic(k)==3, italic(O)[max]==33.58)),
                  expression(list(italic(k)==4, italic(O)[max]==34.51))),
       bty = "n", lty = c(3, 2, 4, 5, 6), col = "orange2", lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
dev.off()

# Samuelson Normalized Decay ----
svg("SamuelsonNormalizedDecay.svg", 7.35, 7.35)
par(mfrow = c(2,2), pty = "s", mar = c(4, 3, 3, .5))

alpha <- .2
q0 <- 25
Pmax <- 1/(alpha)
Pmax
curve((q0 * exp(-alpha* x)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "blue", log = "xy", 
      main = expression("Non-Normalized\n     Consumption"~":"~italic(Q)[0] %*%~ italic(e) ^ (-ring(alpha)*italic(C))),
      ylab = NA, xlab = NA, from = .0001, to = 10000, cex.main = 1.35,
      frame = FALSE, axes = FALSE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = (q0 * exp(-alpha * Pmax)))
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * exp(-alpha* Pmax))+ c(log(100000), -log(100000))), lty = 2, type = "l")
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


q0 <- 500
curve((q0 * exp(-alpha*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
Pmax <- 1/(alpha)
Pmax
points(x = Pmax, y = (q0 * exp(-alpha* Pmax)))
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * exp(-alpha* Pmax))+ c(log(100000), -log(100000))), lty = 2)

q0 <- 1
curve((q0 * exp(-alpha* x)), col = "blue",
      add = TRUE, from = .0001, to = 10000, lwd = 2)
Pmax <- 1/(alpha)
points(x = Pmax, y = (q0 * exp(-alpha* Pmax)))
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * exp(-alpha* Pmax))+ c(log(100000), -log(100000))), lty = 2)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500, italic(P)[max]==5)),
                                expression(list(italic(Q)[0]==25, italic(P)[max]==5)),
                                expression(list(italic(Q)[0]==1, italic(P)[max]==5))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "blue", cex = .9)



q0 <- 25
Pmax <- 1/(alpha*q0)
Pmax
curve((q0 * exp(-alpha*q0*x)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "blue", log = "xy",
      main = expression("   Normalized\nConsumption"~":"~italic(Q)[0] %*%~ italic(e) ^ (-alpha*italic(Q)[0]*italic(C))),
      ylab = NA, xlab = NA, 
      from = .0001, to = 10000, cex.main = 1.35, axes = FALSE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * exp(-alpha*q0*Pmax))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)



q0 <- 500
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
Pmax <- 1/(alpha*q0)
Pmax
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 1
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 2)
Pmax <- 1/(alpha*q0)
Pmax
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

legend("topright", legend = c(expression(list(italic(Q)[0]==500, italic(P)[max]==0.01)),
                                expression(list(italic(Q)[0]==25, italic(P)[max]==.2)),
                                expression(list(italic(Q)[0]==1, italic(P)[max]==5))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "blue", cex = .9)


q0 <- 25
Pmax <- 1/(alpha)
curve(x * (q0 * exp(-alpha* x)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "blue", log = "xy", 
      main = expression("Non-Normalized\n       Expenditure"~":"~italic(C)%*%~italic(Q)[0] %*%~ italic(e) ^ (-ring(alpha)*italic(C))),
      ylab = NA, xlab = NA, from = .0001, to = 10000, cex.main = 1.35,
      frame = TRUE, axes = FALSE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax *(q0 * exp(-alpha * Pmax)), pch = 5)
Pmax *(q0 * exp(-alpha * Pmax))
abline(h = Pmax *(q0 * exp(-alpha * Pmax)), lty = 3)


q0 <- 500
curve(x * (q0 * exp(-alpha*x)), col = "blue",
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
Pmax <- 1/(alpha)
points(x = Pmax, y = Pmax * (q0 * exp(-alpha* Pmax)), pch = 5)
Pmax *(q0 * exp(-alpha * Pmax))
abline(h = Pmax *(q0 * exp(-alpha * Pmax)), lty = 3)


q0 <- 1
curve(x * (q0 * exp(-alpha* x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 2)
Pmax <- 1/(alpha)
points(x = Pmax, y = Pmax *(q0 * exp(-alpha* Pmax)), pch = 5)
Pmax *(q0 * exp(-alpha * Pmax))
abline(h = Pmax *(q0 * exp(-alpha * Pmax)), lty = 3)

legend("bottomleft", legend = c(expression(list(italic(Q)[0]==500, italic(O)[max]==919.7)),
                             expression(list(italic(Q)[0]==25, italic(O)[max]==45.98)),
                             expression(list(italic(Q)[0]==1, italic(O)[max]==1.84))
), box.col = "white", bg = "white", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "blue", cex = .9)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-4:3), labels = as.character(c(NA, .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


q0 <- 25
Pmax <- 1/(alpha*q0)
curve(x * (q0 * exp(-alpha*q0*x)), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "blue", log = "xy",
      main = expression(" Normalized\nExpenditure"~":"~italic(C) %*%~italic(Q)[0] %*%~ italic(e) ^ (-alpha*italic(Q)[0]*italic(C))),
      ylab = NA, xlab = NA, from = .0001, to = 10000, cex.main = 1.35,
      axes = FALSE, frame = TRUE, lty = 5, lwd = 1.5, line = 1.25)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax * (q0 * exp(-alpha*q0*Pmax)), pch = 5)
Pmax * (q0 * exp(-alpha*q0*Pmax))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)


q0 <- 500
curve(x * (q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lty = 4, lwd = 1.5)
Pmax <- 1/(alpha*q0)
points(x = Pmax, y = Pmax * (q0 * exp(-alpha*q0*Pmax)), pch = 5)
Pmax * (q0 * exp(-alpha*q0*Pmax))

q0 <- 1
curve(x * (q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 2)
Pmax <- 1/(alpha*q0)
points(x = Pmax, y = Pmax * (q0 * exp(-alpha*q0*Pmax)), pch = 5)
Pmax * (q0 * exp(-alpha*q0*Pmax))
abline(h = Pmax * (q0 * exp(-alpha*q0*Pmax)), lty = 3)
legend("topright", legend = c(expression(list(italic(Q)[0]==500, italic(O)[max]==1.84)),
                                expression(list(italic(Q)[0]==25, italic(O)[max]==1.84)),
                                expression(list(italic(Q)[0]==1, italic(O)[max]==1.84))
), bty = "n", lty = c(4, 5, 1), lwd = c(1.5, 1.5, 2), col = "blue", cex = .9)

dev.off()
# SND standardized----
svg("SNDStandardized.svg", 7.35, 7.35)
par(mfrow = c(2,2), pty = "s", mar = c(4, 3, 3, .5))
alpha <- .08
x <- 10^seq(-3, 3, by = .1)
q01 <- 100
y1 <- q01*exp(-alpha*q01*x)
q02 <- 10
y2 <- q02*exp(-alpha*q02*x)
q03 <- 1
y3 <- q03*exp(-alpha*q03*x)

plot(x, y1, log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), pch = 21,
     main = expression("Untransformed Consumption"), cex.main = 1.5, 
     bg = "blue", cex = 1.66, xlab = NA, ylab = NA, axes = FALSE)
title(ylab = "Consumption", xlab = "Cost per Unit", line = 2.5)
points(x, y2, pch = 21, bg = "lightblue", cex = 1.33)
points(x, y3, pch = 21, bg = "white", cex = 1)
1/(alpha*c(q01,q02,q03))
legend("topright", c(expression(list(italic(Q)[0]==100, italic(P)[max]==.13)),
                     expression(list(italic(Q)[0]==10, italic(P)[max]==1.25)),
                     expression(list(italic(Q)[0]==1, italic(P)[max]==12.5))), pch = 21,
       pt.bg = c("blue", "lightblue", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
lines(x = exp(log(1/alpha) +c(-log(100000), log(100000))), 
      y = exp(log(q03 * exp(-alpha* q03 * (1/alpha)))+ c(log(100000), -log(100000))), lty = 2, type = "l")



plot(x*q01, y1/q01, log = "xy", ylim = c(.00001, 10), xlim = c(.001, 1000), pch = 21,
     main = expression("Standardized Consumption"), cex.main = 1.5, bg = "blue", cex = 2,
     xlab = NA, ylab = NA, axes = FALSE)
title(ylab = "Standardized Consumption", xlab = "Standardized Cost", line = 2.5)
points(x*q02, y2/q02, pch = 21, bg = "lightblue", cex = 1.43)
points(x*q03, y3/q03, pch = 21, bg = "white", cex = .75)
legend("bottomleft", c(expression(list(italic(Q)[0]==100, italic(P)[max]==12.5)),
                       expression(list(italic(Q)[0]==10, italic(P)[max]==12.5)),
                       expression(list(italic(Q)[0]==1, italic(P)[max]==12.5))), pch = 21,
       pt.bg = c("blue", "lightblue", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(.000002, .00001, jumps = 5), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:0), labels = c(as.character(c(.001, .01, .1, 1))), 
     cex.axis = .9)
axis(2, at = 10^(-5:-4), 
     cex.axis = .9)
axis(2, at = 10^(-4), 
     cex.axis = .9)
lines(x = exp(log(1/alpha) +c(-log(100000), log(100000))), 
      y = exp(log(q03 * exp(-alpha* q03 * (1/alpha)))+ c(log(100000), -log(100000))), lty = 2, type = "l")

plot(x, x*y1, log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), pch = 23, 
     main = expression("Untransformed Expenditure"), cex.main = 1.5,
     bg = "blue", cex = 1.66, xlab = NA, ylab = NA,
     frame = TRUE, axes = FALSE)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x, x*y2, pch = 23, bg = "lightblue", cex = 1.33)
points(x, x*y3, pch = 23, bg = "white", cex = 1)
abline(h = 1/(alpha*exp(1)), lty = 3)
1/(alpha*exp(1))
legend("topleft", c(expression(list(italic(Q)[0]==100, italic(O)[max]==4.6)),
                    expression(list(italic(Q)[0]==10, italic(O)[max]==4.6)),
                    expression(list(italic(Q)[0]==1, italic(O)[max]==4.6))), pch = 23,
       pt.bg = c("blue", "lightblue", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

plot(x*q01, (x*q01)*(y1/q01), log = "xy", ylim = c(.001, 1000), xlim = c(.001, 1000), 
     pch = 23, main = expression("Standardized Expenditure"), cex.main = 1.5,
     bg = "blue", cex = 2,
     xlab = NA,
     ylab = NA,
     frame = TRUE, axes = FALSE)
title(ylab = "Expenditure*", xlab = "Standardized Cost", line = 2.5)
points(x*q02, (x*q02)*(y2/q02), pch = 23, bg = "lightblue", cex = 1.43)
points(x*q03, (x*q03)*(y3/q03), pch = 23, bg = "white", cex = .75)
abline(h = 1/(alpha*exp(1)), lty = 3)
legend("topleft", 
       c(expression(list(italic(Q)[0]==100, italic(O)[max]==4.6)),
         expression(list(italic(Q)[0]==10, italic(O)[max]==4.6)),
         expression(list(italic(Q)[0]==1, italic(O)[max]==4.6))),
       pch = 23,
       pt.bg = c("blue", "lightblue", "white"), bty = "n", pt.cex = c(2, 1.5, 1))
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)
axis(2, at = 10^(-3:3), labels = as.character(c( .001, .01, .1, 1, 10, 100, 1000)), cex.axis = .9)

dev.off()

# Mathematical Conversions for between EXPD and SND ----
svg("EXPDtoSNDConversions.svg", 9, 6.5)
par(mfrow = c(2,3), pty = "s", mar = c(4, 4, 4, 1))

q0 <- 500
k <- 7
alpha <- .1
Pmax <- 1/(alpha*q0*(EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), 
      col = "orange2", log = "xy", lwd = 1.25,
      main = expression(list(alpha[HS]==.1,alpha[SND]==.1)), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, from = .0001, to = 10000, lty = 4)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))

q0 <- 25
k <- 3
Pmax <- 1/(alpha*q0 *(EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), 
      col = "orange2", add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 20)

q0 <- 1
k <- 1
Pmax <- 1/(alpha*q0 *(EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))), pch = 20)

q0 <- 500
alpha <- .1
Pmax <- 1/(q0 * alpha)
curve((q0 * exp(-alpha*q0*x)), col = "blue",
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 25
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 1
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue",
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

legend("topright", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                              expression(list(italic(Q)[0]==25,italic(k)==3)),
                              expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)



q0 <- 500
k <-7
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0*(EXPDToSND(k)))

curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression(list(alpha[HS]==.1/italic(F)[italic(b)], alpha[SND]==.1)), cex.main = 1.5,
      ylab = NA, xlab = NA, 
      axes = FALSE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20, lwd = 1.25)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25

k <- 3
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0 * (EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)


q0 <- 1
k <- 1
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0 * (EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))),
      col = "orange2", add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)


q0 <- 500
alpha <- .1
Pmax <- 1/(q0 * alpha)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 25
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 1
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

legend("topright", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                              expression(list(italic(Q)[0]==25,italic(k)==3)),
                              expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)

q0 <- 500
k <-7
alpha <- .1

Pmax <- 1/(alpha*q0*(EXPDToSND(k)))

curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression(list(alpha[HS]==.1, alpha[SND]==.1%*%~italic(F)[italic(b)])), cex.main = 1.5,
      ylab = NA, xlab = NA, 
      from = .0001, to = 10000, axes = FALSE, lwd = 1.25, lty = 4)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
k <- 3
Pmax <- 1/(alpha*q0 *(EXPDToSND(k)))

curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)

q0 <- 1
k <- 1
Pmax <- 1/(alpha*q0 *(EXPDToSND(k)))
curve((q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 20)

q0 <- 500
k <-7
alpha <- .1 * (EXPDToSND(k))
Pmax <- 1/(q0 * alpha)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))


q0 <- 25
k <- 3

alpha <- .1 * (EXPDToSND(k))
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))

q0 <- 1
k <- 1
alpha <- .1 * (EXPDToSND(k))
Pmax <- 1/(alpha*q0)
curve((q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = (q0 * exp(-alpha*q0*Pmax)))
legend("topright", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                              expression(list(italic(Q)[0]==25,italic(k)==3)),
                              expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)


q0 <- 500
k <- 7
alpha <- .1
Pmax <- 1/(alpha*q0*(EXPDToSND(k)))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression(list(alpha[HS]==.1, alpha[SND]==.1)), cex.main = 1.5,
      ylab = NA, xlab = NA,
      from = .0001, to = 10000, frame = TRUE, axes = FALSE, lwd = 1.25, lty = 4)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
k <- 3
alpha <- .1
Pmax <- 1/(alpha*q0 *(EXPDToSND(k)))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))),
      col = "orange2", add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)

q0 <- 1
k <- 1
alpha <- .1
Pmax <- 1/(alpha*q0 * (EXPDToSND(k)))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)

q0 <- 500
alpha <- .1
Pmax <- 1/(q0 * alpha)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)

q0 <- 25
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)


q0 <- 1
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)

legend("topleft", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                              expression(list(italic(Q)[0]==25,italic(k)==3)),
                              expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)



q0 <- 500
k <-7
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0*EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy",
      main = expression(list(alpha[HS]==.1/italic(F)[italic(b)], alpha[SND]==.1)), cex.main = 1.5,
      ylab = NA, xlab = NA,
      from = .0001, to = 10000, frame = TRUE, axes = FALSE, lwd = 1.25, lty = 4)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))


q0 <- 25
k <- 3
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0 * EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)

q0 <- 1
k <- 1
alpha <- .1/(EXPDToSND(k))
Pmax <- 1/(alpha*q0 * EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)

q0 <- 500
alpha <- .1
Pmax <- 1/(q0 * alpha)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue",
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)


q0 <- 25
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)


q0 <- 1
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)),  col = "blue", 
      add = TRUE, from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)

legend("topleft", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                             expression(list(italic(Q)[0]==25,italic(k)==3)),
                             expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)




q0 <- 500
k <-7
alpha <- .1

Pmax <- 1/(alpha*q0*EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), ylim = c(.001, 1000), xlim = c(.001, 1000), col = "orange2", log = "xy", 
      main = expression(list(alpha[HS]==.1, alpha[SND]==.1%*%~italic(F)[italic(b)])), cex.main = 1.5,
      ylab = NA, xlab = NA, 
      from = .0001, to = 10000, axes = FALSE, frame = TRUE, lwd = 1.25, lty = 4)
title(ylab = "Expenditure", xlab = "Cost per Unit", line = 2.5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))
axis(2, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)))



q0 <- 25
k <- 3
alpha <- .1

Pmax <- 1/(alpha*q0 * EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), col = "orange2",
      add = TRUE, from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)


q0 <- 1
k <- 1
alpha <- .1

Pmax <- 1/(alpha*q0 * EXPDToSND(k))
curve(x*(q0 * 10^(k*(exp(-alpha*q0*x)- 1))), 
      add = TRUE, col = "orange2", from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1))),pch = 18)


q0 <- 500
k <-7
alpha <- .1 * EXPDToSND(k)
Pmax <- 1/(q0 * alpha)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", add = TRUE, 
      from = .0001, to = 10000, lwd = 1.25, lty = 4)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)

q0 <- 25
k <- 3

alpha <- .1 * EXPDToSND(k)
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", add = TRUE,
      from = .0001, to = 10000, lwd = 1.25, lty = 5)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)


q0 <- 1
k <- 1
alpha <- .1 * EXPDToSND(k)
Pmax <- 1/(alpha*q0)
curve(x*(q0 * exp(-alpha*q0*x)), col = "blue", add = TRUE,
      from = .0001, to = 10000, lwd = 1.25)
points(x = Pmax, y = Pmax*(q0 * exp(-alpha*q0*Pmax)), pch = 5)

legend("topleft", legend = c(expression(list(italic(Q)[0]==500,italic(k)==7)),
                             expression(list(italic(Q)[0]==25,italic(k)==3)),
                             expression(list(italic(Q)[0]==1,italic(k)==1))
), bty = "n", lty = c(4, 5, 1), lwd = 1.25, col = "black", cex = 1)


dev.off()


# Model fits on Data from Both Models ---- 

# Koffarnus et al. 2015 Simulations ----
set.seed(09264)
koffarnus2015SimsData <- read.csv("koffarnus2015Sims.csv")

koffarnus2015SimsClean <- koffarnus2015SimsData |> 
    pivot_longer(cols = -id, names_to = "x", values_to = "y") |> na.omit()

koffarnus2015SimsClean$x <- as.numeric(gsub("X", "", koffarnus2015SimsClean$x))
koffarnus2015SimsClean$k <- GetK(koffarnus2015SimsClean)

koffarnus2015SimsEXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = koffarnus2015SimsClean,
                                          iter = c(10, 10),
                                          start_lower = c(q0 = .1, alpha = -5),
                                          start_upper = c(q0 = 2, alpha = 0))

koffarnus2015SimsEXPDNLME <- nlme(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), 
                                  data = koffarnus2015SimsClean,
                                  fixed = list(q0 ~ 1,  
                                               alpha ~ 1),
                                  random = list(pdSymm(q0 + alpha ~ 1)),
                                  start = list(fixed = coef(koffarnus2015SimsEXPDNLS)), 
                                  groups = ~id,
                                  method = "ML",
                                  verbose = 2,
                                  control = list(msMaxIter = 5000,
                                                 niterEM = 5000,
                                                 maxIter = 5000,
                                                 pnlsTol = .0001,
                                                 tolerance = .001,
                                                 apVar = T,
                                                 minScale = .0000001,
                                                 opt = "optim"))


koffarnus2015SimsSNDNLS <- nls_multstart(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                         data = koffarnus2015SimsClean,
                                         iter = c(10, 10),
                                         start_lower = c(q0 = .1, alpha = -5),
                                         start_upper = c(q0 = 2, alpha = 0))

koffarnus2015SimsSNDNLME <- nlme(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                 data = koffarnus2015SimsClean,
                                 fixed = list(q0 ~ 1,  
                                              alpha ~ 1),
                                 random = list(pdSymm(q0 + alpha ~ 1)),
                                 start = list(fixed = coef(koffarnus2015SimsSNDNLS)), 
                                 groups = ~id,
                                 method = "ML",
                                 verbose = 2,
                                 control = list(msMaxIter = 5000,
                                                niterEM = 5000,
                                                maxIter = 5000,
                                                pnlsTol = .0001,
                                                tolerance = .001,
                                                apVar = T,
                                                minScale = .0000001,
                                                opt = "optim"), na.action = na.omit)


anova(koffarnus2015SimsEXPDNLME, koffarnus2015SimsSNDNLME)

# Koffarnus et al. (2015) CPT ----
set.seed(09264)
koffarnus2015CPTData <- read.csv("koffarnus2015CPT.csv")

koffarnus2015CPTClean <- koffarnus2015CPTData |> 
    pivot_longer(cols = -subject_id, names_to = "x", values_to = "y") |> na.omit()
koffarnus2015CPTClean$x <- as.numeric(gsub("X", "", koffarnus2015CPTClean$x))
koffarnus2015CPTClean$k <- GetK(koffarnus2015CPTClean)

koffarnus2015CPTEXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), 
                                         data = koffarnus2015CPTClean,
                                         iter = c(10, 10),
                                         start_lower = c(q0 = .1, alpha = -5),
                                         start_upper = c(q0 = 2, alpha = 0))

koffarnus2015CPTEXPDNLME <- nlme(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), 
                                 data = koffarnus2015CPTClean,
                                 fixed = list(q0 ~ 1,  
                                              alpha ~ 1),
                                 random = list(pdSymm(q0 + alpha ~ 1)),
                                 start = list(fixed = coef(koffarnus2015CPTEXPDNLS)), 
                                 groups = ~subject_id,
                                 method = "ML",
                                 verbose = 2,
                                 control = list(msMaxIter = 5000,
                                                niterEM = 5000,
                                                maxIter = 5000,
                                                pnlsTol = .0001,
                                                tolerance = .001,
                                                apVar = T,
                                                minScale = .0000001,
                                                opt = "optim"))

koffarnus2015CPTSNDNLS <- nls_multstart(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                        data = koffarnus2015CPTClean,
                                         iter = c(10, 10),
                                         start_lower = c(q0 = .1, alpha = -5),
                                         start_upper = c(q0 = 2, alpha = 0))

koffarnus2015CPTSNDNLME <- nlme(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                 data = koffarnus2015CPTClean,
                                 fixed = list(q0 ~ 1,  
                                              alpha ~ 1),
                                 random = list(pdSymm(q0 + alpha ~ 1)),
                                 start = list(fixed = c(1.3,-2)),
                                 # Comment the above "start" line and uncomment the below "start" line
                                 # to see the singularity error  
                                 # start = list(fixed = coef(koffarnus2015CPTSNDNLS)),
                                 groups = ~subject_id,
                                 method = "ML",
                                 verbose = 2,
                                control = list(msMaxIter = 5000,
                                               niterEM = 5000,
                                               maxIter = 5000,
                                               pnlsTol = .0001,
                                               tolerance = .001,
                                               apVar = T,
                                               minScale = .0000001,
                                               opt = "optim"))


koffarnus2015CPTSNDNLMEalt <- nlme(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                data = koffarnus2015CPTClean,
                                fixed = list(q0 ~ 1,  
                                             alpha ~ 1),
                                random = list(pdSymm(q0 + alpha ~ 1)),
                                start = list(fixed = c(1.3,-1.8)),
                                groups = ~subject_id,
                                method = "ML",
                                verbose = 2,
                                control = list(msMaxIter = 5000,
                                               niterEM = 5000,
                                               maxIter = 5000,
                                               pnlsTol = .0001,
                                               tolerance = .001,
                                               apVar = T,
                                               minScale = .0000001,
                                               opt = "optim"))

anova(koffarnus2015CPTEXPDNLME, koffarnus2015CPTSNDNLME, koffarnus2015CPTSNDNLMEalt)

# Kaplan & Reed 2018 Data ----
set.seed(09264)
kaplan2018Data <- read.csv("kaplan2018Data.csv")

kaplan2018Clean <- kaplan2018Data |> 
    select(starts_with("x"))  |>  
    rename(id = X) |>  
    pivot_longer(cols = -id, names_to = "x", values_to = "y")

kaplan2018Clean$x <- as.numeric(gsub("X", "", kaplan2018Clean$x))
kaplan2018Clean <- kaplan2018Clean[complete.cases(kaplan2018Clean), ]
kaplan2018Clean$k <- GetK(kaplan2018Clean)

kaplan2018EXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = kaplan2018Clean,
                                   iter = c(10, 10),
                                   start_lower = c(q0 = .1, alpha = -5),
                                   start_upper = c(q0 = 2, alpha = 0))

kaplan2018EXPDNLME <- nlme(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), 
                           data = kaplan2018Clean,
                           fixed = list(q0 ~ 1,  
                                        alpha ~ 1),
                           random = list(pdSymm(q0 + alpha ~ 1)),
                           start = list(fixed = coef(kaplan2018EXPDNLS)), 
                           groups = ~id,
                           method = "ML",
                           verbose = 2,
                           control = list(msMaxIter = 5000,
                                          niterEM = 5000,
                                          maxIter = 5000,
                                          pnlsTol = .0001,
                                          tolerance = .001,
                                          apVar = T,
                                          minScale = .0000001,
                                          opt = "optim"))

kaplan2018SNDNLS <- nls_multstart(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                                  data = kaplan2018Clean,
                                  iter = c(10, 10),
                                  start_lower = c(q0 = .1, alpha = -5),
                                  start_upper = c(q0 = 2, alpha = 0))


kaplan2018SNDNLME <- nlme(y ~ 10^(q0) * exp(-10^(alpha) * 10^(q0) * x),
                          data = kaplan2018Clean,
                          fixed = list(q0 ~ 1,  
                                       alpha ~ 1),
                          random = list(pdSymm(q0 + alpha ~ 1)),
                          start = list(fixed = coef(kaplan2018SNDNLS)), 
                          groups = ~id,
                          method = "ML",
                          verbose = 2,
                          control = list(msMaxIter = 5000,
                                         niterEM = 5000,
                                         maxIter = 5000,
                                         pnlsTol = .0001,
                                         tolerance = .001,
                                         apVar = T,
                                         minScale = .0000001,
                                         opt = "optim"), na.action = na.omit)


anova(kaplan2018EXPDNLME, kaplan2018SNDNLME)

# Koffarnus 2015 and Kaplan 2018 Individual Plots----

svg("combinedIndividualPlots.svg", 8, 8)


layout(matrix(c(rep(1,3), 2, 3, 4, 
         rep(5, 3), 6, 7, 8, 
         rep(9, 3), 10, 11, 12), ncol = 3, byrow = TRUE), heights = c(1, 6, 1, 6, 1, 6))


par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) Simulations", cex = 2.5)

par(mar = c(5, 4, 2, .5), pty = "m")


for(a in unique(koffarnus2015SimsClean$id)[1:3]){
    koffarnus2015SimsEXPDQ0 <- coef(koffarnus2015SimsEXPDNLME)[which(unique(koffarnus2015SimsClean$id) == a),1]
    koffarnus2015SimsEXPDAlpha <- coef(koffarnus2015SimsEXPDNLME)[which(unique(koffarnus2015SimsClean$id) == a),2]
    koffarnus2015SimsSNDQ0 <- coef(koffarnus2015SimsSNDNLME)[which(unique(koffarnus2015SimsClean$id) == a),1]
    koffarnus2015SimsSNDAlpha <- coef(koffarnus2015SimsSNDNLME)[which(unique(koffarnus2015SimsClean$id) == a),2]
    curve(10^(koffarnus2015SimsEXPDQ0) * 10^(3.994711 * (exp(-10^(koffarnus2015SimsEXPDAlpha) * 10^(koffarnus2015SimsEXPDQ0) * x) - 1)), 
          xlim = c(.01, 1200), ylim = c(0, 25), log = "x", from = .001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Unit", ylab = "Simulated Purchases", main = paste("ID:", a), lwd = 1.25)
    curve(10^koffarnus2015SimsSNDQ0 * exp(- 10^koffarnus2015SimsSNDAlpha * 10^koffarnus2015SimsSNDQ0 * x), 
          from = .001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(koffarnus2015SimsClean, id == a)$x == 0, .01, 
                      subset(koffarnus2015SimsClean, id == a)$x),
           y = subset(koffarnus2015SimsClean, id == a)$y)
    axis(1, at = tickFunc(from = .01, to = .1, jumps = 4), tcl = -.25, labels = NA)
    axis(1, at = 10^(-1:3), labels = c("0.1","1","10","100", "1000"))
    axis(1, at = .01, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}

par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) CPT", cex = 2.5)

par(mar = c(5, 4, 2, .5), pty = "m")



for(a in unique(koffarnus2015CPTClean$subject_id)[1:3]){
    koffarnus2015CPTEXPDQ0 <- coef(koffarnus2015CPTEXPDNLME)[which(rownames(coef(koffarnus2015CPTEXPDNLME)) == a),1]
    koffarnus2015CPTEXPDAlpha <- coef(koffarnus2015CPTEXPDNLME)[which(rownames(coef(koffarnus2015CPTEXPDNLME)) == a),2]
    koffarnus2015CPTSNDQ0 <- coef(koffarnus2015CPTSNDNLME)[which(rownames(coef(koffarnus2015CPTSNDNLME)) == a),1]
    koffarnus2015CPTSNDAlpha <- coef(koffarnus2015CPTSNDNLME)[which(rownames(coef(koffarnus2015CPTSNDNLME)) == a),2]
    curve(10^(koffarnus2015CPTEXPDQ0) * 10^(3.683697 * (exp(-10^(koffarnus2015CPTEXPDAlpha) * 10^(koffarnus2015CPTEXPDQ0) * x) - 1)), 
          xlim = c(.01, 1200), ylim = c(0, 25), log = "x", from = .0001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Cigarette", ylab = "Purchases", main = paste("ID:",a), lwd = 1.25)
    curve(10^koffarnus2015CPTSNDQ0 * exp(- 10^koffarnus2015CPTSNDAlpha * 10^koffarnus2015CPTSNDQ0 * x), 
          from = .0001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(koffarnus2015CPTClean, subject_id == a)$x == 0, .01, 
                      subset(koffarnus2015CPTClean, subject_id == a)$x),
           y = subset(koffarnus2015CPTClean, subject_id == a)$y)
    axis(1, at = tickFunc(from = .01, to = .1, jumps = 4), tcl = -.25, labels = NA)
    axis(1, at = 10^(-1:3), labels = c("0.1","1","10","100", "1000"))
    axis(1, at = .01, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}

par(mar = rep(0,4))
plot.new()
text(.5, .5, "Kaplan & Reed (2018) APT", cex = 2.5)

par(mar = c(5, 4, 2, .5), pty = "m")
for(a in unique(kaplan2018Clean$id)[1:3]){
    kaplan2018EXPDQ0 <- coef(kaplan2018EXPDNLME)[which(unique(kaplan2018Clean$id) == a),1]
    kaplan2018EXPDAlpha <- coef(kaplan2018EXPDNLME)[which(unique(kaplan2018Clean$id) == a),2]
    kaplan2018SNDQ0 <- coef(kaplan2018SNDNLME)[which(unique(kaplan2018Clean$id) == a),1]
    kaplan2018SNDAlpha <- coef(kaplan2018SNDNLME)[which(unique(kaplan2018Clean$id) == a),2]
    curve(10^(kaplan2018EXPDQ0) * 10^(1.722925 * (exp(-10^(kaplan2018EXPDAlpha) * 10^(kaplan2018EXPDQ0) * x) - 1)), 
          xlim = c(.01, 1200), ylim = c(0, 15), log = "x", from = .001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Drink", ylab = "Purchases", main = paste("ID:",a), lwd = 1.25)
    curve(10^kaplan2018SNDQ0 * exp(- 10^kaplan2018SNDAlpha * 10^kaplan2018SNDQ0 * x), 
          from = .001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(kaplan2018Clean, id == a)$x == 0, .01, 
                      subset(kaplan2018Clean, id == a)$x),
           y = subset(kaplan2018Clean, id == a)$y)
    axis(1, at = tickFunc(from = .01, to = .1, jumps = 4), tcl = -.25, labels = NA)
    axis(1, at = 10^(-1:3), labels = c("0.1","1","10","100", "1000"))
    axis(1, at = .01, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}
dev.off()


# Additional Random Effects Fits----

svg("koffarnus2015sim.svg", 11, 11)
layout(matrix(c(rep(1, 5), 2:26), byrow = TRUE, ncol = 5), heights = c(.05, rep((1-.05)/5, 5)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) Simulations", cex = 3)
par(mar = c(5, 4, 2, .5), pty = "m")

for(a in unique(koffarnus2015SimsClean$id)[4:28]){
    koffarnus2015SimsEXPDQ0 <- coef(koffarnus2015SimsEXPDNLME)[which(unique(koffarnus2015SimsClean$id) == a),1]
    koffarnus2015SimsEXPDAlpha <- coef(koffarnus2015SimsEXPDNLME)[which(unique(koffarnus2015SimsClean$id) == a),2]
    koffarnus2015SimsSNDQ0 <- coef(koffarnus2015SimsSNDNLME)[which(unique(koffarnus2015SimsClean$id) == a),1]
    koffarnus2015SimsSNDAlpha <- coef(koffarnus2015SimsSNDNLME)[which(unique(koffarnus2015SimsClean$id) == a),2]
    curve(10^(koffarnus2015SimsEXPDQ0) * 10^(3.994711 * (exp(-10^(koffarnus2015SimsEXPDAlpha) * 10^(koffarnus2015SimsEXPDQ0) * x) - 1)), 
          xlim = c(.001, 1200), ylim = c(0, 25), log = "x", from = .0001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Cigarette", ylab = "Simulated Consumption", main = paste("ID:",a), lwd = 1.25)
    curve(10^koffarnus2015SimsSNDQ0 * exp(- 10^koffarnus2015SimsSNDAlpha * 10^koffarnus2015SimsSNDQ0 * x), 
          from = .0001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(koffarnus2015SimsClean, id == a)$x == 0, .001, 
                      subset(koffarnus2015SimsClean, id == a)$x),
           y = subset(koffarnus2015SimsClean, id == a)$y)
    axis(1, at = tickFunc(), tcl = -.25, labels = NA)
    axis(1, at = 10^(-2:4))
    axis(1, at = .001, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}
dev.off()

svg("koffarnus2015cpt.svg", 11, 11)
layout(matrix(c(rep(1, 5), 2:26), byrow = TRUE, ncol = 5), heights = c(.05, rep((1-.05)/5, 5)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) CPT", cex = 3)
par(mar = c(5, 4, 2, .5), pty = "m")
for(a in unique(koffarnus2015CPTClean$subject_id)[4:28]){
    koffarnus2015CPTEXPDQ0 <- coef(koffarnus2015CPTEXPDNLME)[which(rownames(coef(koffarnus2015CPTEXPDNLME)) == a),1]
    koffarnus2015CPTEXPDAlpha <- coef(koffarnus2015CPTEXPDNLME)[which(rownames(coef(koffarnus2015CPTEXPDNLME)) == a),2]
    koffarnus2015CPTSNDQ0 <- coef(koffarnus2015CPTSNDNLME)[which(rownames(coef(koffarnus2015CPTSNDNLME)) == a),1]
    koffarnus2015CPTSNDAlpha <- coef(koffarnus2015CPTSNDNLME)[which(rownames(coef(koffarnus2015CPTSNDNLME)) == a),2]
    curve(10^(koffarnus2015CPTEXPDQ0) * 10^(3.683697 * (exp(-10^(koffarnus2015CPTEXPDAlpha) * 10^(koffarnus2015CPTEXPDQ0) * x) - 1)), 
          xlim = c(.001, 1200), ylim = c(0, 25), log = "x", from = .0001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Cigarette", ylab = "Purchases", main = paste("ID:",a), lwd = 1.25)
    curve(10^koffarnus2015CPTSNDQ0 * exp(- 10^koffarnus2015CPTSNDAlpha * 10^koffarnus2015CPTSNDQ0 * x), 
          from = .0001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(koffarnus2015CPTClean, subject_id == a)$x == 0, .001, 
                      subset(koffarnus2015CPTClean, subject_id == a)$x),
           y = subset(koffarnus2015CPTClean, subject_id == a)$y)
    axis(1, at = tickFunc(), tcl = -.25, labels = NA)
    axis(1, at = 10^(-2:4))
    axis(1, at = .001, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}
dev.off()


svg("kaplan2018.svg", 11, 11)
layout(matrix(c(rep(1, 5), 2:26), byrow = TRUE, ncol = 5), heights = c(.05, rep((1-.05)/5, 5)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Kaplan & Reed (2018)", cex = 3)
par(mar = c(5, 4, 2, .5), pty = "m")

for(a in unique(kaplan2018Clean$id)[4:28]){
    kaplan2018EXPDQ0 <- coef(kaplan2018EXPDNLME)[which(unique(kaplan2018Clean$id) == a),1]
    kaplan2018EXPDAlpha <- coef(kaplan2018EXPDNLME)[which(unique(kaplan2018Clean$id) == a),2]
    kaplan2018SNDQ0 <- coef(kaplan2018SNDNLME)[which(unique(kaplan2018Clean$id) == a),1]
    kaplan2018SNDAlpha <- coef(kaplan2018SNDNLME)[which(unique(kaplan2018Clean$id) == a),2]
    curve(10^(kaplan2018EXPDQ0) * 10^(1.722925 * (exp(-10^(kaplan2018EXPDAlpha) * 10^(kaplan2018EXPDQ0) * x) - 1)), 
          xlim = c(.001, 1200), ylim = c(0, 15), log = "x", from = .0001, to = 1200, col = "orange", axes = FALSE,
          xlab = "Cost per Drink", ylab = "Purchases", main = paste("ID:",a), lwd = 1.25)
    curve(10^kaplan2018SNDQ0 * exp(- 10^kaplan2018SNDAlpha * 10^kaplan2018SNDQ0 * x), 
          from = .0001, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(x = ifelse(subset(kaplan2018Clean, id == a)$x == 0, .001, 
                      subset(kaplan2018Clean, id == a)$x),
           y = subset(kaplan2018Clean, id == a)$y)
    axis(1, at = tickFunc(), tcl = -.25, labels = NA)
    axis(1, at = 10^(-2:4))
    axis(1, at = .001, 0)
    axis(2, at = c(0, 5, 10, 15, 20, 25))
}
dev.off()



# Scatter Pairs Data ----


svg("correlationpanels.svg", 5,8)

layout(matrix(c(rep(1,2), 2, 3, 
                rep(4, 2), 5, 6, 
                rep(7, 2), 8, 9), ncol = 2, byrow = TRUE), heights = c(1, 6, 1, 6, 1, 6))


par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) Simulations", cex = 2)
par(mar = c(4, 3, 1, .5), pty = "s")
plot(x = coef(koffarnus2015SimsEXPDNLME)[1][,1], y = coef(koffarnus2015SimsSNDNLME)[1][,1], main = expression(italic(Q)[0]),
     xlab = expression(log[10](italic(Q)[0*EXPD])), ylab = expression(log[10](italic(Q)[0*SND])), pch = 22, xlim = c(-.5, 3), ylim = c(-.5, 3))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor(coef(koffarnus2015SimsEXPDNLME)[1][,1], coef(koffarnus2015SimsSNDNLME)[1][,1]),4))),
       bty = "n")
plot(x = coef(koffarnus2015SimsEXPDNLME)[2][,1], y = coef(koffarnus2015SimsSNDNLME)[2][,1], main = expression(alpha),
     xlab = expression(log[10](italic(alpha)[EXPD])), ylab = expression(log[10](italic(alpha)[SND])), xlim = c(-5, 0), ylim = c(-5, 0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor(coef(koffarnus2015SimsEXPDNLME)[1][,1], coef(koffarnus2015SimsSNDNLME)[1][,1]),4))),
       bty = "n")

par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus et al. (2015) CPT", cex = 2, xpd = NA)
par(mar = c(4, 3, 1, .5), pty = "s")
plot(x = coef(koffarnus2015CPTEXPDNLME)[1][,1], y = coef(koffarnus2015CPTSNDNLME)[1][,1], main = expression(italic(Q)[0]),
     xlab = expression(log[10](italic(Q)[0*EXPD])), ylab = expression(log[10](italic(Q)[0*SND])), pch = 22, xlim = c(-.5, 3), ylim = c(-.5, 3))
legend("topleft", 
       # legend = bquote(italic(r)==.(rd(cor(coef(koffarnus2015CPTEXPDNLME)[1][,1], coef(koffarnus2015CPTSNDNLME)[1][,1]),4))),
       # rounds to 1.0000 otherwise
       legend = expression(italic(r)>".9999"),
       bty = "n")
plot(x = coef(koffarnus2015CPTEXPDNLME)[2][,1], y = coef(koffarnus2015CPTSNDNLME)[2][,1], main = expression(alpha),
     xlab = expression(log[10](italic(alpha)[EXPD])), ylab = expression(log[10](italic(alpha)[SND])), xlim = c(-5, 0), ylim = c(-5, 0))
legend("topleft", 
       # legend = bquote(italic(r)==.(rd(cor(coef(koffarnus2015CPTEXPDNLME)[1][,1], coef(koffarnus2015CPTSNDNLME)[1][,1]),4, max = 4))),
       # rounds to 1.0000 otherwise
       legend = expression(italic(r)>".9999"),
       bty = "n")

par(mar = rep(0,4))
plot.new()
text(.5, .5, "Kaplan & Reed (2018) APT", cex = 2, xpd = NA)
par(mar = c(4, 3, 1, .5), pty = "s")
plot(x = coef(kaplan2018EXPDNLME)[1][,1], y = coef(kaplan2018SNDNLME)[1][,1], main = expression(italic(Q)[0]),
     xlab = expression(log[10](italic(Q)[0*EXPD])), ylab = expression(log[10](italic(Q)[0*SND])), pch = 22, xlim = c(-.5, 3), ylim = c(-.5, 3))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor(coef(kaplan2018EXPDNLME)[1][,1], coef(kaplan2018SNDNLME)[1][,1]),4))),
       bty = "n")
plot(x = coef(kaplan2018EXPDNLME)[2][,1], y = coef(kaplan2018SNDNLME)[2][,1], main = expression(alpha),
     xlab = expression(log[10](italic(alpha)[EXPD])), ylab = expression(log[10](italic(alpha)[SND])), xlim = c(-5, 0), ylim = c(-5, 0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor(coef(kaplan2018EXPDNLME)[1][,1], coef(kaplan2018SNDNLME)[1][,1]),4))),
       bty = "n")



dev.off()
# Koffarnus 2015 Simulation Pairs----
svg("koffarnus2015Simsvariables.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = coef(koffarnus2015SimsEXPDNLME)[1][,1],
    SNDQ0 = coef(koffarnus2015SimsSNDNLME)[1][,1],
    AlphaEXPD = coef(koffarnus2015SimsEXPDNLME)[2][,1],
    SNDAlpha = coef(koffarnus2015SimsSNDNLME)[2][,1]), labels = c(expression(log[10](italic(Q)[0*EXPD])), expression(log[10](italic(Q)[0*SND])), 
                                                                  expression(log[10](alpha[EXPD])), expression(log[10](alpha[SND]))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Koffarnus et al. (2015) Simulations SND and EXPD Correlations", scale =FALSE, pch = ".")
dev.off()

# Koffarnus 2015 CPT Pairs----
svg("koffarnus2015CPTvariables.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = coef(koffarnus2015CPTEXPDNLME)[1][,1],
    SNDQ0 = coef(koffarnus2015CPTSNDNLME)[1][,1],
    AlphaEXPD = coef(koffarnus2015CPTEXPDNLME)[2][,1],
    SNDAlpha = coef(koffarnus2015CPTSNDNLME)[2][,1]), labels = c(expression(log[10](italic(Q)[0*EXPD])), expression(log[10](italic(Q)[0*SND])), 
                                                                 expression(log[10](alpha[EXPD])), expression(log[10](alpha[SND]))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Koffarnus et al. (2015) CPT SND and EXPD Correlations", scale =FALSE, pch = ".")
dev.off()

# Kaplan 2018 pairs----
svg("kaplan2018variables.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = coef(kaplan2018EXPDNLME)[1][,1],
    SNDQ0 = coef(kaplan2018SNDNLME)[1][,1],
    AlphaEXPD = coef(kaplan2018EXPDNLME)[2][,1],
    SNDAlpha = coef(kaplan2018SNDNLME)[2][,1],
    apply(kaplan2018Data[c("binges", "totdrinks", "tothours")],
          2, sqrt)), labels = c(expression(log[10](italic(Q)[0*EXPD])), expression(log[10](italic(Q)[0*SND])), 
                                expression(log[10](alpha[EXPD])), expression(log[10](alpha[SND])),
                                expression(sqrt(Binges)), expression(sqrt(Drinks)), expression(sqrt(Hours))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Kaplan & Reed (2018) SND and EXPD Correlations", scale =FALSE, pch = ".")
dev.off()




# EV and Pmax Conversions and Comparisons ----
koff2015SimsEV <- data.frame(EV = 1/(10^coef(koffarnus2015SimsSNDNLME)[2][,1]), 
           CONEV = 1/(10^coef(koffarnus2015SimsEXPDNLME)[2][,1] * EXPDToSND(unique(koffarnus2015SimsClean$k))))
sapply((koff2015SimsEV), median) |> round(2)
sapply((koff2015SimsEV), quantile, c(.25, .75)) |> round(2)
t.test(log10(koff2015SimsEV)[,1], log10(koff2015SimsEV)[,2])$p.value |> round(4)
cor(log10(koff2015SimsEV)[,1], log10(koff2015SimsEV)[,2]) |> round(4)

koff2015CPTEV <- data.frame(SNDEV = 1/(10^coef(koffarnus2015CPTSNDNLME)[2][,1]), 
           CONEV = 1/(10^coef(koffarnus2015CPTEXPDNLME)[2][,1] * EXPDToSND(unique(koffarnus2015CPTClean$k))))
sapply((koff2015CPTEV), median) |> round(2)
sapply((koff2015CPTEV), quantile, c(.25, .75)) |> round(2)
t.test(log10(koff2015CPTEV)[,1], log10(koff2015CPTEV)[,2])$p.value |> round(4)
cor(log10(koff2015CPTEV)[,1], log10(koff2015CPTEV)[,2]) |> round(4)


kaplan2018APTEV <- data.frame(SNDEV = 1/(10^coef(kaplan2018SNDNLME)[2][,1]), 
           CONEV = 1/(10^coef(kaplan2018EXPDNLME)[2][,1] * EXPDToSND(unique(kaplan2018Clean$k))))
sapply((kaplan2018APTEV), median) |> round(2)
sapply((kaplan2018APTEV), quantile, c(.25, .75)) |> round(2)
t.test(log10(kaplan2018APTEV)[,1], log10(kaplan2018APTEV)[,2])$p.value |> round(4)
cor(log10(kaplan2018APTEV)[,1], log10(kaplan2018APTEV)[,2]) |> round(4)


koff2015SimsPmax <- data.frame(SNDPmax = 1/(10^coef(koffarnus2015SimsSNDNLME)[2][,1] * 10^coef(koffarnus2015SimsSNDNLME)[1][,1]), 
           CONPmax = 1/(10^coef(koffarnus2015SimsEXPDNLME)[2][,1] * 10^coef(koffarnus2015SimsEXPDNLME)[1][,1] *
                            EXPDToSND(unique(koffarnus2015SimsClean$k))))
sapply((koff2015SimsPmax), median) |> round(2)
sapply((koff2015SimsPmax), quantile, c(.25, .75)) |> round(2)
t.test(log10(koff2015SimsPmax)[,1], log10(koff2015SimsPmax)[,2])$p.value |> round(4)
cor(log10(koff2015SimsPmax)[,1], log10(koff2015SimsPmax)[,2]) |> round(4)

koff2015CPTPmax <- data.frame(SNDPmax = 1/(10^coef(koffarnus2015CPTSNDNLME)[2][,1] * 10^coef(koffarnus2015CPTSNDNLME)[1][,1]), 
           CONPmax = 1/(10^coef(koffarnus2015CPTEXPDNLME)[2][,1] * 10^coef(koffarnus2015CPTEXPDNLME)[1][,1] *
                            EXPDToSND(unique(koffarnus2015CPTClean$k))))
sapply((koff2015CPTPmax), median) |> round(2)
sapply((koff2015CPTPmax), quantile, c(.25, .75)) |> round(2)
t.test(log10(koff2015CPTPmax)[,1], log10(koff2015CPTPmax)[,2])$p.value |> round(4)
cor(log10(koff2015CPTPmax)[,1], log10(koff2015CPTPmax)[,2]) |> round(4)

kaplan2018APTPmax <- data.frame(SNDPmax = 1/(10^coef(kaplan2018SNDNLME)[2][,1] * 10^coef(kaplan2018SNDNLME)[1][,1]), 
           CONPmax = 1/(10^coef(kaplan2018EXPDNLME)[2][,1] * 10^coef(kaplan2018EXPDNLME)[1][,1] *
                            EXPDToSND(unique(kaplan2018Clean$k))))
sapply((kaplan2018APTPmax), median) |> round(2)
sapply((kaplan2018APTPmax), quantile, c(.25, .75)) |> round(2)
t.test(log10(kaplan2018APTPmax)[,1], log10(kaplan2018APTPmax)[,2])$p.value |> round(4)
cor(log10(kaplan2018APTPmax)[,1], log10(kaplan2018APTPmax)[,2]) |> round(4)

