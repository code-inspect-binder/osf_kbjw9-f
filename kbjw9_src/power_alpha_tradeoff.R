### All mighty power? Clarifying the relationship between power and false-positives ###
# authors: František Bartoš & Maximilian Maier
# contact: f.bartos96@gmail.com, maximilianmaier0401@gmail.com
# 
# supplementary material for recreating figures in the article and/or trying different settings
# the settings can be changed at the beggining of the ploting function
#
# each plotting function starts with specification of:
# pH0:     proportion of the null hypotheses being true
# pH1:     proportion of the alternative hypotheses being true
# effsize: effect size of true alternative hypotheses, in Cohen's d
# n:       overall sample size (assuming that participants are split equally)
# alpha:   significance level alpha (in same cases calculated to correspond to the power)
# power:   power of the statistical tests (in same cases calculated to correspond to the alpha)

library(nleqslv)
library(MCMCpack)
# helper functions for calculationg FPR, alpha, power and nice rounding
get_FPR     <- function(pH1, alpha, power){
  pH0  <- 1 - pH1
  
  FPR  <- pH0*alpha / (pH1*power + pH0*alpha)
  return(FPR)
}
get_alpha   <- function(pH1, power, FPR){
  pH0  <- 1 - pH1
  
  alpha <- - (FPR*pH1*power / (FPR*pH0 - pH0))
  return(alpha)
}
get_power   <- function(alpha, effsize, n){
  
  se      <- 2/sqrt(n)
  nct     <- effsize/se
  
  1 - pnorm(qnorm(1 - alpha/2), nct) + pnorm(qnorm(alpha/2), nct)
}
solve_alpha <- function(x, effsize, n, power){
  y = numeric(1)
  y = get_power(x, effsize, n) - power
  y
}
.r2d        <- function(x){
  x <- format(round(x, 2), nsmall = 2)
  x <- sapply(x, function(i){
    if(i == "1.00"){
      return("1.00")
    }else{
      substr(i, 2, 4)
    }
  })
}
.r3d        <- function(x){
  x <- format(round(x, 3), nsmall = 3)
  x <- sapply(x, function(i){
    if(i == "1.00"){
      return("1.00")
    }else{
      substr(i, 2, 5)
    }
  })
}
# erf and inv_erf: https://stat.ethz.ch/pipermail/r-help/2006-June/108153.html
erf     <- function(x) 2 * pnorm(x * sqrt(2)) - 1
inv_erf <- function(x) qnorm((x + 1)/2)/sqrt(2)
# erf.prime:     https://mathworld.wolfram.com/Erf.html
erf.prime <- function(x) 2/sqrt(pi)*exp(-x^2) 
# inv_erf.prime: https://mathworld.wolfram.com/InverseErf.html
inv_erf.prime  <- function(x) (1/2) * sqrt(pi) * exp( inv_erf(x)^2 )
FDR_grad.alpha <- function(pH1, alpha, ncp, two.sided = TRUE){
  pH0   <- 1 - pH1
  
  if(two.sided){
    beta1 <-     1/2*(1+erf( (sqrt(2) * inv_erf(alpha-1) - ncp)/sqrt(2) )) # lower tail
    beta2 <- 1 - 1/2*(1+erf( (sqrt(2) * inv_erf(1-alpha) - ncp)/sqrt(2) )) # upper tail
    beta  <- beta1 + beta2
  }else{
    beta  <- 1 - 1/2*(1+erf( (sqrt(2) * inv_erf(1-2*alpha) - ncp)/sqrt(2) )) 
  }
  
  nom   <- pH0 * alpha
  denom <- pH0 * alpha + pH1 * beta
  
  if(two.sided){
    pH1beta1.prime <-  1/2*pH1*erf.prime( (sqrt(2) * inv_erf(alpha-1) - ncp)/sqrt(2) )   *inv_erf.prime(alpha-1)
    pH1beta2.prime <- -1/2*pH1*erf.prime( (sqrt(2) * inv_erf(1-alpha) - ncp)/sqrt(2) )*-1*inv_erf.prime(1-alpha)
    pH1beta.prime  <- pH1beta1.prime + pH1beta2.prime    
  }else{
    pH1beta.prime  <- -1/2*pH1*erf.prime( (sqrt(2) * inv_erf(1-2*alpha) - ncp)/sqrt(2) )*-2*inv_erf.prime(1-2*alpha)
  }
  
  nom.prime   <- pH0
  denom.prime <- - ( denom^-2 * (pH0 + pH1beta.prime) )
  
  grad     <- ( nom.prime * denom^-1 ) + ( nom * denom.prime )
  
  return(grad)
}
fill_gradFDR.alpha   <- function(pH1, alpha, ncp, two.sided = TRUE){
  out <- matrix(NA, nrow = length(alpha), ncol = length(ncp))
  for(i in 1:length(alpha)){
    for(j in 1:length(ncp)){
      out[i, j] <- FDR_grad.alpha(pH1, alpha[i], ncp[j], two.sided)
    }
  }
  return(t(out))
}


pdf("FPR_tradeoff_alpha_scaled.pdf",        width = 7.5,  height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100

alpha   <- seq(0,1,.01)
powers  <- get_power(alpha, effsize, n)

par(mar = c(6, 6, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,pH0))
axis(2, seq(0,pH0,.1), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext("Proportion of false-positives", side = 2, line = 3.5, cex = 1.5)
lines(alpha, get_FPR(pH1, alpha, powers), lwd = 3)
dev.off()


pdf("FPR_tradeoff_power_scaled.pdf",        width = 7.5,  height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100

power   <- seq(0,1,.01)
alphas  <- sapply(power, function(p){
  nleqslv::nleqslv(.1, solve_alpha, effsize = effsize, n = n, power = p, control = list(
    stepmax = .00001, maxit = 100000
  ))$x
})

par(mar = c(6, 6, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,pH0))
axis(2, seq(0,pH0,.1), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(alphas[round(power,2) %in% round(seq(0,1,.1),2)]), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext("Proportion of false-positives", side = 2, line = 3.5, cex = 1.5)
lines(power[-1], get_FPR(pH1, alphas[-1], power[-1]), lwd = 3)
dev.off()


pdf("FPR_across_p(H)_and_effect_size.pdf",  width = 12,   height = 9) #####

n        <- 100
pH1s     <- c(.10, .30, .50, .80)
effsizes <- c(.10, .30, .50, .80)

ly <- rbind(c(1, rep(2, 5)),
            cbind(rep(3, 5), rbind(4:8,
                                   cbind(9:12, matrix(13:28, nrow = 4, byrow = T)))),
            29:34)
{
layout(ly, widths = c(.25, .75, rep(1,4)), heights = c(.5, .5, rep(1,4), .75))
par(mar = c(0, 0, 0, 0))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))

plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
text(.5, .5, "Cohen's d", cex = 3)

par(mar = c(0, 0, 0, 3))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))

plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
for(effsize in effsizes){
  plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
  text(.5, .5, effsize, cex = 2)
}


for(pH1 in pH1s){
  plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
  text(.5, .5, pH1, cex = 2)
}

par(mar = c(1, 1, 1, 1))

for(i.p in 1:length(pH1s)){
  for(i.e in 1:length(effsizes)){
    
    effsize <- effsizes[i.e]
    pH1     <- pH1s[i.p]
    
    pH0     <- 1 - pH1
    alpha   <- seq(0,1,.01)
    n       <- 100
    
    powers  <- get_power(alpha, effsize, n)
    
    plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
    if(i.e == 1){
      axis(2, seq(0,1,.25), las = 1)
      mtext("FPR", side = 2, line = 3.5)
    }else{
      axis(2, seq(0,1,.25), labels = rep("",length(seq(0,1,.25))), las = 1)
    }
    if(i.p == length(effsizes)){
      axis(1, seq(0,1,.1), las = 1)
      axis(1, at = seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3)
      if(i.e == 1){
      mtext("Alpha", side = 1, line = 1, at = -.5)
      mtext("Power", side = 1, line = 4, at = -.5)
      }
    }else{
      axis(1, seq(0,1,.1), labels = rep("",length(seq(0,1,.1))), las = 1)
    }
    lines(alpha, get_FPR(pH1, alpha, powers), lwd = 2)
  }
}
mtext(expression(P(H[1])), outer = T, side = 2, adj = .5, line = -3.5, cex = 2, las = 2)
}
dev.off()


pdf("FPR_sample_size_trade_off.pdf",        width = 7,    height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- c(100, 150, 200, 300, 500)

alpha_start <- .05
power_start <- get_power(alpha_start, effsize, n[1])
power_start
pH0 <- 1 - pH1

par(mar = c(4.5, 6.5, 1, 1))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,.1))
axis(2, seq(0,.1,.02), las = 1, cex.axis = 1)
axis(1, seq(0,1,length.out = length(n)), labels = n,las = 1, cex.axis = 1)
mtext("Sample size (n)", side = 1, line = 3, cex = 1.25)
mtext("Proportion of false positives", side = 2, line = 4.5, cex = 1.25)

points(0, get_FPR(pH1, alpha_start, power_start), pch = 16, cex = 2)
lines(seq(0,1,length.out = length(n)), get_FPR(pH1, alpha_start, get_power(alpha_start, effsize, n)),
      type = "b", cex = 1.5)
change_alpha <- sapply(n, function(i){
  nleqslv::nleqslv(.1, solve_alpha, effsize = effsize, n = i, power = power_start, control = list(
    stepmax = .00001, maxit = 10000
  ))$x
  })
lines(seq(0,1,length.out = length(n)), get_FPR(pH1, change_alpha, power_start),
      type = "b", cex = 1.5, pch = 2)
text("fixed alpha", x =.8, y =  .055, cex = 1)
text("fixed power", x =.8, y =  .008, cex = 1)
dev.off()


# binomial charts
.r3d        <- function(x){
  x <- format(round(x, 3), nsmall = 3)
  substr(x, 2, 5)
}
pdf("tasting_tea.pdf",                      width = 9,    height = 4.5) #####
successes <- 0:4
dist_null <- sapply(successes, function(x)dnoncenhypergeom(x, n1 = 4, n2 = 4, m1 = 4, psi = 1))

par(mar = c(4,4.5,1,1), mfrow = c(1,2))
barplot(dist_null, names.arg = successes, las = 1, width = .9, space = .1, ylim = c(0, .60),
        col = c(rep("grey80", 4), rep("grey50", 1)), xlab = "Successes", ylab = "P(Successes)",cex.lab = 1.25)
for(i in 1:5)text(i-.5, dist_null[i] + .015, .r3d(dist_null[i]))

dist_alt <- sapply(successes, function(x)dnoncenhypergeom(x, n1 = 4, n2 = 4, m1 = 4, psi = (.7/(1-.7)) / (.3/(1-.3))))

barplot(dist_alt, names.arg = successes, las = 1, width = .9, space = .1, ylim = c(0, .60),
        col = c(rep("grey80", 4), rep("grey50", 1)), xlab = "Successes", ylab = "P(Successes)",cex.lab = 1.25)
for(i in 1:5)text(i-.5, dist_alt[i] + .015, .r3d(dist_alt[i]))
dev.off()

# power and non-centrality parameter
pdf("power_alpha.pdf",                      width = 9,    height = 4.5) #####
d  <- .5
n  <- 100
mu <- d*sqrt(n)/2

x <- seq(-6, 6, .01)
d_null <- dnorm(x)
d_alt  <- dnorm(x, mu)

par(mar = c(3.5,4,1,1))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(-6, 6), ylim = c(0, .53))
axis(1, seq(-6,6,2))
axis(2, seq(0,.5,.1), las = 2)
mtext("z-statistics", 1, 2.5, cex = 1.25)
mtext("Density", 2, 3, cex = 1.25)

lines(x, d_null, lwd = 2, lty = 2)
lines(x, d_alt,  lwd = 2)
a <- qt(.975, n - 2)
text(a+.1,  .50, expression(~{Phi^-1}(1-frac(alpha,2))), adj = 0)
text(-a-.1, .50, expression(~{Phi^-1}(frac(alpha,2))),   adj = 1)

text(4.5 , .44, expression(1-~Phi[~mu](~{Phi^-1}(1-frac(alpha,2)))),   adj = .5)
text(-4.4, .44, expression(~Phi[~mu](~{Phi^-1}(frac(alpha,2)))),       adj = .5)
arrows(x = a,  x1 = 6.2, y = .40, y1 = .40,   length = 0.10)
arrows(x = -a, x1 =-6.2, y = .40, y1 = .40, length = 0.10)
abline(v =  a)
abline(v = -a)

polygon(c(a,a,x[x>a]), c(0, dnorm(a, mu), d_alt[x > a]), col = rgb(.5,.5,.5,0.5))

dev.off()



pdf("FPR_tradeoff.pdf",                     width = 15,   height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100

alpha   <- seq(0,1,.01)
powers  <- get_power(alpha, effsize, n)

par(mar = c(6, 6, 1, 2), mfrow = c(1,2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,pH0))
axis(2, seq(0,pH0,.1), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext("Proportion of false positives", side = 2, line = 3.5, cex = 1.5)
lines(alpha, get_FPR(pH1, alpha, powers), lwd = 3)


pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100

power   <- seq(0,1,.01)
alphas  <- sapply(power, function(p){
  nleqslv::nleqslv(.1, solve_alpha, effsize = effsize, n = n, power = p, control = list(
    stepmax = .00001, maxit = 100000
  ))$x
})

plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,pH0))
axis(2, seq(0,pH0,.1), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(alphas[round(power,2) %in% round(seq(0,1,.1),2)]), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext("Proportion of false positives", side = 2, line = 3.5, cex = 1.5)
lines(power[-1], get_FPR(pH1, alphas[-1], power[-1]), lwd = 3)
dev.off()


pdf("FPRd_tradeoff_alpha_scaled.pdf",       width = 7.5,  height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100
se      <- 2/sqrt(n)
ncp     <- effsize/se

alpha   <- seq(0,1,.01)
powers  <- get_power(alpha, effsize, n)

par(mar = c(6, 7, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1.5))
axis(2, seq(0,1.5,.3), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext(expression(frac(partialdiff~FDR, partialdiff~alpha)), side = 2, line = 3.5, cex = 1.5)
lines(alpha, fill_gradFDR.alpha(pH1, alpha, ncp), lwd = 3)
dev.off()


pdf("FPRd_tradeoff_alpha_scaled.pdf",       width = 7.5,  height = 5) #####

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100
se      <- 2/sqrt(n)
ncp     <- effsize/se

alpha   <- seq(0,1,.01)
powers  <- get_power(alpha, effsize, n)

par(mar = c(6, 7, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1.5))
axis(2, seq(0,1.5,.3), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext(expression(frac(partialdiff~FDR, partialdiff~alpha)), side = 2, line = 3.5, cex = 1.5)
lines(alpha, fill_gradFDR.alpha(pH1, alpha, ncp), lwd = 3)
dev.off()

pdf("FPRd_tradeoff_power_scaled.pdf",       width = 7.5,  height = 5) 

pH0 <- .5
pH1 <- 1 - pH0

effsize <- .5
n       <- 100
se      <- 2/sqrt(n)
ncp     <- effsize/se

power   <- seq(0,1,.01)
alphas  <- sapply(power, function(p){
  nleqslv::nleqslv(.1, solve_alpha, effsize = effsize, n = n, power = p, control = list(
    stepmax = .00001, maxit = 100000
  ))$x
})

par(mar = c(6, 7, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,12))
axis(2, seq(0,12,2), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(alphas[round(power,2) %in% round(seq(0,1,.1),2)]), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext(expression(frac(partialdiff~FDR, partialdiff~alpha)), side = 2, line = 3.5, cex = 1.5)
lines(alpha, fill_gradFDR.alpha(pH1, alphas, ncp), lwd = 3)
dev.off()

power   <- seq(0,1,.01)
alphas  <- sapply(power, function(p){
  nleqslv::nleqslv(.1, solve_alpha, effsize = effsize, n = n, power = p, control = list(
    stepmax = .00001, maxit = 100000
  ))$x
})

par(mar = c(6, 7, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,12))
axis(2, seq(0,12,2), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(alphas[round(power,2) %in% round(seq(0,1,.1),2)]), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext(expression(frac(partialdiff~FDR, partialdiff~alpha)), side = 2, line = 3.5, cex = 1.5)
lines(alpha, fill_gradFDR.alpha(pH1, alphas, ncp), lwd = 3)
dev.off()
