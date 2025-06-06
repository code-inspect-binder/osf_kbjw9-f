#run the code and then use the function plotFDR to plot the FDR as a function of alpha for
#"t-test" ("one.sample"; "two sample"), "anova" ("k = <number of groups>), and "corelation" 
#also specifiy the directionality of the test with alternative = "greater", "two.sided", or "less"
library(pwr)



get_power   <- function(alpha, effsize, n, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4){ #get_power for different alpha levels, two tail independent samples t-test
        if(test == "t-test") {
                power <- pwr.t.test(n, effsize, alpha, type = type, alternative = alternative)
        }
        if(test == "anova") {
                power <- pwr.anova.test(k, n, effsize, alpha)
        }
        if(test == "correlation") {
                power <- pwr.r.test(n, effsize, alpha, alternative = alternative)
        }
        return(power$power)
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

get_FPR <- function(pH1, alpha, power){
        pH0  <- 1 - pH1
        
        FPR  <- pH0*alpha / (pH1*power + pH0*alpha)
        return(FPR)
}

plotFDR <- function(pH0, effsize, n,  test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
{
pH1 <- 1 - pH0

alpha   <- seq(0,1,.01)
powers  <- get_power(alpha, effsize, n,  test = test, type = type, alternative = alternative, k = k)

par(mar = c(6, 6, 1, 2))
plot(NA,type='n',axes=FALSE,ann=FALSE, xlim = c(0,1), ylim = c(0,1))
axis(2, seq(0,1,.1), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(seq(0,1,.1)), las = 1, cex.axis = 1.15)
axis(1, seq(0,1,.1), labels = .r2d(get_power(seq(0,1,.1), effsize, n)), las = 1, line = 3, cex.axis = 1.15)
mtext("Alpha", side = 1, line = 1, at = -.15, cex = 1.5)
mtext("Power", side = 1, line = 4, at = -.15, cex = 1.5)
mtext("Proportion of false-positives", side = 2, line = 3.5, cex = 1.5)
lines(alpha, get_FPR(pH1, alpha, powers), lwd = 3)
}
