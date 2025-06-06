library(plotly) 
get_FDR     <- function(pH1, alpha, power){
  pH0  <- 1 - pH1
  
  FDR  <- pH0*alpha / (pH1*power + pH0*alpha)
  return(FDR)
}
get_alpha   <- function(pH1, power, FDR){
  pH0  <- 1 - pH1
  
  alpha <- - (FDR*pH1*power / (FDR*pH0 - pH0))
  return(alpha)
}
get_power  <- function(alpha, ncp, two.sided = TRUE){
  if(two.sided){
    return(1 - pnorm(qnorm(1 - alpha/2), ncp) + pnorm(qnorm(alpha/2), ncp))
  }else{
    return(1 - pnorm(qnorm(1 - alpha), ncp))
  }
}
solve_alpha <- function(x, ncp, two.sided, power){
  y = numeric(1)
  y = get_power(x, ncp, two.sided) - power
  y
}

fill_power  <- function(alpha, ncp, two.sided = TRUE){
  out <- matrix(NA, nrow = length(alpha), ncol = length(ncp))
  for(i in 1:length(alpha)){
    for(j in 1:length(ncp)){
      out[i, j] <- get_power(alpha[i], ncp[j], two.sided)
    }
  }
  return(t(out))
}
fill_FDR    <- function(pH1, alpha, ncp, two.sided = TRUE){
  out <- matrix(NA, nrow = length(alpha), ncol = length(ncp))
  for(i in 1:length(alpha)){
    for(j in 1:length(ncp)){
      out[i, j] <- get_FDR(pH1, alpha[i], get_power(alpha[i], ncp[j], two.sided))
    }
  }
  return(t(out))
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



# Relationship between alpha, ncp and power for a two-sided z-test

alpha   <- seq(0, 1, length.out = 500)
ncp     <- seq(0, 5, length.out = 500)
power   <- fill_power(alpha, ncp)

plot_ly(x = ~alpha, y = ~ncp, z = ~power) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = c(0, .30, .50, .80, .90, 1)
    ))
  )


# Relationship between alpha, ncp and power for a one-sided z-test

alpha   <- seq( 0, 1, length.out = 500)
ncp     <- seq(-5, 5, length.out = 500)
power.o <- fill_power(alpha, ncp, two.sided = FALSE)

plot_ly(x = ~alpha, y = ~ncp, z = ~power.o) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = c(0, .30, .50, .80, .90, 1)
    ))
  )



# Relationship between alpha (and its implied power), ncp, and FDR a two-sided z-test
# green corresponds to $P(H_0) = .2$, blue to $P(H_0) = .5$, and red to $P(H_0) = .8$
  
alpha <- seq(0, 1, length.out = 500)
ncp   <- seq(0, 5, length.out = 500)
FDR.2 <- fill_FDR(.2, alpha, ncp)
FDR.5 <- fill_FDR(.5, alpha, ncp)
FDR.8 <- fill_FDR(.8, alpha, ncp)


plot_ly(x = ~alpha, y = ~ncp, z = ~FDR.5) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = seq(0, 1, .2),
      title    = list(text = "FDR")
    )
  )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.8,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "red"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "red"
                )
              )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.2,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "green"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "green"
                )
              ))


# Relationship between alpha (and its implied power), ncp, and FDR a one-sided z-test
# green corresponds to $P(H_0) = .2$, blue to $P(H_0) = .5$, and red to $P(H_0) = .8$
  
alpha <- seq( 0, 1, length.out = 500)
ncp   <- seq(-5, 5, length.out = 500)
FDR.2o <- fill_FDR(.2, alpha, ncp, two.sided = FALSE)
FDR.5o <- fill_FDR(.5, alpha, ncp, two.sided = FALSE)
FDR.8o <- fill_FDR(.8, alpha, ncp, two.sided = FALSE)


plot_ly(x = ~alpha, y = ~ncp, z = ~FDR.5o) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = seq(0, 1, .2),
      title    = list(text = "FDR")
    )
  )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.8o,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "red"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "red"
                )
              )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.2o,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "green"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "green"
                )
              ))



# Gradient of FDR with respect to alpha (and its implied power) across ncp for a two-sided z-test
# Note that the gradient is increasing towards $\infty$ from both side at ncp = 0 and alpha 0, but is equal to 0 at that point. The figure does not allow to depict it properly and it seems like as if it peaked at ncp around .5 and was decreasing towards 0.
# Green corresponds to $P(H_0) = .2$, blue to $P(H_0) = .5$, and red to $P(H_0) = .8$.

alpha <- seq(0, 1, length.out = 500)
ncp   <- seq(0, 5, length.out = 500)
FDR.2d  <- fill_gradFDR.alpha(.2, alpha, ncp, TRUE)
FDR.5d  <- fill_gradFDR.alpha(.5, alpha, ncp, TRUE)
FDR.8d  <- fill_gradFDR.alpha(.8, alpha, ncp, TRUE)


plot_ly(x = ~alpha, y = ~ncp, z = ~FDR.5d) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = seq(0, 1, .2),
      title    = list(text = "FDR'")
    )
  )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.8d,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "red"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "red"
                )
              )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.2d,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "green"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "green"
                )
              ))


# Gradient of FDR with respect to alpha (and its implied power) across ncp for a one-sided z-test
# Note that the gradient is increasing towards $\infty$ from both side at ncp = 0 and alpha 0, but is equal to 0 at that point. The figure does not allow to depict it properly and it seems like as if it peaked at ncp around .5 and was decreasing towards 0.
# Green corresponds to $P(H_0) = .2$, blue to $P(H_0) = .5$, and red to $P(H_0) = .8$.

alpha <- seq( 0, 1, length.out = 500)
ncp   <- seq(-5, 5, length.out = 500)
FDR.2od  <- fill_gradFDR.alpha(.2, alpha, ncp, FALSE)
FDR.5od  <- fill_gradFDR.alpha(.5, alpha, ncp, FALSE)
FDR.8od  <- fill_gradFDR.alpha(.8, alpha, ncp, FALSE)


plot_ly(x = ~alpha, y = ~ncp, z = ~FDR.5od) %>%
  add_surface(contours = list(
    y = list(
      show = TRUE,
      project = list(y = TRUE),
      usecolormap = FALSE,
      color = "blue"
    ),
    x = list(
      show = TRUE,
      project = list(x = TRUE),
      usecolormap = FALSE,
      color = "blue"
    )
  ),
  opacity = .4) %>%
  layout(scene = list(
    xaxis = list(
      tickvals = c(0, .05, .1, .2, .5, 1)
    ),
    yaxis = list(
      tickvals = pretty(ncp)
    ),
    zaxis = list(
      tickvals = seq(0, 1, .2),
      title    = list(text = "FDR'")
    )
  )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.8od,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "red"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "red"
                )
              )) %>%
  add_surface(x = ~alpha,
              y = ~ncp,
              z = ~FDR.2od,
              opacity = .4,
              contours = list(
                y = list(
                  show = TRUE,
                  project = list(y = TRUE),
                  color = "green"
                ),
                x = list(
                  show = TRUE,
                  project = list(x = TRUE),
                  color = "green"
                )
              ))

