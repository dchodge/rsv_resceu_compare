logit <- function(x, a, b)
{
  log((x-a)/(b-x))
} 

inv.logit <- function(x, a, b)
{
  a + (b-a)*exp(x)/(exp(x)+1)
} 

D.inv.logit <- function(x, a, b)
{
  (b-a)*exp(x)/(exp(x)+1)^2
} 

logit.prior <- function (x, a, b, FUNC, bool, ...)
{
    if (bool == TRUE)
        out = log(FUNC(inv.logit(x, a, b), ...)*D.inv.logit(x, a, b))
    else
        out = FUNC(inv.logit(x, a, b), ...)*D.inv.logit(x, a, b)
    
    out
}

logit.prior.sample <- function(a, b, FUNC, ...)
{
  x <- FUNC(1, ...)
  while ((x < a) | (x > b)){
    x <- FUNC(1, ...)
  }
  
  logit(x, a, b)
}

prior <- function (x, a, b, FUNC, bool, ...)
{
    if (bool == TRUE)
        out = log(FUNC(x, ...))
    else
        out = FUNC(x, ...)
    
    out
}

prior.sample <- function(a, b, FUNC, ...)
{
    x <- FUNC(1, ...)
    while ((x < a) | (x > b)){
        x <- FUNC(1, ...)
    }
    
    x
}

fit_vals <- function(mean_ci_vals) {
  fitdist_ci <- function(pars, data, dist) {
      qs <- c(0.025, 0.5, 0.975); 1:3 %>% map_dbl(~ (dist(qs[.x], pars[1], pars[2]) - data[.x])^2) %>% sum
  }

  getop_dist <- c(qweibull, qgamma, qlnorm) %>% map(~optim(c(1, 1), fitdist_ci, data = mean_ci_vals, dist = .x)$value) %>% unlist %>% which.min
  dist <- c(qweibull, qgamma, qlnorm)[getop_dist]
  pars <- (c(qweibull, qgamma, qlnorm) %>% map(~optim(c(1, 1), fitdist_ci, data = mean_ci_vals, dist = .x)$par))[getop_dist][[1]]
  list(dist = dist, pars = pars)
}

fit_vals_beta <- function(mean_ci_vals) {
  fitdist_ci_beta <- function(pars, data) {
      qs <- c(0.025, 0.5, 0.975); 1:3 %>% map_dbl(~ (qbeta(qs[.x], pars[1], pars[2]) - data[.x])^2) %>% sum
  }


  getop_dist <- optim(c(2, 2), fitdist_ci_beta, data = mean_ci_vals)$par
  getop_dist
}