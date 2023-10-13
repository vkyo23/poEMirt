# 0. Setup -----
## Initialize workspace
rm(list = ls())
set.seed(1)

## Load packages
library(tidyverse)

# 1. Simulation -----
## Define sample size
I <- 200
J <- 400
T <- 40
maxK <- 5

## Timemaps
timemap <- matrix(1, I, T)
item_timemap <- rep(0:(T-1), each = J/T)
item_match <- rep(NA, J)

## Number of choice
if(maxK == 2) {
  K <- rep(maxK, J)
} else {
  K <- sample(2:maxK, J, replace = TRUE)
}

## Item parameters
alpha <- beta <- matrix(NA, J, maxK - 1)
for (j in 1:J) {
  Kj <- K[j]
  baseline_probs <- runif(Kj)
  baseline_probs <- exp(baseline_probs)/sum(exp(baseline_probs))
  baseline_sb <- cumprod(1-baseline_probs)/(1-baseline_probs) * baseline_probs
  baseline_sb <- baseline_sb[1:(Kj-1)] 
  baseline_sb <- qlogis(baseline_sb)
  x1_probs <- runif(Kj)
  x1_probs <- exp(x1_probs)/sum(exp(x1_probs))
  x1_sb <- cumprod(1-x1_probs)/(1-x1_probs) * x1_probs
  x1_sb <- x1_sb[1:(Kj-1)]    
  x1_sb <- qlogis(x1_sb)
  alpha[j, 1:(Kj-1)] <- baseline_sb
  beta[j, 1:(Kj-1)] <- x1_sb - baseline_sb
}

# Latent traits
theta <- matrix(NA, I, T)
theta[, 1] <- rnorm(I)
for (t in 2:T) theta[, t] <- rnorm(I, theta[, t-1], 0.1)
con <- apply(theta, 2, which.max)

## Responses
Y <- array(NA, c(I, J, maxK))
n <- sample(500:1000, size = I * J, replace = TRUE) %>% 
  matrix(I, J)
#n <- matrix(1, I, J)
for (i in 1:I) {
  for (j in 1:J) {
    Kj <- K[j]
    t <- item_timemap[j] + 1
    if (Kj > 2) {
      psi <- alpha[j, 1:(Kj-1)] + beta[j, 1:(Kj-1)] * theta[i, t]
      psi <- plogis(psi)
      p <- psi * c(1, cumprod(1 - psi))[-length(psi)]
      p <- c(p, 1 - sum(p))
    } else {
      psi <- alpha[j, 1] + beta[j, 1] * theta[i, t]
      psi <- plogis(psi)
      p <- c(psi, 1 - psi)
    }
    Y[i, j, 1:Kj] <- rmultinom(1, n[i, j], p)
  }
}

# Check
categories <- apply(colSums(Y, na.rm = TRUE), 1, function(x) which(x != 0), simplify = FALSE)
lapply(
  categories,
  function(x) {
    all(sort(x) == order(sort(x)))
  }
) %>% 
  unlist() %>% 
  sum()

# Converting into dataframe 
d <- rep()
for (k in 1:maxK) {
  if (k == 1) {
    d <- Y[, , k] %>% 
      as_tibble() %>% 
      mutate(i = row_number()) %>% 
      pivot_longer(
        cols = -i,
        names_to = 'j',
        values_to = paste0('y', k)
      ) %>% 
      mutate(
        j = j %>% 
          str_remove('V') %>% 
          as.integer()
      ) %>% 
      bind_rows(d, .)
  } else {
    d <- Y[, , k] %>% 
      as_tibble() %>% 
      mutate(i = row_number()) %>% 
      pivot_longer(
        cols = -i,
        names_to = 'j',
        values_to = paste0('y', k)
      ) %>% 
      select(-i, -j) %>% 
      bind_cols(d, .)
  }
}
d <- d %>% 
  left_join(
    tibble(
      j = 1:J,
      t = item_timemap + 1
    ),
    by = 'j'
  ) %>% 
  relocate(t, .before = y1) %>% 
  arrange(j)

# Load funcitons
rm(list = ls()[!ls() %in% c('d', 'con', 'theta', 'alpha', 'beta')])
source('R/read_poEMirt.R')
source('R/poEMirt.R')
source('R/utils.R')
source('R/poEMirt_uncertainty.R')
source('R/predict.R')
source('R/summary.R')
Rcpp::sourceCpp('src/poEMirtdynamic_fit.cpp')
Rcpp::sourceCpp('src/poEMirtdynamic_gibbs_fit.cpp')
Rcpp::sourceCpp('src/poEMirtbase_fit.cpp')
Rcpp::sourceCpp('src/prediction.cpp')

input <- read_poEMirt(
  data = d,
  responses = names(d)[grep('y', names(d))],
  i = 'i',
  j = 'j',
  t = 't'
)
fit <- poEMirt(
  data = input,
  model = 'dynamic',
  constraint = con,
  control = list(
    verbose = 10,
    std = TRUE,
    compute_ll = TRUE,
    tol = 1e-6
  )
)
plot(fit$log_likelihood)

unc <- poEMirt_uncertainty(
  fit = fit,
  method = 'gibbs',
  seed = 1,
  iter = 500,
  control = list(
    verbose = 50,
    warmup = 100,
    thin = 5,
    PG_approx = TRUE,
    save_item_parameters = TRUE
  )
)
gibbs <- summary(unc, parameter = c('alpha', 'beta', 'theta'))

object <- poEMirt_uncertainty(
  fit = fit,
  method = 'bootstrap',
  seed = 1,
  iter = 100,
  control = list(
    verbose = 10,
    warmup = 10,
    PG_approx = TRUE,
    save_item_parameters = TRUE
  )
)
boot <- summary(object, parameter = 'theta')

gibbs_theta <- gibbs %>% 
  filter(parameter == 'theta')
theta_vec <- theta %>% 
  as_tibble() %>% 
  mutate(i = row_number()) %>% 
  pivot_longer(-i) %>% 
  .$value
gibbs_theta %>% 
  mutate(true = as.vector(theta_vec)) %>% 
  ggplot(aes(x = mean, y = true)) +
  geom_point()
cor(gibbs_theta$median, theta_vec)
cor(gibbs_theta$mean, theta_vec)

boot_theta <- boot %>% 
  filter(parameter == 'theta')
boot_theta %>% 
  mutate(true = theta_vec) %>% 
  ggplot(aes(x = estimate, y = true)) +
  geom_point()
cor(boot_theta$estimate, theta_vec)
