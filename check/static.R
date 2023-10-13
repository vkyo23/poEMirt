# 0. Setup -----
## Initialize workspace
rm(list = ls())
set.seed(1)

## Load packages and function
library(tidyverse)

# 1. Simulation -----
## Define sample size
I <- 500
J <- 100
maxK <- 5

## Number of choice
K <- sample(2:maxK, J, replace = TRUE)
#K <- rep(maxK, J)

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
theta <- rnorm(I)
con <- which.max(theta)
theta[con]

## Responses
Y <- array(NA, c(I, J, maxK))
n <- sample(500:100, size = I * J, replace = TRUE) %>% 
  matrix(I, J)
#n <- matrix(1, I, J)
for (i in 1:I) {
  for (j in 1:J) {
    Kj <- K[j]
    if (Kj > 2) {
      psi <- alpha[j, 1:(Kj-1)] + beta[j, 1:(Kj-1)] * theta[i]
      psi <- plogis(psi)
      p <- psi * c(1, cumprod(1 - psi))[-length(psi)]
      p <- c(p, 1 - sum(p))
    } else {
      psi <- alpha[j, 1] + beta[j, 1] * theta[i]
      psi <- plogis(psi)
      p <- c(psi, 1 - psi)
    }
    Y[i, j, 1:Kj] <- rmultinom(1, n[i, j], p)
  }
}

categories <- apply(colSums(Y, na.rm = TRUE), 1, function(x) which(x != 0), simplify = FALSE)
lapply(
  categories,
  function(x) {
    all(sort(x) == order(sort(x)))
  }
) %>% 
  unlist()

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

# Load funcitons
rm(list = ls()[!ls() %in% c('d', 'con', 'theta', 'alpha', 'beta')])
source('R/read_poEMirt.R')
source('R/poEMirt.R')
source('R/utils.R')
source('R/poEMirt_uncertainty.R')
source('R/predict.R')
source('R/summary.R')
Rcpp::sourceCpp('src/poEMirtbase_fit.cpp')
Rcpp::sourceCpp('src/poEMirtbase_gibbs_fit.cpp')
Rcpp::sourceCpp('src/prediction.cpp')

data <- read_poEMirt(
  data = d,
  responses = names(d)[grep('y', names(d))],
  i = 'i',
  j = 'j'
)
fit <- poEMirt(
  data = data,
  model = 'static',
  constraint = con,
  control = list(
    verbose = 10,
    std = TRUE,
    compute_ll = TRUE
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
gibbs <- summary(unc)

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
boot <- summary(object)

gibbs_theta <- gibbs %>% 
  filter(parameter == 'theta')
gibbs_theta %>% 
  mutate(true = theta) %>% 
  ggplot(aes(x = mean, y = true)) +
  geom_point()

boot_theta <- boot %>% 
  filter(parameter == 'theta')
boot_theta %>% 
  mutate(true = theta) %>% 
  ggplot(aes(x = estimate, y = true)) +
  geom_point()
