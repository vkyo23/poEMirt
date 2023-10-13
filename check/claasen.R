#' To do: Estimating latent democracy support with pgIRT (replication of Claasen (2019))
#' Date: 2023-09-14
#' Author: Ukyo Kanetaka

# 0. Setup -----
## Initialize workspace
rm(list = ls())
set.seed(1)

## Load packages and function
library(tidyverse)
library(tidylog)
library(countrycode)
library(patchwork)
Rcpp::sourceCpp('src/poEMirtdynamic_fit.cpp')
Rcpp::sourceCpp('src/poEMirtbase_fit.cpp')
Rcpp::sourceCpp('src_old/dyn_EMstep.cpp')
source('R/poEMirt.R')
source('R/read_poEMirt.R')

# 1. Data processing -----
## Load data
df <- readxl::read_xlsx('check/data/dem_support.xlsx')

## Add country and item_year ID
df <- df |> 
  filter(Year <= 2020) |> 
  mutate(
    item_year = str_c(Year, Item, sep = '_'),
    item_year_id = factor(item_year) |> 
      as.integer()
  ) |> 
  arrange(item_year) 

## Drop country with less than 2 years
ob <- df |> 
  distinct(Country, Year) |> 
  group_by(Country) |> 
  summarise(obs = n()) |> 
  filter(obs > 1) %>%
  .$Country
df <- df |> 
  filter(Country %in% ob) |> 
  mutate(
    cid = factor(Country) |> 
      as.integer()
  )

## Reference data
### Country
cd <- df |> 
  distinct(cid, country = Country) |> 
  arrange(cid)

### Item Year
it <- df |> 
  distinct(item_year_id, item_year, year = Year, item = Item) 

# Old solution -----
## 1.1. Data for pgIRT ----
I <- nrow(cd)
J <- nrow(it)
T <- length(unique(df$Year))
### Y
Y <- array(NA, dim = c(I, J, 2))
df <- df |> 
  mutate(y0 = Sample - Response)
Y[, , 1] <- df |> 
  distinct(cid, item_year_id, Response) |> 
  spread(key = item_year_id, value = Response) %>%
  .[, -1] |> 
  as.matrix() 
Y[, , 2] <- df |> 
  distinct(cid, item_year_id, y0) |> 
  spread(key = item_year_id, value = y0) %>%
  .[, -1] |> 
  as.matrix() 
rownames(Y) <- cd$country
colnames(Y) <- it$item_year

### N
N <- apply(Y, c(1, 2), sum)

### Item timemap
it <- it |> 
  mutate(
    year_id = factor(year) |> 
      as.integer()
  )
item_timemap <- it$year_id - 1
item_timemap

### Timemap
timemap <- matrix(1, I, T)
timemap2 <- df |> 
  distinct(cid, Year) |> 
  mutate(dum = 1) |> 
  spread(key = Year, value = dum) |> 
  as.matrix() %>%
  .[, -1]
timemap2[is.na(timemap2)] <- 0

### Identifying the first year
for (i in 1:I) {
  first <- min(which(timemap2[i, ] == 1)) - 1
  if (first != 0) {
    timemap[i, 1:first] <- 0
  }
}

### Categories
maxcat <- rep(2, J)

### Matched item
item_match <- rep(NA, nrow(it))

# 2. Initial values -----
## theta
theta_init <- rowMeans(Y[, , 1] / N, na.rm = TRUE) |> 
  scale()
theta_init <- matrix(theta_init, I, T)
constraint <- rep(which(cd$country == 'Denmark'), I) - 1

## alpha and beta
alpha_init <- beta_init <- matrix(NA, J, 1)
for (j in 1:J) {
  cc <- glm(cbind(Y[, j, 1], Y[, j, 2]) ~ theta_init[, 1], family = 'binomial') |> 
    coef()
  alpha_init[j, ] <- cc[1]
  beta_init[j, ] <- cc[2]
}
alpha_init[is.na(alpha_init)] <- beta_init[is.na(beta_init)] <- 0.1

# 3. Run -----
maxit <- 500
tol <- 1e-6
verbose <- 10
std <- FALSE

## 3.1. Fit -----
fit_old <- dyn_EMstep(
  Y = Y,
  N = N,
  lambda = alpha_init,
  gamma = beta_init,
  theta = theta_init,
  maxcat = maxcat,
  constraint = constraint,
  l0 = 0,
  L0 = 25,
  g0 = 1,
  G0 = 0.04,
  t0 = rep(0, I),
  D0 = rep(1, I),
  Delta = 0.0003,
  ind_timemap = timemap,
  item_timemap = item_timemap,
  item_match = item_match,
  maxit = maxit,
  tol = tol,
  verbose = verbose,
  std = std
)
fit_old$theta <- scale(fit_old$theta)
fit_old$theta[109, ]

# 2. poEMirt -----
df2 <- df %>% 
  mutate(y0 = Sample - Response,
         year_id = factor(Year) %>% 
           as.integer()) %>% 
  select(cid, item_year_id, year_id, `1` = Response, `2` = y0) %>% 
  pivot_longer(
    cols = -(cid:year_id),
    names_to = 'k',
    values_to = 'resp'
  ) %>% 
  mutate(
    k = as.integer(k)
  ) %>% 
  arrange(cid)

## Read data
dat <- read_poEMirt(
  data = df2,
  response = 'resp',
  i = 'cid',
  j = 'item_year_id',
  k = 'k',
  t = 'year_id'
)

## theta
fit <- poEMirt(
  data = dat,
  model = 'dynamic',
  priors = list(
    a0 = 0,
    A0 = 25,
    b0 = 0,
    B0 = 25,
    theta0 = 0,
    Delta0 = 1,
    Delta = 0.00005
  )
)
fit$info$init$theta[109, ]
fit$parameter$theta <- scale(fit$parameter$theta)
fit$parameter$theta[109, ]

plot(fit$parameter$theta, fit_old$theta)

# 4. Result -----
## theta
theta_est <- fit$parameter$theta |> # MEAN
  scale() %>% 
  as_tibble() |> 
  mutate(cid = row_number()) |> 
  pivot_longer(
    cols = -cid,
    names_to = 'year_id',
    values_to = 'theta'
  ) |> 
  mutate(
    year_id = year_id |> 
      str_extract('\\d+') |> 
      as.integer()
  ) 

## Merge
theta_est <- theta_est |> 
  left_join(
    cd,
    by = 'cid'
  ) |> 
  left_join(
    it |> 
      distinct(year) %>% 
      arrange(year) %>% 
      mutate(year_id = row_number()),
    by = 'year_id'
  ) |> 
  select(country, year, theta) 

## 4.1. Compare with Claasen's estimates -----
Claasen <- read_csv('check/data/mood_est_v5.csv')

### Add ISO3 code
Claasen <- Claasen |> 
  mutate(
    iso = countrycode(Country, 'country.name', 'iso3n')
  )
theta_est <- theta_est |> 
  mutate(
    iso = countrycode(country, 'country.name', 'iso3n')
  )

### Merge
comp <- theta_est |> 
  left_join(
    Claasen |> 
      select(iso, year = Year, SupDem),
    by = c('iso', 'year')
  )

### Visualization
corr <- comp |> 
  with(
    cor(theta, SupDem, use = 'complete.obs')
  ) |> 
  round(2) 
p1 <- comp |> 
  drop_na() |> 
  ggplot(aes(x = theta, y = SupDem)) +
  geom_abline(slope = 1, intercept = 0,
              linewidth = 1, alpha = .6) +
  geom_point(color = '#0F8766',
             size = 3, alpha = .6) +
  annotate(geom = 'text', x = 2, y = -1, 
           label = paste0('r = ', corr),
           size = 5) +
  theme_bw() +
  xlab('Estimated latent democracy support with pgIRT') +
  ylab("Claassen (2019)'s estimates") +
  ggtitle('Correlation between pgIRT and Claassen (2019) estimates') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(-4, 4)) + ylim(c(-4, 4))
p1

p2 <- comp |> 
  filter(country %in% c('Sweden', 'Denmark', 'Russia', 'United States of America')) |> 
  ggplot(aes(x = year, y = theta)) + 
  geom_line(aes(y = SupDem, color = 'Claassen (2019)'),
            linewidth = .7) +
  geom_line(aes(color = 'pgIRT'), linewidth = .7) +
  scale_color_manual(values = c('black', '#c71585')) +
  facet_wrap(~ country) +
  theme_bw() +
  xlab('Year') + ylab('pgIRT estimates of latent democracy support') +
  ggtitle('Yearly estimates of latent democracy support') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        plot.title = element_text(hjust = .5))
p2
