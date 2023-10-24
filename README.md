
<!-- README.md is generated from README.Rmd. Please edit that file -->

# poEMirt

<!-- badges: start -->

[![R-CMD-check](https://github.com/vkyo23/poEMirt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vkyo23/poEMirt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`poEMirt` is an `R` package that implements fast EM item response theory
for public opinion data analysis to estimate cross-national latent
public opinion over years. This package utilizes the Polya-Gamma data
augmentation scheme (Polson, Scott & Windle 2013) and exact EM algorithm
(Goplerud 2019) for fast computation. This package incorporates direct
C++ programming and does not depend on ready-made other statistical
modeling languages, thus the computation is more efficient. It can be
applied to not only binomial and multinomial outcomes but also bernoulli
and categorial ones. Users can estimate statistical uncertainty by
parametric bootstrap and Gibbs Sampling.

## Installation

You can install the development version of `poEMirt` from
[GitHub](https://github.com) with:

``` r
remotes::install_github("vkyo23/poEMirt")
```

## Usage

### Multinomial outcome

``` r
library(poEMirt)
data("sim_data_dynamic")

head(sim_data_dynamic)
#> # A tibble: 6 × 8
#>       i     j     t    y1    y2    y3    y4    y5
#>   <int> <int> <dbl> <int> <int> <int> <int> <int>
#> 1     1     1     1   366   296    NA    NA    NA
#> 2     2     1     1   451   260    NA    NA    NA
#> 3     3     1     1   346   259    NA    NA    NA
#> 4     4     1     1   538   163    NA    NA    NA
#> 5     5     1     1   286   267    NA    NA    NA
#> 6     6     1     1   497   137    NA    NA    NA
tail(sim_data_dynamic)
#> # A tibble: 6 × 8
#>       i     j     t    y1    y2    y3    y4    y5
#>   <int> <int> <dbl> <int> <int> <int> <int> <int>
#> 1   195   400    40   272    95   454    NA    NA
#> 2   196   400    40   162    84   360    NA    NA
#> 3   197   400    40   168   134   445    NA    NA
#> 4   198   400    40   285    97   532    NA    NA
#> 5   199   400    40   284    99   465    NA    NA
#> 6   200   400    40   282   130   539    NA    NA

# Convert into poEMirt-readable data
data <- read_poEMirt(
  data = sim_data_dynamic,
  responses = paste0('y', 1:5),
  i = "i",
  j = "j",
  t = "t"
)
summary(data)
#> ----- poEMirtData summary -----
#> * Data size:
#>   - I =  200 
#>   - J =  400 
#>   - T =  40 
#>   - min(Kj) =  2 
#>   - max(Kj) =  5 
#>   - mean(Kj) = 3.44 
#> * NA rate: 0%

# Fit the model
fit <- poEMirt(
  data = data,
  model = "dynamic",
  control = list(
    verbose = 10
  )
)
#> === poEMirt starts! ===
#> * Setting priors.....DONE!
#> * Finding best initial values.....DONE!
#> 
#> === Expectation-Maximization ===
#> * Model converged at iteration 10 : 0.4 sec.

# Summarize the result
summary(fit, parameter = "theta")
#> # A tibble: 8,000 × 4
#>    parameter index  reference estimate
#>    <chr>     <chr>  <chr>        <dbl>
#>  1 theta     [1,1]  [i,t]       -0.525
#>  2 theta     [1,2]  [i,t]       -0.582
#>  3 theta     [1,3]  [i,t]       -0.413
#>  4 theta     [1,4]  [i,t]       -0.524
#>  5 theta     [1,5]  [i,t]       -0.400
#>  6 theta     [1,6]  [i,t]       -0.339
#>  7 theta     [1,7]  [i,t]       -0.398
#>  8 theta     [1,8]  [i,t]       -0.361
#>  9 theta     [1,9]  [i,t]       -0.205
#> 10 theta     [1,10] [i,t]       -0.242
#> # ℹ 7,990 more rows
```

### Estimation of statistical uncertainty

To estimate statistical uncertainty (e.g., standard deviation), you can
use `poEMirt_uncertainty()`. In this function, you can choose
*parametric bootstrap* (`"bootstrap"`) or *Gibbs sampling* (`"gibbs"`).

#### Bootstrap

``` r
fit_boot <- poEMirt_uncertainty(
  fit = fit,
  method = "bootstrap",
  seed = 1,
  iter = 100,
  control = list(
    verbose = 10
  )
)
#> === Parametric bootstrap to estimate statistical uncertainty for poEMirt ===
#> * Bootstrap 1 / 100 
#> * Bootstrap 10 / 100 
#> * Bootstrap 20 / 100 
#> * Bootstrap 30 / 100 
#> * Bootstrap 40 / 100 
#> * Bootstrap 50 / 100 
#> * Bootstrap 60 / 100 
#> * Bootstrap 70 / 100 
#> * Bootstrap 80 / 100 
#> * Bootstrap 90 / 100 
#> * Bootstrap 100 / 100 
#> DONE!

summary(fit_boot, parameter = "theta", ci = 0.95)
#> # A tibble: 8,000 × 7
#>    parameter index  reference estimate     sd ci_lwr ci_upr
#>    <chr>     <chr>  <chr>        <dbl>  <dbl>  <dbl>  <dbl>
#>  1 theta     [1,1]  [i,t]       -0.525 0.0380 -0.599 -0.450
#>  2 theta     [1,2]  [i,t]       -0.582 0.0319 -0.645 -0.520
#>  3 theta     [1,3]  [i,t]       -0.413 0.0267 -0.465 -0.360
#>  4 theta     [1,4]  [i,t]       -0.524 0.0316 -0.586 -0.463
#>  5 theta     [1,5]  [i,t]       -0.400 0.0340 -0.466 -0.333
#>  6 theta     [1,6]  [i,t]       -0.339 0.0350 -0.408 -0.271
#>  7 theta     [1,7]  [i,t]       -0.398 0.0319 -0.461 -0.336
#>  8 theta     [1,8]  [i,t]       -0.361 0.0341 -0.428 -0.294
#>  9 theta     [1,9]  [i,t]       -0.205 0.0274 -0.259 -0.152
#> 10 theta     [1,10] [i,t]       -0.242 0.0328 -0.307 -0.178
#> # ℹ 7,990 more rows
```

#### Gibbs sampling

Using EM estimates as starting values of Gibbs sampler.

``` r
fit_gibbs <- poEMirt_uncertainty(
  fit = fit,
  method = "gibbs",
  seed = 1,
  iter = 500,
  control = list(
    verbose = 50,
    warmup = 100,
    thin = 5
  )
)
#> === Gibbs Sampling to estimate statistical uncertainty for poEMirt ===
#> * Warmup 1 / 600
#> * Warmup 50 / 600
#> * Warmup 100 / 600
#> * Sampling 150 / 600
#> * Sampling 200 / 600
#> * Sampling 250 / 600
#> * Sampling 300 / 600
#> * Sampling 350 / 600
#> * Sampling 400 / 600
#> * Sampling 450 / 600
#> * Sampling 500 / 600
#> * Sampling 550 / 600
#> * Sampling 600 / 600

summary(fit_gibbs, parameter = "theta", ci = 0.95)
#> # A tibble: 8,000 × 9
#>    parameter index  reference   mean median     sd ci_lwr ci_upr  rhat
#>    <chr>     <chr>  <chr>      <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#>  1 theta     [1,1]  [i,t]     -0.537 -0.541 0.0549 -0.635 -0.433 0.997
#>  2 theta     [1,2]  [i,t]     -0.557 -0.553 0.0400 -0.641 -0.488 1.00 
#>  3 theta     [1,3]  [i,t]     -0.471 -0.469 0.0360 -0.534 -0.395 0.991
#>  4 theta     [1,4]  [i,t]     -0.547 -0.550 0.0394 -0.605 -0.471 0.992
#>  5 theta     [1,5]  [i,t]     -0.463 -0.464 0.0711 -0.594 -0.331 1.00 
#>  6 theta     [1,6]  [i,t]     -0.374 -0.374 0.0695 -0.507 -0.249 1.00 
#>  7 theta     [1,7]  [i,t]     -0.349 -0.349 0.0596 -0.465 -0.245 1.02 
#>  8 theta     [1,8]  [i,t]     -0.313 -0.309 0.0465 -0.402 -0.222 1.00 
#>  9 theta     [1,9]  [i,t]     -0.203 -0.203 0.0445 -0.290 -0.123 1.01 
#> 10 theta     [1,10] [i,t]     -0.244 -0.245 0.0501 -0.344 -0.151 0.997
#> # ℹ 7,990 more rows
```

### Binary, categorical outcome

Notably, `poEMirt` can also be applied to roll-call vote,
individual-level survey, and exam data analysis where outcomes are
binary or categorical since these are the special case of `poEMirt`
where n = 1 for each observation.

#### Binary

Using the roll-call vote data of U.S. Supreme Court `Rehnquist` from
[`{MCMCpack}`](https://github.com/cran/MCMCpack) by Martin, Quinn & Park
(2011). The outcome of interest is 1 if a judge votes yea and 0 if nay.

``` r
library(tidyr)
library(dplyr)
data(Rehnquist, package = "MCMCpack")

df <- Rehnquist |> 
  mutate(rcid = row_number()) |> 
  select(-term) |> 
  mutate(
    time = factor(time) |> 
      as.integer()
  ) |> 
  pivot_longer(
    cols = -c(time, rcid),
    names_to = "judge",
    values_to = "vote"
  ) 
head(df)
#> # A tibble: 6 × 4
#>    time  rcid judge      vote
#>   <int> <int> <chr>     <dbl>
#> 1     1     1 Rehnquist     0
#> 2     1     1 Stevens       1
#> 3     1     1 O.Connor      0
#> 4     1     1 Scalia        0
#> 5     1     1 Kennedy       1
#> 6     1     1 Souter        1

#For binary data, users must create the columns of `outcome = 1` and `outcome = 0` data respectively. 
df <- df |> 
  mutate(
    y1 = 1L * (vote == 1),
    y0 = 1L * (vote == 0)
  )
head(df)
#> # A tibble: 6 × 6
#>    time  rcid judge      vote    y1    y0
#>   <int> <int> <chr>     <dbl> <int> <int>
#> 1     1     1 Rehnquist     0     0     1
#> 2     1     1 Stevens       1     1     0
#> 3     1     1 O.Connor      0     0     1
#> 4     1     1 Scalia        0     0     1
#> 5     1     1 Kennedy       1     1     0
#> 6     1     1 Souter        1     1     0

# Read data
data <- read_poEMirt(
  data = df,
  responses = c("y1", "y0"),
  i = "judge",
  j = "rcid",
  t = "time"
)
summary(data)
#> ----- poEMirtData summary -----
#> * Data size:
#>   - I =  9 
#>   - J =  485 
#>   - T =  11 
#>   - min(Kj) =  2 
#>   - max(Kj) =  2 
#>   - mean(Kj) = 2 
#> * NA rate: 0.5%

# Constraint
con <- which(rownames(data$response) == "Thomas")
con <- rep(con, max(df$time))

# Fit
fit <- poEMirt(
  data = data,
  model = "dynamic",
  constraint = con,
  control = list(
    verbose = 10
  )
)
#> === poEMirt starts! ===
#> * Setting priors.....DONE!
#> * Finding best initial values.....DONE!
#> 
#> === Expectation-Maximization ===
#> * Model converged at iteration 8 : 0 sec.

summary(fit, parameter = "theta")
#> # A tibble: 99 × 4
#>    parameter index       reference estimate
#>    <chr>     <chr>       <chr>        <dbl>
#>  1 theta     [Breyer,1]  [i,t]        -1.43
#>  2 theta     [Breyer,2]  [i,t]        -1.44
#>  3 theta     [Breyer,3]  [i,t]        -1.47
#>  4 theta     [Breyer,4]  [i,t]        -1.47
#>  5 theta     [Breyer,5]  [i,t]        -1.47
#>  6 theta     [Breyer,6]  [i,t]        -1.47
#>  7 theta     [Breyer,7]  [i,t]        -1.50
#>  8 theta     [Breyer,8]  [i,t]        -1.50
#>  9 theta     [Breyer,9]  [i,t]        -1.49
#> 10 theta     [Breyer,10] [i,t]        -1.47
#> # ℹ 89 more rows
```

#### Categorical

Using the roll-call vote data of United Nations General Assembly
(session 1 to 74) from [`{unvotes}`](https://github.com/dgrtwo/unvotes)
by Robinson & Goguen-Compagnoni (2021). Note that the original data is
from Voeten (2013). The outcome of interest consists of “yes”,
“abstain”, and “no”.

``` r
data(un_roll_calls, package = "unvotes")
data(un_votes, package = "unvotes")

# Add sessions
df <- un_votes |> 
  left_join(
    un_roll_calls |> 
      select(rcid, session),
    by = "rcid"
  ) 
head(df)
#> # A tibble: 6 × 5
#>    rcid country            country_code vote  session
#>   <dbl> <chr>              <chr>        <fct>   <dbl>
#> 1     3 United States      US           yes         1
#> 2     3 Canada             CA           no          1
#> 3     3 Cuba               CU           yes         1
#> 4     3 Haiti              HT           yes         1
#> 5     3 Dominican Republic DO           yes         1
#> 6     3 Mexico             MX           yes         1

# Create response columns
df <- df |> 
  mutate(tmp = 1) |> 
  spread(
    key = vote, 
    value = tmp
  ) |> 
  mutate(
    across(
      .cols = yes:no,
      .fns = function(x) ifelse(is.na(x), 0, x)
    )
  )
head(df)
#> # A tibble: 6 × 7
#>    rcid country            country_code session   yes abstain    no
#>   <dbl> <chr>              <chr>          <dbl> <dbl>   <dbl> <dbl>
#> 1     3 United States      US                 1     1       0     0
#> 2     3 Canada             CA                 1     0       0     1
#> 3     3 Cuba               CU                 1     1       0     0
#> 4     3 Haiti              HT                 1     1       0     0
#> 5     3 Dominican Republic DO                 1     1       0     0
#> 6     3 Mexico             MX                 1     1       0     0

# Read data
data <- read_poEMirt(
  data = df,
  responses = c("yes", "abstain", "no"), # "no" is the baseline 
  i = "country",
  j = "rcid",
  t = "session"
)
#> * Remove following items due to no variation in responses
#>   - 11 122 182 368 399 428 440 441 443 444 447 448 449 451 452 460 461 462 484 487 503 505 642 733 799 876 913 943 944 946 1177 1196 1211 1251 1288 1299 1348 1349 1364 1378 1393 1400 1408 1417 1450 1472 1519 1542 1565 1595 1673 1736 1754 1755 1759 1763 1765 1770 1771 1846 1851 1856 1907 1960 2068 2363 2453 3009 3142 3245 3680 4024 4025 4385 4534 4852 4997 5034 5342 5496 5574 5893
summary(data)
#> ----- poEMirtData summary -----
#> * Data size:
#>   - I =  200 
#>   - J =  6120 
#>   - T =  74 
#>   - min(Kj) =  2 
#>   - max(Kj) =  3 
#>   - mean(Kj) = 2.75 
#> * NA rate: 29.7%

# Constraint
con <- which(rownames(data$response) == "United States")
con <- rep(con, max(df$session))

# Fit
fit <- poEMirt(
  data = data,
  model = "dynamic",
  constraint = con,
  control = list(
    verbose = 10
  )
)
#> === poEMirt starts! ===
#> * Setting priors.....DONE!
#> * Finding best initial values.....DONE!
#> 
#> === Expectation-Maximization ===
#> Iteration 10: eval = 0.000136153
#> Iteration 20: eval = 4.45304e-05
#> Iteration 30: eval = 5.49363e-06
#> * Model converged at iteration 40 : 16.4 sec.

summary(fit, parameter = "theta")
#> # A tibble: 10,919 × 4
#>    parameter index            reference estimate
#>    <chr>     <chr>            <chr>        <dbl>
#>  1 theta     [Afghanistan,1]  [i,t]       -0.688
#>  2 theta     [Afghanistan,2]  [i,t]       -0.735
#>  3 theta     [Afghanistan,3]  [i,t]       -0.792
#>  4 theta     [Afghanistan,4]  [i,t]       -0.860
#>  5 theta     [Afghanistan,5]  [i,t]       -0.897
#>  6 theta     [Afghanistan,6]  [i,t]       -0.913
#>  7 theta     [Afghanistan,7]  [i,t]       -0.933
#>  8 theta     [Afghanistan,8]  [i,t]       -0.928
#>  9 theta     [Afghanistan,9]  [i,t]       -0.928
#> 10 theta     [Afghanistan,10] [i,t]       -0.942
#> # ℹ 10,909 more rows
```

# References

- Goplerud, M. (2019). “A Multinomial Framework for Ideal Point
  Estimation”. *Political Analysis*, 27(1), 69-89.
- Martin A.D., Quinn K.M. & Park J.H. (2011). “MCMCpack: Markov Chain
  Monte Carlo in R.” *Journal of Statistical Software*, 42(9), 22.
- Polson, N. G., Scott, J. G., & Windle, J. (2013). “Bayesian inference
  for logistic models using Pólya–Gamma latent variables”. *Journal of
  the American statistical Association*, 108(504), 1339-1349.
- Robinson, D., & Goguen-Compagnoni, N. (2021) “unvotes: United Nations
  General Assembly Voting Data”. R package.
- Voeten, E. (2013). *Data and analyses of voting in the United Nations:
  General Assembly*. Routledge handbook of international organization,
  54-66.
