
<!-- README.md is generated from README.Rmd. Please edit that file -->

# poEMirt

<!-- badges: start -->

[![R-CMD-check](https://github.com/vkyo23/poEMirt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vkyo23/poEMirt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`poEMirt` is an `R` package that implements fast EM item response theory
models for public opinion data analysis. Models can be applied to
binary, categorical, binomial and multinomial outcome data. Users can
also estimate statistical uncertainty by parametric bootstrap and Gibbs
Sampling.

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
#> * poEMirt (2023): A fast item response theory models for public opinion data analysis
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
#> * NA rate: 31.2%

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
#> * Model converged at iteration 10 : 0.3 sec.

# Summarize the result
summary(fit, parameter = "theta")
#> * Summarizing following parameters: theta.
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
#> * Summarizing following parameters: theta.
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

#### Gibbs Sampling

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
#> * Summarizing following parameters: theta.
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
