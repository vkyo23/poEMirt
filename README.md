
<!-- README.md is generated from README.Rmd. Please edit that file -->

# poEMirt

<!-- badges: start -->
<!-- badges: end -->

`poEMirt` is an `R` pacakge that implements fast EM item response theory
models for public opinion data analysis. Models include binary,
categorical, binomial and multinomial outcomes along with static and
dynamic estimation. Users can also estimate statistical uncertainty by
parametric bootstrap or Gibbs Sampling.

## Installation

You can install the development version of `poEMirt` from
[GitHub](https://github.com) with:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
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
#> * Dimension of the response array: 200 x 400 x 5

# Fit the model
fit <- poEMirt(
  data = data,
  model = "dynamic",
  control = list(
    verbose = 10,
    constrant = 1
  )
)
#> ================
#> poEMirt starts! 
#> ================
#> * Setting priors.....DONE!
#> * Finding best initial values.....DONE!
#> 
#> === Expectation-Maximization ===
#> * Model converged at iteration 10 : 0.2 sec.

# Summarize the result
summary(fit, parameter = "theta")
#> * Summarizing following parameters: theta.
#> # A tibble: 8,000 × 4
#>    parameter index   reference estimate
#>    <chr>     <chr>   <chr>        <dbl>
#>  1 theta     [1, 1]  [i, t]      -0.525
#>  2 theta     [1, 2]  [i, t]      -0.583
#>  3 theta     [1, 3]  [i, t]      -0.413
#>  4 theta     [1, 4]  [i, t]      -0.525
#>  5 theta     [1, 5]  [i, t]      -0.400
#>  6 theta     [1, 6]  [i, t]      -0.340
#>  7 theta     [1, 7]  [i, t]      -0.398
#>  8 theta     [1, 8]  [i, t]      -0.361
#>  9 theta     [1, 9]  [i, t]      -0.205
#> 10 theta     [1, 10] [i, t]      -0.243
#> # … with 7,990 more rows
```

### Estimation of statistical uncertainty

To estimate statistical uncertainty (e.g., standard deviation), you can
use `poEMirt_uncertainty()`. In this function, you can choose
*parametric bootstrap* (“bootstrap”) or *Gibbs Sampling* (“gibbs”).

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

summary(fit_boot, parameter = "theta", ci = 0.95)
#> * Summarizing following parameters: theta.
#> # A tibble: 8,000 × 7
#>    parameter index   reference estimate     sd ci_lwr ci_upr
#>    <chr>     <chr>   <chr>        <dbl>  <dbl>  <dbl>  <dbl>
#>  1 theta     [1, 1]  [i, t]      -0.525 0.0335 -0.590 -0.459
#>  2 theta     [1, 2]  [i, t]      -0.583 0.0303 -0.642 -0.523
#>  3 theta     [1, 3]  [i, t]      -0.413 0.0270 -0.466 -0.360
#>  4 theta     [1, 4]  [i, t]      -0.525 0.0294 -0.582 -0.467
#>  5 theta     [1, 5]  [i, t]      -0.400 0.0348 -0.468 -0.332
#>  6 theta     [1, 6]  [i, t]      -0.340 0.0321 -0.402 -0.277
#>  7 theta     [1, 7]  [i, t]      -0.398 0.0319 -0.461 -0.336
#>  8 theta     [1, 8]  [i, t]      -0.361 0.0322 -0.424 -0.298
#>  9 theta     [1, 9]  [i, t]      -0.205 0.0286 -0.262 -0.149
#> 10 theta     [1, 10] [i, t]      -0.243 0.0349 -0.311 -0.174
#> # … with 7,990 more rows
```

#### Gibbs Sampling

Using EM estimates as starting values of Gibbs Sampler.
`control$PG_approx = TRUE` accelerates the implementation but is not
recommended for small
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
data (e.g.,
![n \< 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3C%2010 "n < 10")).

``` r
fit_gibbs <- poEMirt_uncertainty(
  fit = fit,
  method = "gibbs",
  seed = 1,
  iter = 500,
  control = list(
    verbose = 50,
    PG_approx = TRUE, 
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
#>    parameter index   reference   mean median     sd ci_lwr ci_upr  rhat
#>    <chr>     <chr>   <chr>      <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#>  1 theta     [1, 1]  [i, t]    -0.527 -0.528 0.0629 -0.667 -0.400 0.995
#>  2 theta     [1, 2]  [i, t]    -0.555 -0.557 0.0337 -0.620 -0.488 0.999
#>  3 theta     [1, 3]  [i, t]    -0.458 -0.459 0.0310 -0.521 -0.405 0.994
#>  4 theta     [1, 4]  [i, t]    -0.544 -0.540 0.0391 -0.631 -0.475 0.993
#>  5 theta     [1, 5]  [i, t]    -0.470 -0.470 0.0635 -0.607 -0.361 1.02 
#>  6 theta     [1, 6]  [i, t]    -0.361 -0.367 0.0823 -0.505 -0.208 1.00 
#>  7 theta     [1, 7]  [i, t]    -0.345 -0.347 0.0572 -0.451 -0.250 0.994
#>  8 theta     [1, 8]  [i, t]    -0.307 -0.300 0.0387 -0.381 -0.230 1.01 
#>  9 theta     [1, 9]  [i, t]    -0.195 -0.196 0.0433 -0.272 -0.117 1.02 
#> 10 theta     [1, 10] [i, t]    -0.245 -0.244 0.0528 -0.340 -0.140 1.04 
#> # … with 7,990 more rows
```
