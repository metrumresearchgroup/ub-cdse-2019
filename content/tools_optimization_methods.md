Estimate parameters in a PBPK model
================
Metrum Research Group

  - [Packages and setup](#packages-and-setup)
  - [Reference](#reference)
  - [Data](#data)
  - [PBPK model: pitavastatin / CsA
    DDI](#pbpk-model-pitavastatin-csa-ddi)
  - [Objective function](#objective-function)
      - [Prediction function](#prediction-function)
  - [Data grooming](#data-grooming)
  - [Optimize](#optimize)
      - [`nloptr::newuoa`: minimization without
        derivatives](#nloptrnewuoa-minimization-without-derivatives)
          - [The final objective function value and
            estimates](#the-final-objective-function-value-and-estimates)
      - [`optim`: Nelder-Mead](#optim-nelder-mead)
      - [`neldermead`: Alternate
        Nelder-Mead](#neldermead-alternate-nelder-mead)
      - [`DEoptim`: differential evolution
        algorithm](#deoptim-differential-evolution-algorithm)
          - [DA for the plot](#da-for-the-plot)
      - [`GenSA`: simulated annealing](#gensa-simulated-annealing)
      - [`hydroPSO`: particle swarm
        optimization](#hydropso-particle-swarm-optimization)
      - [`nloptr::CRS`: controlled random
        search](#nloptrcrs-controlled-random-search)
  - [Compare optimization methods](#compare-optimization-methods)

# Packages and setup

``` r
library(tidyverse)
library(PKPDmisc)
library(mrgsolve)
library(nloptr)
library(DEoptim)
library(GenSA)
library(hydroPSO)
source("script/functions.R")
source("script/global.R")
```

``` r
set.seed(10101)
```

``` r
theme_set(theme_bw())
theme_update(legend.position = "top")
scale_colour_discrete <- function(...) scale_color_brewer(palette="Set2")
```

Models are located here:

``` r
model_dir <- "model"
```

# Reference

**Quantitative Analyses of Hepatic OATP-Mediated Interactions Between
Statins and Inhibitors Using PBPK Modeling With a Parameter Optimization
Method**

  - T Yoshikado, K Yoshida, N Kotani, T Nakada, R Asaumi, K Toshimoto, K
    Maeda, H Kusuhara and Y Sugiyama

  - CLINICAL PHARMACOLOGY & THERAPEUTICS | VOLUME 100 NUMBER 5 |
    NOVEMBER 2016

  - <https://www.ncbi.nlm.nih.gov/pubmed/27170342>

# Data

  - Example taken from figure 4a from the publication
  - Using this as example data to fit

<!-- end list -->

``` r
data.file <- "data/fig4a.csv"

data <-
  data.file %>% 
  read_csv() %>% 
  mutate(
    profile = NULL, 
    type=ID, 
    typef=factor(ID, labels = c("Statin", "Statin+CsA")), 
    DV = ifelse(DV==-1, NA_real_, DV)
  )
```

  - The goal is to fit the pitavastatin data either alone or in
    combination with cyclosporin administered 1 hour before the
    pitavastatin

<!-- end list -->

``` r
ggplot(data=data,aes(time,DV)) + 
  geom_point(aes(col = typef), size = 3) + 
  geom_line(col = "darkgrey", aes(group = typef)) + 
  scale_y_continuous(trans="log", limits=c(0.1,300), breaks=logbr()) 
```

![](tools_optimization_methods_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# PBPK model: pitavastatin / CsA DDI

  - Check out the model / data with a quick simulation

<!-- end list -->

``` r
mod <- mread_cache("yoshikado", model_dir)
```

Make some persistent updates to the model

  - Simulate out to 14 hours
  - Only interested in `CP`, the pitavastatin concentration

<!-- end list -->

``` r
mod <- mod %>% update(end=14, delta=0.1) %>% Req(CP) 
```

A practice simulation

``` r
dose <- filter(data, evid==1) %>% mutate(typef=NULL)

sims <- 
  mod %>% 
  mrgsim_d(dose, obsaug=TRUE) %>% 
  mutate(type = typef(ID))

ggplot(sims, aes(time,CP,col=type)) + 
  geom_line(lwd = 1) + 
  scale_x_continuous(breaks = seq(0,12,2)) + 
  scale_y_log10(name = "Pitavastatin concentration")
```

![](tools_optimization_methods_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
sims %>% 
  group_by(type) %>% 
  summarise(auc = auc_partial(time,CP)) %>% 
  mutate(fold_increase = auc /first(auc))
```

    ## # A tibble: 2 x 3
    ##   type                 auc fold_increase
    ##   <fct>              <dbl>         <dbl>
    ## 1 Pitavastatin alone  44.1          1   
    ## 2 Pitavastatin + CsA 161.           3.65

# Objective function

  - Least squares objective function
  - Weighted by the observations

Arguments:

  - `dv` the observed data
  - `pred` the predicted data

<!-- end list -->

``` r
wss <- function(dv, pred, weight = 1/dv) {
  sum(((dv-pred)*weight)^2,na.rm=TRUE) 
}
```

#### Prediction function

  - Let’s go through step by step what each line is doing for us

Arguments:

  - `p` the parameters proposed by the optimizer
  - `.data` the simulation template (doses and observation records)
  - `yobs` a vector of observed data which matches observations in
    `.data`
  - `pred` logical; if `TRUE`, just return predicted data

<!-- end list -->

``` r
sim_ofv <- function(p, data, pred = FALSE) {
  
  names(p) <- names(theta)
  
  p <- lapply(p,exp)
  
  out <- mod %>% param(p) %>% mrgsim_q(data, output="df")
  
  if(pred) return(out)
  
  ofv <- wss(data[["DV"]], out[["CP"]])
  
  return(ofv)
  
  #return(-1*sum(dnorm(log(yobs),log(out$CP),.par$sigma,log=TRUE)))
  
}
```

What this function does:

1.  Take in arguments; focus is on a new set of parameters `p` proposed
    by the optimizer; other arguments are just fixed data that we need
2.  Get the parameters out of log scale
3.  Also, put names on the list of parameters; this is crutial
4.  Update the model object with the new parameters
5.  (optionally simulate and return)
6.  Simulate from the data set, taking only observed values
7.  Calculate and return the objective function value

# Data grooming

  - Pick out the observations
  - Drop the non-numeric columns

<!-- end list -->

``` r
data <-  dplyr::select(data, -typef)
```

# Optimize

First, set up the initial estimates

``` r
theta <- c(
  fbCLintall = 1.2, 
  ikiu = 1.2, 
  fbile = 0.9, 
  ka = 0.1, 
  ktr = 0.1
) %>% log()
```

## `nloptr::newuoa`: minimization without derivatives

``` r
fit <- nloptr::newuoa(x0 = theta, fn = sim_ofv, data = data)
```

``` r
fit
```

    ## $par
    ## [1] -0.20418269 -4.51446274 -1.06752471 -0.01113143 -0.37149033
    ## 
    ## $value
    ## [1] 0.6860763
    ## 
    ## $iter
    ## [1] 406
    ## 
    ## $convergence
    ## [1] 1
    ## 
    ## $message
    ## [1] "NLOPT_SUCCESS: Generic success return value."

#### The final objective function value and estimates

``` r
sim_ofv(fit$par,data=data)
```

    ## [1] 0.6860763

``` r
exp(fit$par) %>% set_names(names(theta))
```

    ## fbCLintall       ikiu      fbile         ka        ktr 
    ## 0.81531341 0.01094949 0.34385861 0.98893030 0.68970568

## `optim`: Nelder-Mead

``` r
fit1b <- optim(theta, sim_ofv, data=data, control = list(maxit = 1000))
```

## `neldermead`: Alternate Nelder-Mead

``` r
fit1c <- nloptr::neldermead(x0=theta, fn=sim_ofv, data = data )
```

## `DEoptim`: differential evolution algorithm

<https://en.wikipedia.org/wiki/Differential_evolution>

“Performs evolutionary global optimization via the Differential
Evolution algorithm.”

``` r
lower <- rep(-6,length(theta)) %>% setNames(names(theta))
upper <- rep(5, length(theta)) %>% setNames(names(theta))

set.seed(330303)

decontrol <- DEoptim.control(
  trace = 10,
  NP=10*length(theta), 
  CR=0.925, 
  F=0.85,
  itermax=90, 
  storepopfrom=0
)

fit2 <- DEoptim(
  fn=sim_ofv, 
  lower=lower,
  upper=upper, 
  control=decontrol,
  data=data
)
```

    ## Iteration: 10 bestvalit: 2.359438 bestmemit:    0.090828   -4.810051   -1.032300   -0.596372   -0.748306
    ## Iteration: 20 bestvalit: 2.359438 bestmemit:    0.090828   -4.810051   -1.032300   -0.596372   -0.748306
    ## Iteration: 30 bestvalit: 0.745879 bestmemit:   -0.220189   -4.444982   -1.045812   -0.135432   -0.454597
    ## Iteration: 40 bestvalit: 0.733471 bestmemit:   -0.201757   -4.500855   -1.075498   -0.170296   -0.414804
    ## Iteration: 50 bestvalit: 0.691808 bestmemit:   -0.220169   -4.512342   -1.101670   -0.043517   -0.409125
    ## Iteration: 60 bestvalit: 0.688542 bestmemit:   -0.210497   -4.516279   -1.078079    0.011272   -0.380982
    ## Iteration: 70 bestvalit: 0.686352 bestmemit:   -0.203172   -4.509802   -1.066623   -0.004658   -0.362483
    ## Iteration: 80 bestvalit: 0.686121 bestmemit:   -0.203278   -4.514187   -1.064788   -0.010720   -0.370142
    ## Iteration: 90 bestvalit: 0.686081 bestmemit:   -0.204019   -4.514555   -1.067274   -0.011397   -0.372682

#### DA for the plot

``` r
pops <- lapply(fit2$member$storepop, as.data.frame)
hx <- bind_rows(pops)
hx <- mutate(hx, iteration=rep(1:decontrol$itermax,each=decontrol$NP))
hx <- mutate(hx, pop = rep(1:decontrol$NP, time=decontrol$itermax))
hxm <- gather(hx, variable, value, 1:5) %>% mutate(value = exp(value))
best <- as_tibble(fit2$member$bestmemit) %>% 
  mutate(iteration = 1:decontrol$itermax)
bestm <- gather(best,variable,value,1:5) %>% mutate(value = exp(value))
```

``` r
ggplot(data=hxm) + 
  geom_line(aes(iteration,value,group=pop),col="darkslateblue") + 
  geom_line(data=bestm,aes(iteration,value),col="orange",lwd=1) + 
  scale_y_continuous(trans="log", breaks=10^seq(-4,4), name="Parameter value") + 
  facet_wrap(~variable, ncol=2, scales="free_y") + theme_bw()
```

![](tools_optimization_methods_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## `GenSA`: simulated annealing

``` r
set.seed(11001)

sacontrol <- list(maxit = 100, nb.stop.improvement = 20, verbose = TRUE)

fit3 <- GenSA(
  NULL, sim_ofv, lower=lower+1, upper=upper-1, data = data, control = sacontrol
)
```

    ## Initializing par with random data inside bounds
    ## It: 1, obj value: 3.482240586
    ## It: 37, obj value: 0.6860966824

## `hydroPSO`: particle swarm optimization

<https://en.wikipedia.org/wiki/Particle_swarm_optimization>

``` r
set.seed(22022013)

fit4 <- hydroPSO(
  theta, fn = "sim_ofv", lower = lower, upper = upper, 
  control = list(maxit = 100, REPORT = 5),
  data = data
)
```

    ## 

    ## ================================================================================

    ## [                                Initialising  ...                             ]

    ## ================================================================================

    ## 

    ## [npart=40 ; maxit=100 ; method=spso2011 ; topology=random ; boundary.wall=absorbing2011]

    ## 

    ## [ user-definitions in control: maxit=100 ; REPORT=5 ]

    ## 

    ## 

    ## ================================================================================

    ## [ Writing the 'PSO_logfile.txt' file ...                                       ]

    ## ================================================================================

    ## 

    ## ================================================================================

    ## [                                 Running  PSO ...                             ]

    ## ================================================================================

    ## 

    ## iter:  5  Gbest: 6.606E+00  Gbest_rate:  0.00%  Iter_best_fit: 9.789E+00  nSwarm_Radius: 3.31E-01  |g-mean(p)|/mean(p): 98.64%

    ## iter: 10  Gbest: 3.531E+00  Gbest_rate:  0.00%  Iter_best_fit: 5.648E+00  nSwarm_Radius: 2.33E-01  |g-mean(p)|/mean(p): 71.63%

    ## iter: 15  Gbest: 3.024E+00  Gbest_rate:  0.00%  Iter_best_fit: 3.737E+00  nSwarm_Radius: 2.16E-01  |g-mean(p)|/mean(p): 65.56%

    ## iter: 20  Gbest: 3.024E+00  Gbest_rate:  0.00%  Iter_best_fit: 3.957E+00  nSwarm_Radius: 1.13E-01  |g-mean(p)|/mean(p): 52.06%

    ## iter: 25  Gbest: 1.834E+00  Gbest_rate:  0.00%  Iter_best_fit: 2.150E+00  nSwarm_Radius: 8.61E-02  |g-mean(p)|/mean(p): 64.19%

    ## iter: 30  Gbest: 1.814E+00  Gbest_rate:  1.07%  Iter_best_fit: 1.814E+00  nSwarm_Radius: 3.51E-02  |g-mean(p)|/mean(p): 56.16%

    ## iter: 35  Gbest: 1.427E+00  Gbest_rate: 13.37%  Iter_best_fit: 1.427E+00  nSwarm_Radius: 3.85E-02  |g-mean(p)|/mean(p): 58.39%

    ## iter: 40  Gbest: 1.266E+00  Gbest_rate:  0.00%  Iter_best_fit: 1.295E+00  nSwarm_Radius: 3.71E-02  |g-mean(p)|/mean(p): 52.37%

    ## iter: 45  Gbest: 1.119E+00  Gbest_rate:  0.00%  Iter_best_fit: 1.131E+00  nSwarm_Radius: 2.81E-02  |g-mean(p)|/mean(p): 44.35%

    ## iter: 50  Gbest: 8.196E-01  Gbest_rate:  0.00%  Iter_best_fit: 8.938E-01  nSwarm_Radius: 2.16E-02  |g-mean(p)|/mean(p): 47.29%

    ## iter: 55  Gbest: 7.554E-01  Gbest_rate:  1.52%  Iter_best_fit: 7.554E-01  nSwarm_Radius: 1.31E-02  |g-mean(p)|/mean(p): 37.85%

    ## iter: 60  Gbest: 7.201E-01  Gbest_rate:  0.00%  Iter_best_fit: 7.227E-01  nSwarm_Radius: 9.71E-03  |g-mean(p)|/mean(p): 27.18%

    ## iter: 65  Gbest: 6.983E-01  Gbest_rate:  0.00%  Iter_best_fit: 7.069E-01  nSwarm_Radius: 4.99E-03  |g-mean(p)|/mean(p): 21.55%

    ## iter: 70  Gbest: 6.887E-01  Gbest_rate:  0.40%  Iter_best_fit: 6.887E-01  nSwarm_Radius: 3.40E-03  |g-mean(p)|/mean(p):  3.88%

    ## iter: 75  Gbest: 6.871E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.881E-01  nSwarm_Radius: 2.32E-03  |g-mean(p)|/mean(p):  2.02%

    ## iter: 80  Gbest: 6.866E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.869E-01  nSwarm_Radius: 1.29E-03  |g-mean(p)|/mean(p):  1.15%

    ## iter: 85  Gbest: 6.864E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.865E-01  nSwarm_Radius: 1.03E-03  |g-mean(p)|/mean(p):  0.36%

    ## iter: 90  Gbest: 6.862E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.862E-01  nSwarm_Radius: 5.75E-04  |g-mean(p)|/mean(p):  0.22%

    ## iter: 95  Gbest: 6.861E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.861E-01  nSwarm_Radius: 3.39E-04  |g-mean(p)|/mean(p):  0.17%

    ## iter:100  Gbest: 6.861E-01  Gbest_rate:  0.00%  Iter_best_fit: 6.861E-01  nSwarm_Radius: 2.73E-04  |g-mean(p)|/mean(p):  0.05%

    ## 

    ## [ Writing output files... ]

    ## 

    ##                                     |

    ## ================================================================================

    ## [                          Creating the R output ...                           ]

    ## ================================================================================

## `nloptr::CRS`: controlled random search

``` r
set.seed(11000222)

crs <- crs2lm(
  x0 = theta, 
  fn=sim_ofv, 
  lower = lower, 
  upper = upper, 
  data=data,
  maxeval=2500
)
```

# Compare optimization methods

``` r
results <- list(theta, fit$par, fit1b$par, fit1c$par, fit2$optim$bestmem, fit3$par, fit4$par, crs$par)

results <- map(results, set_names, nm = names(theta))

results <- map(results, exp)

tibble(
  method = c("initial", "newuoa", "nelder", "nelder2", "DEoptim", "SA", "PSO","CRS"),
  fbCLintall = map_dbl(results, "fbCLintall"), 
  ikiu = map_dbl(results, "ikiu"), 
  fbile = map_dbl(results, "fbile"), 
  ka = map_dbl(results, "ka"), 
  ktr = map_dbl(results, "ktr")
) %>% mutate_if(is.double, list(signif), digits = 4)
```

    ## # A tibble: 8 x 6
    ##   method  fbCLintall   ikiu fbile    ka   ktr
    ##   <chr>        <dbl>  <dbl> <dbl> <dbl> <dbl>
    ## 1 initial      1.2   1.2    0.9   0.1   0.1  
    ## 2 newuoa       0.815 0.0110 0.344 0.989 0.690
    ## 3 nelder       0.814 0.0110 0.343 0.976 0.691
    ## 4 nelder2      0.815 0.0110 0.344 0.989 0.690
    ## 5 DEoptim      0.815 0.0110 0.344 0.989 0.689
    ## 6 SA           0.815 0.0110 0.344 0.989 0.688
    ## 7 PSO          0.815 0.0110 0.344 0.988 0.689
    ## 8 CRS          0.815 0.0110 0.344 0.988 0.690

``` r
value0 <- sim_ofv(theta,data)
results <- c(value0, fit$value, fit1b$value, fit1c$value, fit2$optim$bestval, fit3$value, fit4$value, crs$value)

tibble(
  method = c("initial", "newuoa", "nelder", "nelder2", "DEoptim", "SA", "PSO","CRS"),
  value = results
) 
```

    ## # A tibble: 8 x 2
    ##   method   value
    ##   <chr>    <dbl>
    ## 1 initial 11.0  
    ## 2 newuoa   0.686
    ## 3 nelder   0.687
    ## 4 nelder2  0.686
    ## 5 DEoptim  0.686
    ## 6 SA       0.686
    ## 7 PSO      0.686
    ## 8 CRS      0.686
