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
          - [Get some predictions to look at how the fit
            went](#get-some-predictions-to-look-at-how-the-fit-went)
          - [Make some plots](#make-some-plots)
          - [A nicer plot](#a-nicer-plot)
          - [The final objective function value and
            estimates](#the-final-objective-function-value-and-estimates)

# Packages and setup

``` r
library(tidyverse)
library(PKPDmisc)
library(mrgsolve)
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


data
```

    . # A tibble: 23 x 8
    .       ID  time     DV  evid   amt   cmt  type typef     
    .    <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     
    .  1     2  0     NA        1  2000     2     2 Statin+CsA
    .  2     2  1     NA        1    30     1     2 Statin+CsA
    .  3     2  1.49  73.7      0     0     0     2 Statin+CsA
    .  4     2  1.99 102.       0     0     0     2 Statin+CsA
    .  5     2  2.49  59.9      0     0     0     2 Statin+CsA
    .  6     2  3.00  37.6      0     0     0     2 Statin+CsA
    .  7     2  3.97  15.7      0     0     0     2 Statin+CsA
    .  8     2  5.01   9.24     0     0     0     2 Statin+CsA
    .  9     2  6.99   3.54     0     0     0     2 Statin+CsA
    . 10     2  9.01   2.22     0     0     0     2 Statin+CsA
    . # … with 13 more rows

``` r
data %>% filter(evid==1)
```

    . # A tibble: 3 x 8
    .      ID  time    DV  evid   amt   cmt  type typef     
    .   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     
    . 1     2     0    NA     1  2000     2     2 Statin+CsA
    . 2     2     1    NA     1    30     1     2 Statin+CsA
    . 3     1     1    NA     1    30     1     1 Statin

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

<img src="tools_optimization_pbpk_ddi_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

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

<img src="tools_optimization_pbpk_ddi_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
sims %>% 
  group_by(type) %>% 
  summarise(auc = auc_partial(time,CP)) %>% 
  mutate(fold_increase = auc /first(auc))
```

    . # A tibble: 2 x 3
    .   type                 auc fold_increase
    .   <fct>              <dbl>         <dbl>
    . 1 Pitavastatin alone  44.1          1   
    . 2 Pitavastatin + CsA 161.           3.65

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

### Prediction function

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
  
  out <- mod %>% param(p) %>% mrgsim_d(data, output="df")
  
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

    . $par
    . [1] -0.20418269 -4.51446274 -1.06752471 -0.01113143 -0.37149033
    . 
    . $value
    . [1] 0.6860763
    . 
    . $iter
    . [1] 406
    . 
    . $convergence
    . [1] 1
    . 
    . $message
    . [1] "NLOPT_SUCCESS: Generic success return value."

### Get some predictions to look at how the fit went

Recall that our parameters are

``` r
fit$par
```

    . [1] -0.20418269 -4.51446274 -1.06752471 -0.01113143 -0.37149033

We can generate a prediction that matches our data like this

``` r
sim_ofv(fit$par,data = dose, pred = TRUE) %>% filter(time >= 1)
```

    .     ID time         CP
    . 1    2  1.0  0.0000000
    . 2    2  1.0  0.0000000
    . 3    2  1.1 18.1930017
    . 4    2  1.2 28.5554696
    . 5    2  1.3 36.2290752
    . 6    2  1.4 42.1357832
    . 7    2  1.5 46.6374148
    . 8    2  1.6 49.9684588
    . 9    2  1.7 52.3044311
    . 10   2  1.8 53.7866262
    . 11   2  1.9 54.5344715
    . 12   2  2.0 54.6522186
    . 13   2  2.1 54.2326959
    . 14   2  2.2 53.3594832
    . 15   2  2.3 52.1082195
    . 16   2  2.4 50.5474332
    . 17   2  2.5 48.7391048
    . 18   2  2.6 46.7390827
    . 19   2  2.7 44.5974159
    . 20   2  2.8 42.3586469
    . 21   2  2.9 40.0620870
    . 22   2  3.0 37.7420867
    . 23   2  3.1 35.4283119
    . 24   2  3.2 33.1460283
    . 25   2  3.3 30.9163956
    . 26   2  3.4 28.7567722
    . 27   2  3.5 26.6810267
    . 28   2  3.6 24.6998552
    . 29   2  3.7 22.8210988
    . 30   2  3.8 21.0500592
    . 31   2  3.9 19.3898070
    . 32   2  4.0 17.8414807
    . 33   2  4.1 16.4045707
    . 34   2  4.2 15.0771874
    . 35   2  4.3 13.8563108
    . 36   2  4.4 12.7380180
    . 37   2  4.5 11.7176908
    . 38   2  4.6 10.7901997
    . 39   2  4.7  9.9500664
    . 40   2  4.8  9.1916048
    . 41   2  4.9  8.5090405
    . 42   2  5.0  7.8966122
    . 43   2  5.1  7.3486535
    . 44   2  5.2  6.8596597
    . 45   2  5.3  6.4243391
    . 46   2  5.4  6.0376514
    . 47   2  5.5  5.6948353
    . 48   2  5.6  5.3914257
    . 49   2  5.7  5.1232634
    . 50   2  5.8  4.8864976
    . 51   2  5.9  4.6775832
    . 52   2  6.0  4.4932737
    . 53   2  6.1  4.3306106
    . 54   2  6.2  4.1869104
    . 55   2  6.3  4.0597502
    . 56   2  6.4  3.9469513
    . 57   2  6.5  3.8465633
    . 58   2  6.6  3.7568471
    . 59   2  6.7  3.6762582
    . 60   2  6.8  3.6034311
    . 61   2  6.9  3.5371628
    . 62   2  7.0  3.4763990
    . 63   2  7.1  3.4202189
    . 64   2  7.2  3.3678227
    . 65   2  7.3  3.3185192
    . 66   2  7.4  3.2717140
    . 67   2  7.5  3.2268993
    . 68   2  7.6  3.1836442
    . 69   2  7.7  3.1415857
    . 70   2  7.8  3.1004210
    . 71   2  7.9  3.0598997
    . 72   2  8.0  3.0198179
    . 73   2  8.1  2.9800116
    . 74   2  8.2  2.9403519
    . 75   2  8.3  2.9007399
    . 76   2  8.4  2.8611028
    . 77   2  8.5  2.8213895
    . 78   2  8.6  2.7815678
    . 79   2  8.7  2.7416213
    . 80   2  8.8  2.7015467
    . 81   2  8.9  2.6613514
    . 82   2  9.0  2.6210514
    . 83   2  9.1  2.5806701
    . 84   2  9.2  2.5402360
    . 85   2  9.3  2.4997819
    . 86   2  9.4  2.4593433
    . 87   2  9.5  2.4189578
    . 88   2  9.6  2.3786640
    . 89   2  9.7  2.3385011
    . 90   2  9.8  2.2985078
    . 91   2  9.9  2.2587223
    . 92   2 10.0  2.2191814
    . 93   2 10.1  2.1799208
    . 94   2 10.2  2.1409741
    . 95   2 10.3  2.1023731
    . 96   2 10.4  2.0641477
    . 97   2 10.5  2.0263255
    . 98   2 10.6  1.9889318
    . 99   2 10.7  1.9519896
    . 100  2 10.8  1.9155198
    . 101  2 10.9  1.8795409
    . 102  2 11.0  1.8440694
    . 103  2 11.1  1.8091193
    . 104  2 11.2  1.7747028
    . 105  2 11.3  1.7408302
    . 106  2 11.4  1.7075096
    . 107  2 11.5  1.6747476
    . 108  2 11.6  1.6425488
    . 109  2 11.7  1.6109167
    . 110  2 11.8  1.5798527
    . 111  2 11.9  1.5493575
    . 112  2 12.0  1.5194299
    . 113  2 12.1  1.4900680
    . 114  2 12.2  1.4612686
    . 115  2 12.3  1.4330277
    . 116  2 12.4  1.4053404
    . 117  2 12.5  1.3782009
    . 118  2 12.6  1.3516028
    . 119  2 12.7  1.3255392
    . 120  2 12.8  1.3000025
    . 121  2 12.9  1.2749847
    . 122  2 13.0  1.2504776
    . 123  2 13.1  1.2264725
    . 124  2 13.2  1.2029604
    . 125  2 13.3  1.1799321
    . 126  2 13.4  1.1573785
    . 127  2 13.5  1.1352900
    . 128  2 13.6  1.1136572
    . 129  2 13.7  1.0924706
    . 130  2 13.8  1.0717208
    . 131  2 13.9  1.0513982
    . 132  2 14.0  1.0314934
    . 133  1  1.0  0.0000000
    . 134  1  1.0  0.0000000
    . 135  1  1.1  9.0207778
    . 136  1  1.2 12.3918173
    . 137  1  1.3 14.1500160
    . 138  1  1.4 15.0023896
    . 139  1  1.5 15.2638810
    . 140  1  1.6 15.1293726
    . 141  1  1.7 14.7289994
    . 142  1  1.8 14.1531818
    . 143  1  1.9 13.4663922
    . 144  1  2.0 12.7153077
    . 145  1  2.1 11.9339314
    . 146  1  2.2 11.1469867
    . 147  1  2.3 10.3722744
    . 148  1  2.4  9.6223734
    . 149  1  2.5  8.9059013
    . 150  1  2.6  8.2284675
    . 151  1  2.7  7.5934033
    . 152  1  2.8  7.0023254
    . 153  1  2.9  6.4555725
    . 154  1  3.0  5.9525432
    . 155  1  3.1  5.4919589
    . 156  1  3.2  5.0720667
    . 157  1  3.3  4.6907960
    . 158  1  3.4  4.3458784
    . 159  1  3.5  4.0349393
    . 160  1  3.6  3.7555664
    . 161  1  3.7  3.5053621
    . 162  1  3.8  3.2819803
    . 163  1  3.9  3.0831542
    . 164  1  4.0  2.9067152
    . 165  1  4.1  2.7506048
    . 166  1  4.2  2.6128828
    . 167  1  4.3  2.4917306
    . 168  1  4.4  2.3854527
    . 169  1  4.5  2.2924751
    . 170  1  4.6  2.2113427
    . 171  1  4.7  2.1407154
    . 172  1  4.8  2.0793636
    . 173  1  4.9  2.0261626
    . 174  1  5.0  1.9800876
    . 175  1  5.1  1.9402073
    . 176  1  5.2  1.9056789
    . 177  1  5.3  1.8757417
    . 178  1  5.4  1.8497115
    . 179  1  5.5  1.8269757
    . 180  1  5.6  1.8069874
    . 181  1  5.7  1.7892604
    . 182  1  5.8  1.7733646
    . 183  1  5.9  1.7589212
    . 184  1  6.0  1.7455982
    . 185  1  6.1  1.7331068
    . 186  1  6.2  1.7211966
    . 187  1  6.3  1.7096531
    . 188  1  6.4  1.6982933
    . 189  1  6.5  1.6869632
    . 190  1  6.6  1.6755346
    . 191  1  6.7  1.6639024
    . 192  1  6.8  1.6519821
    . 193  1  6.9  1.6397075
    . 194  1  7.0  1.6270284
    . 195  1  7.1  1.6139089
    . 196  1  7.2  1.6003255
    . 197  1  7.3  1.5862653
    . 198  1  7.4  1.5717249
    . 199  1  7.5  1.5567085
    . 200  1  7.6  1.5412270
    . 201  1  7.7  1.5252972
    . 202  1  7.8  1.5089402
    . 203  1  7.9  1.4921807
    . 204  1  8.0  1.4750466
    . 205  1  8.1  1.4575678
    . 206  1  8.2  1.4397759
    . 207  1  8.3  1.4217034
    . 208  1  8.4  1.4033836
    . 209  1  8.5  1.3848496
    . 210  1  8.6  1.3661348
    . 211  1  8.7  1.3472716
    . 212  1  8.8  1.3282921
    . 213  1  8.9  1.3092271
    . 214  1  9.0  1.2901065
    . 215  1  9.1  1.2709588
    . 216  1  9.2  1.2518113
    . 217  1  9.3  1.2326896
    . 218  1  9.4  1.2136180
    . 219  1  9.5  1.1946192
    . 220  1  9.6  1.1757144
    . 221  1  9.7  1.1569231
    . 222  1  9.8  1.1382634
    . 223  1  9.9  1.1197518
    . 224  1 10.0  1.1014035
    . 225  1 10.1  1.0832321
    . 226  1 10.2  1.0652498
    . 227  1 10.3  1.0474675
    . 228  1 10.4  1.0298949
    . 229  1 10.5  1.0125404
    . 230  1 10.6  0.9954114
    . 231  1 10.7  0.9785141
    . 232  1 10.8  0.9618537
    . 233  1 10.9  0.9454344
    . 234  1 11.0  0.9292597
    . 235  1 11.1  0.9133321
    . 236  1 11.2  0.8976535
    . 237  1 11.3  0.8822250
    . 238  1 11.4  0.8670470
    . 239  1 11.5  0.8521196
    . 240  1 11.6  0.8374420
    . 241  1 11.7  0.8230132
    . 242  1 11.8  0.8088318
    . 243  1 11.9  0.7948957
    . 244  1 12.0  0.7812028
    . 245  1 12.1  0.7677506
    . 246  1 12.2  0.7545362
    . 247  1 12.3  0.7415568
    . 248  1 12.4  0.7288090
    . 249  1 12.5  0.7162895
    . 250  1 12.6  0.7039949
    . 251  1 12.7  0.6919216
    . 252  1 12.8  0.6800658
    . 253  1 12.9  0.6684238
    . 254  1 13.0  0.6569919
    . 255  1 13.1  0.6457662
    . 256  1 13.2  0.6347428
    . 257  1 13.3  0.6239181
    . 258  1 13.4  0.6132880
    . 259  1 13.5  0.6028489
    . 260  1 13.6  0.5925970
    . 261  1 13.7  0.5825286
    . 262  1 13.8  0.5726400
    . 263  1 13.9  0.5629276
    . 264  1 14.0  0.5533878

We can also get the predictions under the initial conditions by passing
in `theta` rather than `fit$par`

In the next block, generate

1.  Predictions with the final estimates
2.  Predications with the initial estimates
3.  Observed data to
overlay

<!-- end list -->

``` r
df_pred <- sim_ofv(fit$par, dose, pred=TRUE) %>% mutate(type = typef(ID))
df_init <- sim_ofv(theta,   dose, pred=TRUE) %>% mutate(type = typef(ID))
df_obs <-  mutate(data, type=typef(ID))
```

### Make some plots

``` r
ggplot(df_pred, aes(time,CP)) + 
  geom_line(lwd=1) + 
  geom_point(data = df_obs, aes(time,DV),col="firebrick",size=2) + 
  facet_wrap(~type) + scale_y_log10() 
```

<img src="tools_optimization_pbpk_ddi_files/figure-gfm/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

### A nicer plot

``` r
ggplot(data=df_pred) + 
  geom_line(data=df_init,aes(time,CP,lty="A"), col="black", lwd=0.7) +
  geom_line(aes(time,CP,lty="B"),col="black",lwd=0.7) + 
  geom_point(data=df_obs,aes(time,DV,col=type),size=3) + 
  facet_wrap(~type) + 
  scale_y_continuous(trans="log",breaks=10^seq(-4,4), 
                     limits=c(0.1,100),
                     "Pitavastatin concentration (ng/mL)") +
  scale_x_continuous(name="Time (hours)", breaks=seq(0,14,2)) +
  scale_linetype_manual(values= c(2,1), guide = FALSE,
                        labels=c("Initial estimates", "Final estimates"), name="") +
  theme_bw() + theme(legend.position="top") 
```

<img src="tools_optimization_pbpk_ddi_files/figure-gfm/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

### The final objective function value and estimates

``` r
sim_ofv(fit$par,data=data)
```

    . [1] 0.6860763

``` r
exp(fit$par)
```

    . [1] 0.81531341 0.01094949 0.34385861 0.98893030 0.68970568
