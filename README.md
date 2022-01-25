Clinical Trial Simulation
================
Adam Lang
1/24/2022

-   [Overview](#overview)
-   [Installation](#installation)
-   [SampleSizeEstimation Function](#samplesizeestimation-function)
-   [Example](#example)
    -   [Data](#data)
    -   [Running Model](#running-model)
-   [Type I Error](#type-i-error)

Overview
========

Clinical Trial Simulation allows a user to estimate the statistical
power of a pre-defined treatment effect at a set of sample sizes. This
may be useful in informing recruitment of future clinical trials. Using
a Linear-Mixed-Effects Model (LME) trained on pilot data, a user can
determine the power of a treatment effect at a sample size via Monte
Carlo Simulation.

<br>

Installation
============

Both ***SampleSizeHelpers.R*** and ***SampleSizeEstimationFunction.R***
must be downloaded. The following packages must be installed:

``` r
install.packages("plyr")
install.packages("dplyr")
install.packages("lme4")
install.packages("nlme")
install.packages("simr")
install.packages("stringr")
install.packages("matrixStats")
install.packages("lmerTest")
install.packages("MASS") 
install.packages("splitstackshape") 
install.packages("cramer")
```

<br>

SampleSizeEstimation Function
=============================

**SampleSizeEstimation** requires the following arguments:

-   ***model*** A fitted LME model (must be a lmer of package
    ***lme4***). This model is assumed to have a subject specific random
    intercept or a subject specific random intercept and slope specified
    with parentheses. <code> Ex: “(1\|Subject)” “(1+Time\|Subject)”
    </code>

-   ***parameter*** The name of the rate of change parameter of
    interest. <code> Ex: “time”, “t”, “Months” </code>

-   ***pct.change*** The percent change in the parameter of interest
    <code> Ex: .05 (50% improvement) (Either pct.change or delta must be
    specified, but not both, defaults to NULL) </code>

-   ***delta*** The change in the pilot estimate of the parameter of
    interest. <code> (Either pct.change or delta must be specified, but
    not both, defaults to NULL) </code>

-   ***time*** Numeric vector of timepoints. <code> Ex: c(0, .5, 1,
    1.5, 2) </code>

-   ***data*** Pilot data used to fit pilot model

-   ***sample\_sizes*** A numeric vector of sample sizes per arm to
    calculate power. <code> Ex: c(100, 150, 200, 250) </code>

-   ***nsim*** Number of iterations to run at each sample size <code>
    (defaults to 500) </code>

-   ***sig.level*** Type I Error rate <code> defaults to 0.05 </code>

-   ***verbose*** Print model fit progress <code> defaults to TRUE
    </code>

-   ***balance.covariates*** Character vector of additional covariates
    to balance while running simulation <code> defaults to NULL </code>

Example
=======

Data
----

``` r
#Example Data
library(simstudy)

set.seed(123) #for reproducibility

#Creating simulated covariates
def <- defData(     varname = "Int",  dist = "normal",      formula = 25 ,  variance = 7.5, id = "id")
def <- defData(def, varname = "Num1", dist = "normal",      formula = 8  ,  variance = 2, id = "id")
def <- defData(def, varname = "Num2", dist = "normal",      formula = 4   , variance = 1, id = "id")
def <- defData(def, varname = "Cat1", dist = "categorical", formula = ".25 ;.25;.5",      id = "id")
def <- defData(def, varname = "Cat2", dist = "categorical", formula = ".5;.5",            id = "id")
def <- genData(500, def)

#Adding random effects
def <- addCorData(def, idname   = "id", mu  = c(0, 0), sigma   = c(10, 1), rho = .5, corstr = "cs", cnames = c("a0", "a1"))
def <- addPeriods(def, nPeriods = 5, idvars = "id", timevarName = "Time")
def$timeID <- NULL

#Create simulated outcome based on covariates
def2 <- defDataAdd(varname  = "Outcome", 
                   formula  = "(Int + a0) + 5*Num1 + 3.5*Num2 + .5*Cat1 - 5*Cat2 + (-2.5 + a1)*period", 
                   variance = 10)
def  <- addColumns(def2, def)


def           <- def[,c("id", "period", "Outcome", "Num1", "Num2", "Cat1", "Cat2")]
colnames(def) <- c("Subject", "Time", "Outcome", "Num1", "Num2", "Cat1", "Cat2")
def$Subject   <- factor(def$Subject)
def$Cat1      <- factor(def$Cat1)
def$Cat2      <- factor(def$Cat2)
def           <- as.data.frame(def)
```

``` r
#fit model
library(lme4)
library(lmerTest)
library(nlme)
library(ggplot2)

simulated.model <- lme4::lmer(Outcome ~ Time + Num1 + Num2 + Cat1 + Cat2 + (1 + Time|Subject), def)
```

<br>

[![plot-Model.png](https://i.postimg.cc/xdnZ30J5/plot-Model.png)](https://postimg.cc/219TzNxL)

Population model (black) vs subject specific curves <br>

Fixed Effects

    ##               Estimate Std. Error       df     t value      Pr(>|t|)
    ## (Intercept) 22.5984099 3.43375705 495.3410   6.5812489  1.190092e-10
    ## Time        -2.4840735 0.06267994 498.9967 -39.6310751 2.978914e-156
    ## Num1         4.8168181 0.33663157 493.9729  14.3088722  4.203944e-39
    ## Num2         3.4474025 0.48413949 493.9726   7.1206803  3.812191e-12
    ## Cat12        0.5841370 1.32697022 493.9725   0.4402036  6.599822e-01
    ## Cat13       -0.2204504 1.18164873 493.9725  -0.1865617  8.520809e-01
    ## Cat22       -4.8025258 0.96057223 493.9725  -4.9996508  7.990883e-07

<br>

Random Effects

    ##  Groups   Name        Std.Dev. Corr 
    ##  Subject  (Intercept) 10.8248       
    ##           Time         1.0106  0.542
    ##  Residual              3.0709

<br>

Running Model
-------------

``` r
set.seed(123)
ss.out <- SampleSizeEstimation(model        = simulated.model,
                               parameter    = "Time",
                               pct.change   = 0.3,
                               time         = c(0, 1, 2, 3, 4),
                               data         = def,
                               nsim         = 500,
                               sample_sizes = c(20, 60, 100),
                               verbose      = FALSE)
```

``` r
ss.out$Power_Per_Sample
```

    ##    mean    ci.low   ci.high
    ## 1 0.344 0.3023909 0.3874665
    ## 2 0.814 0.7770718 0.8471649
    ## 3 0.968 0.9485532 0.9816008

<br>

Type I Error
============

Assessing Type I Error is important in bayesian inference as it may not
be controlled at 5%. In order to estimate Type I error rates, a user can
run the simulation with <code> pct.change </code> or <code> delta
</code> set at 0

``` r
set.seed(123)
ss.out <- SampleSizeEstimation(model        = simulated.model,
                               parameter    = "Time",
                               pct.change   = 0,
                               time         = c(0, 1, 2, 3, 4),
                               data         = def,
                               nsim         = 500,
                               sample_sizes = c(20, 60, 100),
                               verbose      = FALSE)
```

``` r
ss.out$Power_Per_Sample
```

    ##    mean     ci.low    ci.high
    ## 1 0.062 0.04251106 0.08685260
    ## 2 0.066 0.04586336 0.09144290
    ## 3 0.046 0.02938018 0.06822543
