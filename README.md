Clinical Trial Simulation
================
Adam Lang
1/24/2022

-   [Overview](#overview)
-   [Installation](#installation)
-   [Running SampleSizeEstimation
    Function](#running-samplesizeestimation-function)
-   [Example](#example)

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

Running SampleSizeEstimation Function
=====================================

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

-   ***verbose*** Logical print model fit progress <code> defaults to
    TRUE </code>

-   ***balance.covariates*** Character vector of additional covariates
    to balance while running simulation <code> defaults to NULL </code>

Example
=======

``` r
#Example Data
library(simstudy)

set.seed(123) #for reproducibility

#Creating simulated covariates
def <- defData(varname = "Int", dist = "normal", formula = 50 , variance = 1, id = "id")
def <- defData(def, varname = "Num1", dist = "normal", formula = 8  , variance = 2, id = "id")
def <- defData(def, varname = "Num2", dist = "normal", formula = 4   , variance  = 1, id = "id")
def <- defData(def, varname = "Cat1", dist = "categorical", formula = ".25 ;.25;.5", id = "id")
def <- defData(def, varname = "Cat2", dist = "categorical", formula = ".5;.5", id = "id")
def <- genData(500, def)

#Adding random effects
def <- addCorData(def, idname   = "id", mu  = c(0, 0), sigma   = c(10, 1), rho = -.5, corstr = "cs", cnames = c("a0", "a1"))
def <- addPeriods(def, nPeriods = 5, idvars = "id", timevarName = "Time")
def$timeID <- NULL

#Create simulated outcome based on covariates
def2 <- defDataAdd(varname = "Outcome", 
                   formula = "(Int + a0) + 5*Num1 + 3.5*Num2 + .5*Cat1 - 5*Cat2 + (-5 + a1)*period", variance = 5)
def <- addColumns(def2, def)


def <- def[,c("id", "period", "Outcome", "Num1", "Num2", "Cat1", "Cat2")]
colnames(def) <- c("Subject", "Time", "Outcome", "Num1", "Num2", "Cat1", "Cat2")
def$Subject   <- factor(def$Subject)
def$Cat1      <- factor(def$Cat1)
def$Cat2      <- factor(def$Cat2)
def <- as.data.frame(def)
```

Population model (black) vs subject specific curves

``` r
#fit model
library(lme4)
library(lmerTest)
library(nlme)
library(ggplot2)
simulated.model <-lme4::lmer(Outcome ~ Time + Num1 + Num2 + Cat1 + Cat2 + (1 + Time|Subject), def)
```

<br>

[![plot-Model.png](https://i.postimg.cc/7LjRdFBG/plot-Model.png)](https://postimg.cc/kVFwtpPq)

Fixed Effects

    ##               Estimate Std. Error       df     t value      Pr(>|t|)
    ## (Intercept) 44.9209545 3.00215851 499.1830  14.9628857  4.372318e-42
    ## Time        -4.9927456 0.05265967 499.0013 -94.8115563 2.707480e-321
    ## Num1         5.0563904 0.29374314 494.0047  17.2136460  2.295992e-52
    ## Num2         3.7347077 0.42245785 494.0043   8.8404266  1.685706e-17
    ## Cat12       -0.7677797 1.15790802 494.0039  -0.6630748  5.075919e-01
    ## Cat13       -0.7673643 1.03110116 494.0040  -0.7442182  4.570983e-01
    ## Cat22       -4.3884154 0.83819085 494.0039  -5.2355802  2.437315e-07

<br>

Random Effects

    ##  Groups   Name        Std.Dev. Corr  
    ##  Subject  (Intercept) 10.40816       
    ##           Time         0.95656 -0.484
    ##  Residual              2.17144

``` r
set.seed(123)
ss.out <- SampleSizeEstimation(model        = simulated.model,
                               parameter    = "Time",
                               pct.change   = 0.1,
                               time         = c(0, 1, 2, 3, 4),
                               data         = def,
                               nsim         = 100,
                               sample_sizes = c(20, 60, 100),
                               verbose      = FALSE)
```

``` r
ss.out$Power_Per_Sample
```

    ##   mean    ci.low   ci.high
    ## 1 0.20 0.1266556 0.2918427
    ## 2 0.66 0.5584667 0.7517765
    ## 3 0.81 0.7193020 0.8815568
