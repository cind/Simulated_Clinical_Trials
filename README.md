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
def <- defData(def, varname = "Num1", dist = "normal", formula = .25  , variance = .1, id = "id")
def <- defData(def, varname = "Num2", dist = "normal", formula = -1  , variance  = .4, id = "id")
def <- defData(def, varname = "Cat1", dist = "categorical", formula = ".25 ;.25 ;.5", id = "id")
def <- defData(def, varname = "Cat2", dist = "categorical", formula = ".5;.5", id = "id")
def <- genData(500, def)

#Adding random effects
def <- addCorData(def, idname   = "id", mu  = c(0, 0), sigma   = c(4, .2), rho = -.5, corstr = "cs", cnames = c("a0", "a1"))
def <- addPeriods(def, nPeriods = 5, idvars = "id", timevarName = "Time")
def$timeID <- NULL

#Create simulated outcome based on covariates
def2 <- defDataAdd(varname = "Outcome", 
                   formula = "(Int + a0) + 5*Num1 + 3.5*Num2 + .5*Cat1 - 5*Cat2 + (-5 + a1)*period", variance = 8)
def <- addColumns(def2, def)


def <- def[,c("id", "period", "Outcome", "Num1", "Num2", "Cat1", "Cat2")]
colnames(def) <- c("Subject", "Time", "Outcome", "Num1", "Num2", "Cat1", "Cat2")
def$Subject   <- factor(def$Subject)
def$Cat1      <- factor(def$Cat1)
def$Cat2      <- factor(def$Cat2)
```

Population model (black) vs subject specific curves

<br>

<img src="https://imgur.com/a/j5DJuhD">

Fixed Effects

    ##                Estimate Std. Error       df       t value      Pr(>|t|)
    ## (Intercept) 45.97099594 0.55068975 520.8599   83.47893849 1.120665e-303
    ## Time        -4.97883439 0.03907265 499.1385 -127.42504884  0.000000e+00
    ## Num1         4.77189369 0.60698716 494.0018    7.86160561  2.397772e-14
    ## Num2         3.57630051 0.30863855 494.0018   11.58734220  1.246759e-27
    ## Cat12        0.04234756 0.53502101 494.0018    0.07915121  9.369444e-01
    ## Cat13        0.07876923 0.47642885 494.0018    0.16533262  8.687498e-01
    ## Cat22       -4.73484111 0.38729304 494.0018  -12.22547430  3.268544e-30

<br>

Random Effects

    ##  Groups   Name        Std.Dev. Corr  
    ##  Subject  (Intercept) 4.324835       
    ##           Time        0.094427 -0.901
    ##  Residual             2.746670
