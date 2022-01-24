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
