
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flatheadresp

<!-- badges: start -->
<!-- badges: end -->

The goal of flatheadresp is to streamline working with respirometry data
that is produced from the AquaResp software (i.e., used in the TAF lab
for flathead metabolic rate experiments), including importing AquaResp
metadata and experimental data files into R, correcting mass-specific
MO2 for post-experiment body mass measurements or other fixes, and
calculating MO_2 from the linear regression slope of raw O2 data and
correlation of O_2 ~ time. This allows calculating accurate coefficients
of determination (R^2) to overcome the bug in the AquaResp v3.0 version
where any missing O_2 values resulted in erroneous R^2 values for MO2
measurements

## Installation

You can install the development version of flatheadresp from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("bwwolfe/flatheadresp")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(flatheadresp)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
