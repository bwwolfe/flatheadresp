
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flatheadresp

<!-- badges: start -->
<!-- badges: end -->

The goal of flatheadresp is to streamline working with respirometry data
that is produced from the AquaResp software (i.e., used in the IMAS
Taroona Aquaculture Facility lab for flathead metabolic rate
experiments), including importing AquaResp metadata and experimental
data files into R, correcting mass-specific $\dot M O_2$ for
post-experiment body mass measurements or other fixes, and calculating
$\dot M O_2$ from the linear regression slope of raw O2 data and
correlation of O<sub>2</sub> ~ time.

AquaResp provide mass-specific oxygen consumption ($\dot M O_2$)
estimates for animals from oxygen concentration in each experimental
chamber during the sealed measurement portion of each flush-wait-measure
cycle. $\dot M O_2$ is provided in units of mg O₂/kg/hr (milligrams of
oxygen per kilogram of fish bodyweight per hour). $\dot M O_2$ are
calculated from the slope of a regression line fit to (declining) oxygen
concentration within the sealed chamber over time. Oxygen concentration
is sampled from the probe in the experimental chamber every second. As
is standard for respirometry work, $\dot M O_2$ estimates are quality
controlled with the coefficients of determination ($R^2$) calculated
with Pearson’s correlation, with $\dot M O_2$ values of $R^2 \geq 0.95$
considered precise enough. However, there is a bug in AquaResp where
occasionally, an O<sub>2</sub> concentration measurement is missed and
recorded as a 0 instead of NA, and then this causes an erroneously low
$R^2$ (Another source of error in this calculation is the Firesting
probes used in the IMAS TAF lab also sometimes register a -300 for PO2
if the optical probe is removed from its housing). This bug will also
create a slightly smaller magnitude slope and thus $\dot M O_2$ estimate
as well.

The functionality in this package allows for calculating accurate $R^2$
to overcome the bug in the AquaResp v3.0 version where any missing
O<sub>2</sub> values resulted in erroneous $R^2$ values for MO2
measurements. Doing this manually requires opening each measurement
cycle’s file (sometimes hundreds, e.g. SDA experiments) and
hand-calculating MO2 values.

Another useful function will allow animal masses and densities to be
corrected (AquaResp assume 1 kg/L which would be neutral buoyancy in
freshwater and slightly positively buoyant in saltwater, however benthic
Sand Flathead are around 1.1 kg/L).

`flatheadresp` also provides tools for monitoring experiments in
real-time.

## Installation

You can install the development version of flatheadresp from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("bwwolfe/flatheadresp")
```

## How to use

``` r
library(flatheadresp)
```

### Load Example Data

The package ships an example AquaResp experiment under
`inst/extdata/aquaresp_experiment` that is used in the examples in the
documentation. Use `system.file()` to locate it.

``` r
exp_path <- system.file("extdata", "aquaresp_experiment",
                            package = "flatheadresp")
```

The good news is you only need to provide the file path of the
experiment you want to analyse as the `path =` argument, and the
package’s functions will extract the needed files. Here, we’ve saved the
experiment folder path as `exp_path`, and we can provide this to
functions for analysis. If you want to see where the files are saved and
browse them yourself, you can see the directory path with
`print(exp_path)` or just type `exp_path` into the console.

This demo experiment’s directory is structured in the standard format
used by AquaResp and this structure and file naming convention will be
expected by the functions in this package:

``` r
list.files(exp_path)
#>  [1] "All slopes"                  "Experimental information"   
#>  [3] "notes.txt"                   "Oxygen data raw"            
#>  [5] "Summary data ABS resp 1.txt" "Summary data ABS resp 2.txt"
#>  [7] "Summary data ABS resp 3.txt" "Summary data ABS resp 4.txt"
#>  [9] "Summary data resp 1.txt"     "Summary data resp 2.txt"    
#> [11] "Summary data resp 3.txt"     "Summary data resp 4.txt"
```

This should be how AquaResp has automatically saved the experimental
data, if e.g. you have four chambers, but worth noting if there are any
differences in the naming conventions or if folders are renamed, it
probably won’t work (also, there is usually an “Experimental
Information” folder created which was excluded with the example
experiment because it was empty). Herein *‘experiment’* will refer to
that which is recorded in the `exp_path`.

### Getting started

A good starting point is the function `summarise_experiment()`, which
will provide an easy to read summary of the experiment:

``` r
summarise_experiment(exp_path)
#> 
#> ── Experiment Summary ──
#> 
#> Started: 2025-04-26 14:18:54
#> Number of chambers: 4
#> 
#> ── Environmental Conditions
#> Salinity: 35.5 ppt
#> Temperature: 15.6 deg C
#> 
#> ── Fish Mass by Chamber
#> Chamber 1: 0.123 kg
#> Chamber 2: 0.058 kg
#> Chamber 3: 0.127 kg
#> Chamber 4: 0.001 kg
#> 
#> ── Respirometer Volumes
#> All chambers: 3.35 L
#> 
#> ── Cycle Timing
#> Flush: 240 s
#> Wait: 60 s
#> Measurement: 1200 s
#> 52 total cycles.
```

Here we can see there were three fish (Sand Flathead, *Platycephalus
bassensis*) were ran with masses from 58 g to 127 g, and Chamber 4 (ran
as a blank to measure background respiration), had an arbitrary 0.001 kg
mass entered.

To see all of the metadata for all chambers that AquaResp has stored,
get_exp_metadata will provide it collated into a dataframe:

``` r
get_exp_metadata(exp_path)
#>   chamber Experiment.start..UNIX.time Flush.time..s Wait.time..s
#> 1       1                  1745641134           240           60
#> 2       2                  1745641134           240           60
#> 3       3                  1745641134           240           60
#> 4       4                  1745641134           240           60
#>   Measurement.time..s Mass.of.fish..kg Volume.respirometer..L
#> 1                1200            0.123                   3.35
#> 2                1200            0.058                   3.35
#> 3                1200            0.127                   3.35
#> 4                1200            0.001                   3.35
#>   Real.volume..vresp...vfish...neutrally.bouyant...L Salinity Temperature
#> 1                                              3.227     35.5        15.6
#> 2                                              3.292     35.5        15.6
#> 3                                              3.223     35.5        15.6
#> 4                                              3.349     35.5        15.6
#>   Oxygen.solubilty..mg.O2...L           exp_start
#> 1                      8.0623 2025-04-26 14:18:54
#> 2                      8.0623 2025-04-26 14:18:54
#> 3                      8.0623 2025-04-26 14:18:54
#> 4                      8.0623 2025-04-26 14:18:54
```

### Load original, uncorrected AquaResp data with `get_exp_mo2s()`

To get the experimental MO_2 data that AquaResp recorded/calculated, use
the function `get_exp_mo2s()`, which returns it as a dataframe collated
for all chambers and cycles. No corrections are applied at this point.
There is `chambers` argument which will allow for subsetting by chamber,
e.g. to exclude the blank respirometer chamber 4 (here just the first
few rows are shown):

``` r

get_exp_mo2s(exp_path, chambers = c(1,2,3)) # or chamber = -4
```

    #>     cycle chamber          Clock.TIME TIME.HOURS  TIME.UNIX        MO2
    #> 1       1       1 2025-04-26 16:02:07  0.3500000 1745648526 124.189246
    #> 53      1       2 2025-04-26 16:02:07  0.3502778 1745648526  88.729289
    #> 105     1       3 2025-04-26 16:02:07  0.3502778 1745648526   4.523504
    #> 2       2       1 2025-04-26 16:27:07  0.7666667 1745650026 132.901748
    #>             SLOPE Intercept  Pearson.R       R.2 P      Std.Err
    #> 1   -0.0163090768  96.79758 -0.9992922 0.9985849 0 1.774512e-05
    #> 53  -0.0053861000 101.88323 -0.9956077 0.9912346 0 1.463943e-05
    #> 105 -0.0006141263 100.63985 -0.9337680 0.8719226 0 6.803116e-06
    #> 2   -0.0174532408 100.40244 -0.9325148 0.8695838 0 1.953616e-04
    #>     Measurement.duration.seconds   avg.po2 median.po2 minimum.po2 max.po2
    #> 1                           1199  87.01214     86.936      77.252  96.350
    #> 53                          1199  98.65157     98.804      95.181 101.435
    #> 105                         1199 100.27138    100.256      99.780 100.747
    #> 2                           1199  89.93050     89.839       0.000 100.297
    #>     delta.po2 oxygen.solubility ratio.vreal.fish
    #> 1      19.098          8.062287         26.23577
    #> 53      6.254          8.062287         56.75862
    #> 105     0.967          8.062287         25.37795
    #> 2     100.297          8.062287         26.23577
    #>     total.experiment.duration.hours  minutes seconds       days
    #> 1                         0.3500000 21.00000    1260 0.01458333
    #> 53                        0.3502778 21.01667    1261 0.01459491
    #> 105                       0.3502778 21.01667    1261 0.01459491
    #> 2                         0.7666667 46.00000    2760 0.03194444

## Recalculate experimental MO_2s with `calc_exp_mo2s()`

Generally speaking, however, you will want to use the main workhorse
function of `flatheadresp`, `calc_exp_mo2s()`. This function
recalculates all of the MO2 and related regression and correlation stats
after removing spurious missing PO2 values from the data, and also fixes
an AquaResp bug that results in spurious MO2 calculations for any
measurements that span midnight on the computer clock time. The function
will also return the original AquaResp data, with `_uncorrected`
appended to their names.

``` r

mo2_data <- calc_exp_mo2s(exp_path, chambers = -4)
#> 
#> ── MO2 Summary for selected chambers: 1, 2, 3 ──
#> 
#> With a po2 and r2 difference tolerance of 0.001:
#> • 52  cycle(s) total.
#> • 10  cycle(s) had corrected values.
#> • 9   cycle(s) had an uncorrected minimum.po2 <= 0.
#> • 1   cycle(s) had an uncorrected minimum.po2 < 0.
#> • 9   cycle(s) R^2 were corrected (changed > 0.001).
#> 
#> Max percent of pO2 measurements missing (<= 0) in a cycle: 1.00%
#> 
#> 0 cycle(s) in which all selected chambers' corrected R.2 < 0.95.
#> 
#> Cycles with corrected R.2 < 0.95 by chamber:
#> Chamber 1: 7 cycle(s) (18, 19, 20, 23, 27, 39, 45)
#> Chamber 2: 16 cycle(s) (22, 29, 30, 32, 35, 37, 38, 40, 41, 42 ...)
#> Chamber 3: 8 cycle(s) (1, 2, 3, 4, 5, 19, 22, 50)
```

`calc_exp_mo2s()` provides a summary of corrections and some metrics
that are useful to flag data quality issues. Details about these metrics
and their use can be found in the `?calc_exp_mo2s` documentation Details
section.

The function itself returns all of the original columns that AquaResp
saves to file and `get_exp_mo2s()` returns, however the recalculated
values take the place in the first few columns, followed by the original
values which are at the end of the dataframe’s columns:

``` r
head(mo2_data, n = 1)
#>   cycle chamber          Clock.TIME TIME.HOURS  TIME.UNIX    MO2       SLOPE
#> 1     1       1 2025-04-26 16:02:07       0.35 1745648526 124.19 -0.01630915
#>   Intercept  Pearson.R       R.2 P      Std.Err Measurement.duration.seconds
#> 1  96.78197 -0.9992926 0.9985857 0 1.774011e-05                         1199
#>    avg.po2 median.po2 minimum.po2 max.po2 delta.po2 oxygen.solubility
#> 1 87.01214     86.936      77.252   96.35    19.098          8.062287
#>   ratio.vreal.fish total.experiment.duration.hours minutes seconds       days
#> 1         26.23577                            0.35      21    1260 0.01458333
#>   pct0 corrected MO2_uncorrected SLOPE_uncorrected Intercept_uncorrected
#> 1    0     FALSE        124.1892       -0.01630908              96.79758
#>   Pearson.R_uncorrected R.2_uncorrected P_uncorrected Std.Err_uncorrected
#> 1            -0.9992922       0.9985849             0        1.774512e-05
#>   avg.po2_uncorrected median.po2_uncorrected minimum.po2_uncorrected
#> 1            87.01214                 86.936                  77.252
#>   max.po2_uncorrected delta.po2_uncorrected
#> 1               96.35                19.098
```

### Updating animal masses and/or densities with `fix_exp_mo2s()`

Sometimes an animal mass is incorrectly entered for an experiment and
needs to be updated. However, this is not as simple as multiplying the
mass-specific MO_2 by the original mass and re-dividing by the new mass,
as the animal mass is also used by AquaResp to calculate the real volume
of the respirometer (assuming neutral buoyancy of the animal in
freshwater e.g. 1 kg of fish = 1 L volume). The `fix_exp_mo2s()`
function will take new masses (in kg) for one or more of the
experiment’s chambers and perform this recalculation.

Also, the implicitly assumed density that AquaResp uses (1 kg/L) would
probably be a safe assumption for freshwater animals that are neutrally
buoyant, but is already slightly off for neutral buoyancy in seawater
(~1.025 kg/L), and will also be incorrect for benthic species that sink.
For example, Sand Flathead are (based on a few tests by the author)
about 1.1 kg/L. So the `fix_exp_mo2()` function allows for both masses
and densities to be updated and will correct the values accordingly.

For example, let’s say after the experiment it was found the 2nd and 3rd
chambers’ masses were accidentally swapped when entered into AquaResp
and need to be swapped back:

``` r

invisible(
  fix_exp_mo2s(path = exp_path,
  new_masses = c(`2` = 0.127, `3` = 0.058)) # name each new mass with the chamber number
  )                                         # in backticks or use NA for chambers that
#> 
#> ── Applying MO2 corrections per chamber ──
#> 
#> Chamber 1: 0.123 kg -> not changed
#> Respirometer ratio (rRespFish): 26.24 -> 26.24 L/kg
#> Mean MO2 % difference: 0%
#> 
#> Chamber 2: 0.058 kg -> 0.1270 kg
#> Respirometer ratio (rRespFish): 56.76 -> 25.38 L/kg
#> Mean MO2 % difference: -55.29%
#> 
#> Chamber 3: 0.127 kg -> 0.0580 kg
#> Respirometer ratio (rRespFish): 25.38 -> 56.76 L/kg
#> Mean MO2 % difference: 123.65%
#> 
#> Chamber 4: 0.001 kg -> not changed
#> Respirometer ratio (rRespFish): 3349 -> 3349 L/kg
#> Mean MO2 % difference: 0%
                                            # don't need to be updated
                                            # i.e. this works the same: c(NA, 0.127, 0.058, NA)
```

Or, we can update all of the chambers’ values with the appropriate
density for the species (the blank chamber 4 mass is changed to 0.1 kg
as well for illustration):

``` r
fixed_mo2s <-
  fix_exp_mo2s(path = exp_path,
               new_masses = c(`4` = 0.1), #single named value in vector - only applies to that chamber
               new_densities = 1.1) #if only scalar value is provided, it is applied to all chambers 
#> 
#> ── Applying MO2 corrections per chamber ──
#> 
#> Chamber 1: 0.123 kg -> not changed
#> Respirometer ratio (rRespFish): 26.24 -> 26.33 L/kg
#> Density: original = 1.000 kg/L, new = 1.100 kg/L
#> Mean MO2 % difference: 0.35%
#> 
#> Chamber 2: 0.058 kg -> not changed
#> Respirometer ratio (rRespFish): 56.76 -> 56.85 L/kg
#> Density: original = 1.000 kg/L, new = 1.100 kg/L
#> Mean MO2 % difference: 0.16%
#> 
#> Chamber 3: 0.127 kg -> not changed
#> Respirometer ratio (rRespFish): 25.38 -> 25.47 L/kg
#> Density: original = 1.000 kg/L, new = 1.100 kg/L
#> Mean MO2 % difference: 0.36%
#> 
#> Chamber 4: 0.001 kg -> 0.1000 kg
#> Respirometer ratio (rRespFish): 3349 -> 32.59 L/kg
#> Density: original = 1.000 kg/L, new = 1.100 kg/L
#> Mean MO2 % difference: -99.03%
```

As you can see, the result of including an appropriate density value
varies with the mass of the animal, so it could be an important thing to
correct for if the species is much different from the default 1 kg/L.
Both mass and density can be updated simultaneously.

Note in the below output of `fix_exp_mo2s()` with most of the columns
omitted for brevity — the **fixed** MO2 is now in the `MO2` column,
while the original value returned by `calc_exp_mo2s()` is retained but
renamed `MO2_orig`. Also, new columns for the original and new mass and
density (if applicable) are appended at the end of the columns.

``` r
head(fixed_mo2s[,c(1,2,6,7,40:43)], n = 5)
#>     cycle chamber        MO2    MO2_orig Mass_orig Mass_new Density_orig
#> 1       1       1 124.620356   124.19003     0.123       NA            1
#> 53      1       2  88.871749    88.72963     0.058       NA            1
#> 105     1       3   4.539734     4.52353     0.127       NA            1
#> 157     1       4 -25.554630 -2625.96096     0.001      0.1            1
#> 2       2       1 131.060852   130.60828     0.123       NA            1
#>     Density_new
#> 1           1.1
#> 53          1.1
#> 105         1.1
#> 157         1.1
#> 2           1.1
```

## Plotting:

Sometimes for diagnosing issues, it is useful to plot the raw oxygen
measurements across a cycle’s measurement window:

``` r
plot_cycle_po2(cycle_number = 9, path = exp_path)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

Alternatively, you can read in a single cycle’s data with `read_cycle()`
and plot with ggplot2. Convert to long format for plotting:

``` r
library(tidyr)

cycle9 <- read_cycle(cycle_number = 9, path = exp_path) 

cycle_long <- cycle9 %>%
  tidyr::pivot_longer(cols = starts_with("ch"),
               names_to = "chamber",
               names_pattern = "ch(\\d)\\.po2",
               values_to = "po2")
```

Plot partial pressure of oxygen (PO<sub>2</sub>) in the chambers during
the cycle:

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.5.1
ggplot(cycle_long,
       aes(Unix.Time - min(Unix.Time), po2, colour = chamber)) +
  geom_path() +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic(12)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

### Plotting MO<sub>2</sub> over time

Visualize MO<sub>2</sub> across cycles with `plot_exp()`:

``` r

plot_exp(path = exp_path, chambers = -4) #chamber 4 (blank) omitted
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the flatheadresp package.
#>   Please report the issue at <https://github.com/bwwolfe/flatheadresp/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

The setting `show_cycle_labels = FALSE` will remove the cycle number
labels if they are getting too busy. Note that the output is a ggplot,
and so can be customised further by adding ggplot2 layers to it,
e.g. theme() to change theme elements.

## Monitoring

For quick quality control checks and monitoring ongoing experiments, the
`print_exp_mo2s()` function will provide an overview of both MO2 values
and R^2s.

``` r
print_exp_mo2s(exp_path)
#>                  Ch1        Ch2        Ch3        Ch4   |     Start Time             
#> Cycle 1   |     124.19      88.73       4.52   -2625.93 | 2025-04-26 16:02:07      
#> Cycle 2   |     132.90     113.66       4.07    -157.10 | 2025-04-26 16:27:07      
#> Cycle 3   |     120.47     109.80       0.07     272.57 | 2025-04-26 16:52:07      
#> Cycle 4   |      98.28     108.07       4.96     564.54 | 2025-04-26 17:17:08      
#> Cycle 5   |      86.05      99.45      51.46    1127.53 | 2025-04-26 17:42:07      
#> Cycle 6   |      86.68      93.45      79.90     853.11 | 2025-04-26 18:07:07      
#> Cycle 7   |     103.29      95.38      62.93     427.43 | 2025-04-26 18:32:07      
#> Cycle 8   |      83.76      98.63      54.61     323.19 | 2025-04-26 18:57:07      
#> Cycle 9   |      59.79     100.46      51.18     204.35 | 2025-04-26 19:22:07      
#> Cycle 10  |      66.60     104.34      50.60     275.28 | 2025-04-26 19:47:07      
#> Cycle 11  |      46.23      98.80      63.93     828.97 | 2025-04-26 20:12:07      
#> Cycle 12  |      59.40      95.90      47.71     381.05 | 2025-04-26 20:37:07      
#> Cycle 13  |      70.63      97.60      42.39     113.98 | 2025-04-26 21:02:07      
#> Cycle 14  |      81.84      97.16      67.07     -24.64 | 2025-04-26 21:27:07      
#> Cycle 15  |      52.27      95.83      41.19      25.87 | 2025-04-26 21:52:07      
#> Cycle 16  |      36.34      87.05      39.57      21.00 | 2025-04-26 22:17:08      
#> Cycle 17  |      54.46      95.28      36.79     610.08 | 2025-04-26 22:42:07      
#> Cycle 18  |      38.77      61.33      34.70     178.55 | 2025-04-26 23:07:08      
#> Cycle 19  |      41.29      62.66      48.51     -55.98 | 2025-04-26 23:32:08      
#> Cycle 20  |      36.23      52.42      74.19     -61.83 | 2025-04-26 23:57:08      
#> Cycle 21  |      42.74      52.49      47.01    -157.19 | 2025-04-27 00:22:08      
#> Cycle 22  |      30.78      62.27      49.37    -426.20 | 2025-04-27 00:47:08      
#> Cycle 23  |      55.81      39.32      85.85     271.24 | 2025-04-27 01:12:08      
#> Cycle 24  |      47.55      34.09      43.48     184.76 | 2025-04-27 01:37:09      
#> Cycle 25  |      31.84      31.81      45.31    -156.37 | 2025-04-27 02:02:08      
#> Cycle 26  |      34.36      38.60      49.32     103.81 | 2025-04-27 02:27:08      
#> Cycle 27  |      37.68      35.18      45.52    -204.60 | 2025-04-27 02:52:08      
#> Cycle 28  |      33.01      38.64      94.00     -82.42 | 2025-04-27 03:17:08      
#> Cycle 29  |      30.92      35.09      77.92     502.05 | 2025-04-27 03:42:09      
#> Cycle 30  |      22.15      25.79      58.20      54.17 | 2025-04-27 04:07:08      
#> Cycle 31  |      22.89      28.44      66.09    -140.52 | 2025-04-27 04:32:08      
#> Cycle 32  |      21.50      34.29      85.57    -242.99 | 2025-04-27 04:57:08      
#> Cycle 33  |      22.30      73.48     137.25    -211.91 | 2025-04-27 05:22:09      
#> Cycle 34  |      23.73      43.47      92.29    -111.58 | 2025-04-27 05:47:09      
#> Cycle 35  |      22.63      34.42      64.60     416.33 | 2025-04-27 06:12:09      
#> Cycle 36  |      17.68      20.62      54.02    -224.83 | 2025-04-27 06:37:09      
#> Cycle 37  |      19.88      27.68      55.58    -156.69 | 2025-04-27 07:02:09      
#> Cycle 38  |      19.44      33.80      52.36    -225.43 | 2025-04-27 07:27:09      
#> Cycle 39  |      27.77      34.59      47.08    -174.29 | 2025-04-27 07:52:09      
#> Cycle 40  |      23.28      33.75      47.52    -151.73 | 2025-04-27 08:17:09      
#> Cycle 41  |      21.88      31.72      36.04     437.53 | 2025-04-27 08:42:09      
#> Cycle 42  |      20.84      38.02      35.91     438.60 | 2025-04-27 09:07:09      
#> Cycle 43  |      20.00      28.66      33.44    -187.84 | 2025-04-27 09:32:09      
#> Cycle 44  |      23.04      31.08      38.85    -222.35 | 2025-04-27 09:57:09      
#> Cycle 45  |      28.94      31.10      40.71    -233.34 | 2025-04-27 10:22:09      
#> Cycle 46  |      22.40      32.67      36.54    -164.52 | 2025-04-27 10:47:09      
#> Cycle 47  |      20.65      31.69      33.90     139.55 | 2025-04-27 11:12:09      
#> Cycle 48  |      20.42      28.07      33.74     255.32 | 2025-04-27 11:37:10      
#> Cycle 49  |      17.09      21.52      33.67      50.99 | 2025-04-27 12:02:10      
#> Cycle 50  |      18.69      26.25      42.76    -188.53 | 2025-04-27 12:27:10      
#> Cycle 51  |      20.27      28.05      34.39    -235.79 | 2025-04-27 12:52:10      
#> Cycle 52  |      29.02      31.42      34.56    -199.99 | 2025-04-27 13:17:10      
#> 
#>                                Legend                               
#>    min                                                         max
#>   -2625.93 ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 1127.53
#>   (no colour scaling in low-colour mode)
#>   Flags:  R2 < 0.95 → value; missing pO2s → value; both → value
#>   Counts: R2 < 0.95 = 99,     missing pO2s = 33,     both = 33
```

<!-- ## Cycle Summary -->
<!-- Plot PO₂ min/max per cycle: -->
<!-- ```{r, fig.width = 8} -->
<!-- mo2_data <- calc_exp_mo2s(path = exp_path, report = FALSE) -->
<!-- ggplot(mo2_data, -->
<!--        aes(x = cycle, colour = as.factor(chamber))) + -->
<!--   geom_segment(aes(y = minimum.po2, yend = max.po2), -->
<!--                position = position_dodge2(width = .6), -->
<!--                linewidth = 0.8) + -->
<!--   theme_classic(16) + -->
<!--   scale_colour_brewer("Chamber", palette = "Dark2")+ -->
<!--   scale_y_continuous(expression(p*O[2]~"range during cycle"~`(`*'%'~O[2]~"sat"*`)`)) + -->
<!--   labs(x = "Cycle #") -->
<!-- ``` -->
