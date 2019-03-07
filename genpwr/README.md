
<!-- README.md is generated from README.Rmd. Please edit that file -->
### genpwr

The genpwr package for R (&gt;3.5.1) performs power and sample size calculations for genetic association studies and allows for mis-specification of the genetic model. Calculations can be performed for binary (case/control) and continuous outcomes. Power and sample size calculations are possible for genetic effects as well as gene by environment interactions.

### Example

To calculate power to detect an odds ratio of 2 for a 1:1 case control study with 2,000 subjects, assuming an alpha of 0.05, at minor allele frequencies of 0.1, 0.2, and 0.3:

``` r
library(genpwr)
#> Loading required package: ggplot2
#> Loading required package: nleqslv
#> Loading required package: MASS

genpwr.calc(calc = "power", model = "logistic", N = 2000, OR = 2,
            Alpha = 0.05, MAF = c(0.1,0.2,0.3), Case.Rate = 0.5)
#>     Test.Model True.Model MAF OR N_total N_cases N_controls Case.Rate
#> 1     Dominant   Dominant 0.1  2    2000    1000       1000       0.5
#> 3     Dominant   Additive 0.1  2    2000    1000       1000       0.5
#> 5     Dominant  Recessive 0.1  2    2000    1000       1000       0.5
#> 7     Dominant   Dominant 0.2  2    2000    1000       1000       0.5
#> 9     Dominant   Additive 0.2  2    2000    1000       1000       0.5
#> 11    Dominant  Recessive 0.2  2    2000    1000       1000       0.5
#> 13    Dominant   Dominant 0.3  2    2000    1000       1000       0.5
#> 15    Dominant   Additive 0.3  2    2000    1000       1000       0.5
#> 17    Dominant  Recessive 0.3  2    2000    1000       1000       0.5
#> 12   Recessive   Dominant 0.1  2    2000    1000       1000       0.5
#> 31   Recessive   Additive 0.1  2    2000    1000       1000       0.5
#> 51   Recessive  Recessive 0.1  2    2000    1000       1000       0.5
#> 71   Recessive   Dominant 0.2  2    2000    1000       1000       0.5
#> 91   Recessive   Additive 0.2  2    2000    1000       1000       0.5
#> 111  Recessive  Recessive 0.2  2    2000    1000       1000       0.5
#> 131  Recessive   Dominant 0.3  2    2000    1000       1000       0.5
#> 151  Recessive   Additive 0.3  2    2000    1000       1000       0.5
#> 171  Recessive  Recessive 0.3  2    2000    1000       1000       0.5
#> 14    Additive   Dominant 0.1  2    2000    1000       1000       0.5
#> 32    Additive   Additive 0.1  2    2000    1000       1000       0.5
#> 52    Additive  Recessive 0.1  2    2000    1000       1000       0.5
#> 72    Additive   Dominant 0.2  2    2000    1000       1000       0.5
#> 92    Additive   Additive 0.2  2    2000    1000       1000       0.5
#> 112   Additive  Recessive 0.2  2    2000    1000       1000       0.5
#> 132   Additive   Dominant 0.3  2    2000    1000       1000       0.5
#> 152   Additive   Additive 0.3  2    2000    1000       1000       0.5
#> 172   Additive  Recessive 0.3  2    2000    1000       1000       0.5
#> 16         2df   Dominant 0.1  2    2000    1000       1000       0.5
#> 33         2df   Additive 0.1  2    2000    1000       1000       0.5
#> 53         2df  Recessive 0.1  2    2000    1000       1000       0.5
#> 73         2df   Dominant 0.2  2    2000    1000       1000       0.5
#> 93         2df   Additive 0.2  2    2000    1000       1000       0.5
#> 113        2df  Recessive 0.2  2    2000    1000       1000       0.5
#> 133        2df   Dominant 0.3  2    2000    1000       1000       0.5
#> 153        2df   Additive 0.3  2    2000    1000       1000       0.5
#> 173        2df  Recessive 0.3  2    2000    1000       1000       0.5
#>     Power_at_Alpha_0.05
#> 1            0.99997130
#> 3            0.99999117
#> 5            0.06094645
#> 7            0.99999997
#> 9            1.00000000
#> 11           0.12562959
#> 13           0.99999999
#> 15           1.00000000
#> 17           0.26400143
#> 12           0.23736708
#> 31           0.72802319
#> 51           0.32261174
#> 71           0.51913618
#> 91           0.99712300
#> 111          0.84110046
#> 131          0.65907220
#> 151          0.99999669
#> 171          0.99128361
#> 14           0.99994745
#> 32           0.99999535
#> 52           0.09704782
#> 72           0.99999973
#> 92           1.00000000
#> 112          0.39542984
#> 132          0.99999976
#> 152          1.00000000
#> 172          0.83339405
#> 16           0.99987562
#> 33           0.99997633
#> 53           0.24913314
#> 73           0.99999976
#> 93           1.00000000
#> 113          0.75849311
#> 133          0.99999996
#> 153          1.00000000
#> 173          0.97950882
```

"The return object contains information about power for additive, dominant, recessive, and 2df / genotypic tests of association, assuming various true underlying genetic effects (additive, dominant, recessive). "

### Installation instructions

To install genpwr, perform the following steps:

-   Install R version 3.5.1 or higher by following the instructions at <http://www.r-project.org>
-   From the R environment, install and load the "genpwr" package:

``` r
install.packages("genpwr")
```

-   Load the library

``` r
library(genpwr)
```

### Demo

Install the genpwr package as described above.

Run the genpwr demo program

``` r
demo(genpwr_demo)
```
