The genpwr Package
=========================

The informativeDropout package for R (>3.5.1) performs power and sample size calculations     
    for case/control genetic association studies and allows for   
    model mis-specification.

### Installation instructions

The results in the above manuscript were produced using R version 3.5.1. To reproduce the results,
perform the following steps:

* Install R version 3.5.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
install.packages("devtools")
library(devtools)
```
* Install the "informativeDropout" package directly from Github.com
```R
install_github("camillemmoore/Power_Genetics", subdir="genpwr")
```
* Load the library
```R
library(genpwr)
```

### Demo

Install the genpwr package as described above.

Run the genwr demo program

```R
demo(genpwr_demo)
```





