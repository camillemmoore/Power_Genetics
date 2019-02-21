## Overview
This is the GitHub Repository for genpwr, an R package which aids in power and sample size calculations for genetic association studies.  The genpwr package allows for power and sample size calculations for case/control studies, as well as studies of continuous phenotypes, under genetic model mis-specification.  

Power and/or sample size can be calculted for the following logisitic and linear regression models:
- Additive Genetic Effect
- Dominant Genetic Effect
- Recessive Genetic Effect
- Unspecified Genetic Effect (2 degree of freedom test - most flexible)

The following true underlying genetic models can be assumed:
- Additive
- Dominant
- Recessive

Gene x environment interactions with a binary (ex. smoking yes/no) or continuous environmental factor can be included in the models.   

## Install genpwr

To install the package, first, you need to install the devtools package. You can do this from CRAN. Invoke R and then type:

```
install.packages("devtools")
```

Next, load the devtools package:

```
library(devtools)
```

Finally, install the package.

```
install_github("camillemmoore/Power_Genetics", subdir="genpwr")
```

## Vignettes and demos
To see examples of how to use the genpwr package, see our vignette:

```
browseVignettes('genpwr')
```

There is also a demo to reproduce the study designs in our upcoming paper:

```
demo(genpwr_demo)
```
