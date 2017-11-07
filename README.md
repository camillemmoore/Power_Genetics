## Overview
This is Camille's GitHub Repository for power and sample size calculations for genetic association studies.  Currently, the genpwr package allows for power and sample size calcualtions in case/control studies under genetic model mis-specification.  

Power and/or sample size can be calculted for the following logisitic regression models:
- Additive Genetic Effect
- Dominant Genetic Effect
- Recessive Genetic Effect
- Unspecified Genetic Effect (2 degree of freedom test - most flexible)

The following true underlying genetic models can be assumed:
- Additive
- Dominant
- Recessive

Currently, covariates and gene x environment interactions cannot be included in the models. This package will be extended to continuous outcomes.  

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
install_github("camillemmoore/Power_Genetics", subdir="genpwr"))
```
