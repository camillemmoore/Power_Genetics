## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
knitr::opts_chunk$set(fig.width=8, fig.height=5) 
options(tibble.print_min = 4L, tibble.print_max = 4L)
#library(dplyr)
#library(ggplot2)
set.seed(1014)

## ---- eval=F-------------------------------------------------------------
#  install.packages("devtools")

## ---- eval = F-----------------------------------------------------------
#  library(devtools)

## ---- eval = F-----------------------------------------------------------
#  install_github("camillemmoore/Power_Genetics", subdir="genpwr")

## ------------------------------------------------------------------------
library(genpwr)

## ------------------------------------------------------------------------
pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
   N=400, Case.Rate=0.2, k=NULL,
   MAF=seq(0.18, 0.25, 0.01), OR=c(3),Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
head(pw)

