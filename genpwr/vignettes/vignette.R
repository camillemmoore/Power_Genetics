## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
knitr::opts_chunk$set(fig.width=8, fig.height=5) 
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(dplyr)
library(ggplot2)
set.seed(1014)

## ---- eval=F-------------------------------------------------------------
#  install.packages("devtools")

## ------------------------------------------------------------------------
library(devtools)

## ---- eval = F-----------------------------------------------------------
#  install_github("camillemmoore/Power_Genetics", subdir="genpwr")

## ------------------------------------------------------------------------
library(genpwr)

