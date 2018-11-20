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

## ------------------------------------------------------------------------
pw <- power.calc(N=400, Case.Rate=0.4, k=NULL,
   MAF=seq(0.05, 0.1, 0.01), OR=c(3),Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
head(pw)

## ------------------------------------------------------------------------
power.plot(pw)

## ------------------------------------------------------------------------
ss <- ss.calc(OR=4, Case.Rate=0.4, k=NULL,
   MAF=seq(0.3, 0.4, 0.02), power=0.8, Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
ss.plot(ss)

## ------------------------------------------------------------------------
or <- odds_ratio_function(N=1000, Case.Rate=0.4, k=NULL,
   MAF=seq(0.30, 0.45, 0.02), power=0.4, Alpha=0.05,
   True.Model="All", Test.Model="All")

## ------------------------------------------------------------------------
or.plot(or)

## ------------------------------------------------------------------------
pec <- power_envir.calc(N=500, Case.Rate=0.3, MAF=seq(0.2,0.4,0.02), OR_G=3, 
                 OR_E=3.5, OR_GE=4, P_e = 0.4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
power.plot(pec)

## ------------------------------------------------------------------------
sse <- ss_envir.calc(power=0.5, Case.Rate=0.4, MAF=seq(0.3,0.4,0.05), OR_G=1.5, 
                 OR_E=1.3, OR_GE=1.2, P_e = 0.4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
ss.plot(sse)

