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
pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
   N=400, Case.Rate=0.2, k=NULL,
   MAF=seq(0.18, 0.25, 0.01), OR=c(3),Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
head(pw)

## ------------------------------------------------------------------------
power.plot(pw)

## ------------------------------------------------------------------------
pw <- genpwr.calc(calc = "power", model = "linear",
   N=40, sd_y=4, k=NULL,
   MAF=seq(0.18, 0.25, 0.01), ES=c(3),Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
ss <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = NULL,
   OR=4, Case.Rate=0.4, k=NULL,
   MAF=seq(0.3, 0.4, 0.02), Power=0.8, Alpha=0.05,
   True.Model=c("Dominant", "Recessive", "Additive"), 
   Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

## ------------------------------------------------------------------------
ss.plot(ss)

## ------------------------------------------------------------------------
or <- genpwr.calc(calc = "es", model = "logistic", ge.interaction = NULL,
   N=1000, Case.Rate=0.4, k=NULL,
   MAF=seq(0.30, 0.45, 0.02), Power=0.4, Alpha=0.05,
   True.Model="All", Test.Model="All")

## ------------------------------------------------------------------------
or.plot(or)

## ------------------------------------------------------------------------
pec <- genpwr.calc(calc = "power", model = "logistic", 
                 ge.interaction = "binary",
                 N=500, Case.Rate=0.3, MAF=seq(0.2,0.4,0.02), OR_G=3, 
                 OR_E=3.5, OR_GE=4, P_e = 0.4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
power.plot(pec)

## ------------------------------------------------------------------------
pec <- genpwr.calc(calc = "power", model = "linear",
                 ge.interaction = "binary",
                 N=50, sd_y=3, MAF=seq(0.2,0.34,0.02), ES_G=3, 
                 ES_E=1.5, ES_GE=2, P_e = 0.4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
pec <- genpwr.calc(calc = "power", model = "logistic", 
                 ge.interaction = "continuous",
                 N=500, Case.Rate=0.3, MAF=seq(0.25,0.4,0.02), OR_G=3, 
                 OR_E=3.5, OR_GE=4, sd_e = 4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
sse <- genpwr.calc(calc = "ss", model = "logistic", ge.interaction = "binary",
                 Power=0.5, Case.Rate=0.4, MAF=seq(0.3,0.4,0.05), OR_G=1.5, 
                 OR_E=1.3, OR_GE=1.2, P_e = 0.4, 
                 Alpha=0.05, True.Model='All', Test.Model='All')

## ------------------------------------------------------------------------
ss.plot(sse)

