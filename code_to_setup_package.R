#R Script to Create a Power_Genetics Package

#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("H:/CM_Repositories/Power_Genetics/")
#create("genpwr")

setwd("./genpwr")
document()

setwd("..")
install("genpwr")

library(genpwr)
