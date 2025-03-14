

# Demo


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(tidyverse)
library(egg)
library(goFlux)
library(devtools)
library(zoo)
library(pbapply)

source(file = "R/clickflux.R")
source(file = "R/get_dxdy.R")

# Loading data
mydata_all <- NULL
fs <- list.files(path = "data/",pattern = ".RData", full.names = T)
for(f in fs){
  load(file = f)
  mydata$Etime <- as.numeric(mydata$Etime)
  mydata_all <- rbind(mydata_all, mydata)
  rm(mydata)
}

# Loading auxfile table
myauxfile = read.csv("data/myauxfile.csv")
myauxfile$start.time <- as.POSIXct(myauxfile$start.time, tz = 'UTC', format="%d/%m/%Y %H:%M")


# plot incubations overview
p <- plot.incubations(mydata_all)
print(p)
# to save these plots in a dedicated path, do
# gg_save_pdf(list = p, path = , filename = "myfilename.pdf")


# manual inspection of CO2 data and flux calculation
CO2_manID <- clickflux(dataframe = mydata_all, myauxfile = myauxfile,
                          shoulder = 0,
                          gastype = "CO2dry_ppm",
                          plot.lim = c(200,1000),
                       fluxSeparation = F, displayPlots = T)

# take a look at CH4 data and identify possible ebullition events
p <- plot.fluxSeparation(dataframe = mydata_all, gastype = "CH4dry_ppb", kstar = 0.4)
print(p)

CH4_manID <- clickflux(mydata_all = mydata_all, myauxfile = myauxfile,
                       shoulder = 0,
                       gastype = "CH4dry_ppb",
                       plot.lim = c(1800,max(mydata_all$CH4dry_ppb)),
                       fluxSeparation = T, displayPlots = T)




