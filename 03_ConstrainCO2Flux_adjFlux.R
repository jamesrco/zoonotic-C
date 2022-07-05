# 03_ConstrainCO2Flux_adjFlux.R
# Created July 5, 2022
# Purpose: Third in series of scripts used to constrain the crude estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This third script performs the actual adjustment of the Sala et al. benthic
# CO2 flux data

# *** Assumes user has already run 01_ConstrainCO2Flux_IO.R (in current session)
# and that the object "coord.matches.RData" generated using
# 02_ConstrainCO2Flux_coordMatch.R (likely via AWS) is present in
# zoonotic-c/data/derived/output/

# set the working directory; create directory for output

# setwd("~/Code/zoonotic-C/") # for my laptop
setwd("~/zoonotic-C/") # for AWS

# libraries

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed

# we are ready at this point to adjust the fluxes in the Sala et al dataset
# using the appropriate (nearest match) benthic sequestration fractions in the
# Siegel et al dataset

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# zoonotic-c/data/derived/output/

load("data/derived/output/coord.matches.RData")

# define some efflux fractions for each of the timescales (should be 1 - the
# sequestration fraction) and some functions to perform the calculations

effluxFrac.25yr <- 1-as.numeric(unlist(fseq_bottom_25yr.raw))
effluxFrac.50yr <- 1-as.numeric(unlist(fseq_bottom_50yr.raw))
effluxFrac.100yr <- 1-as.numeric(unlist(fseq_bottom_100yr.raw))

constrainFlux.25yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.25yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.50yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.25yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.100yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.25yr[coord.matches$Closest[a]]
  return(adjFlux)
}

# run the calculations; set up for parallel processing (with some benchmarking)

time0 <- Sys.time()
adjCO2efflux.25yr <- unlist(mclapply(1:length(Sala_CO2_efflux.df$co2_efflux), constrainFlux.25yr, mc.cores = 32))
time1 <- Sys.time()
print(time1 - time0)

# save output
save(adjCO2efflux.25yr, file = "data/derived/output/adjCO2efflux.25yr.RData")

time0 <- Sys.time()
adjCO2efflux.50yr <- unlist(mclapply(1:length(Sala_CO2_efflux.df$co2_efflux), constrainFlux.50yr, mc.cores = 32))
time1 <- Sys.time()
print(time1 - time0)

# save output
save(adjCO2efflux.50yr, file = "data/derived/output/adjCO2efflux.50yr.RData")

time0 <- Sys.time()
adjCO2efflux.100yr <- unlist(mclapply(1:length(Sala_CO2_efflux.df$co2_efflux), constrainFlux.100yr, mc.cores = 32))
time1 <- Sys.time()
print(time1 - time0)

# save output
save(adjCO2efflux.100yr, file = "data/derived/output/adjCO2efflux.100yr.RData")

# get some basic statistics; with unadjusted data for comparison

sum(adjCO2efflux.25yr, na.rm=T)
sum(adjCO2efflux.50yr, na.rm=T)
sum(adjCO2efflux.100yr, na.rm=T)
sum(Sala_CO2_efflux.df$co2_efflux, na.rm=T)

# send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.

system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))
