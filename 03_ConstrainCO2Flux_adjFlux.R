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

setwd("~/Code/zoonotic-C/") # for my laptop

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

# load in the file containing the coordinate matches, generated in previous
# script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# zoonotic-c/data/derived/output/

load("data/derived/output/coord.matches.RData")

# define functions to appropriately adjust the Sala et al data

constrainFlux.25yr <- function(rawFlux, closestMatch){
  adjFlux <- rawFlux*(1-fseq_bottom_25yr.raw[closestMatch])
  return adjFlux
}
