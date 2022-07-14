# 03b_ConstrainCO2Flux_adjFlux_extended.R
# Created July 12, 2022
# Purpose: Alternate, more robust version for third in series of scripts used
# to constrain the estimate of benthic CO2 flux in Sala  et al. 2021 using the
# sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This alternate third script performs the actual adjustment of the Sala et al
# benthic CO2 flux data and includes a time-integrated calculation for total
# emissions to the atmosphere over a 100 year time period. This script is
# different than the primary version of the third script (03a_ConstrainCO2Flux_adjFlux.R)
# in that calculations for each year are made explictly, rather than estimated
# using a prediction curve constructed from a subset of values

# *** Assumes user has already run 01_ConstrainCO2Flux_IO.R (in current session)
# and that the object "coord.matches.RData" generated using
# 02_ConstrainCO2Flux_coordMatch.R (likely via AWS) is present in
# zoonotic-c/data/derived/output/

# set the working directory; create directory for output

setwd("~/Code/zoonotic-C/") # for my laptop
# setwd("~/zoonotic-C/") # for AWS

# libraries

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file

# we are ready at this point to adjust the fluxes in the Sala et al dataset
# using the appropriate (nearest match) benthic sequestration fractions in the
# Siegel et al dataset

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# zoonotic-c/data/derived/output/

load("data/derived/output/coord.matches.RData")

# load the benthic sequestration fractions for entire Siegel et al. model domain,
# from 1-200 y, then 200-1000 y in 100 y increments
# also load the years

fseq_bottom_multYears.raw <- readMat("data/derived/benthic_seqfractions/fseq_bottom_multyears.mat")
fseq_bottom.multyears <- fseq_bottom_multYears.raw$fseq.bottom.multyears # clean up a bit
fseq_bottom.multyears[fseq_bottom.multyears>=1] <- 1

seqFracYears.raw <- read.csv(file = "data/derived/benthic_seqfractions/benthic_years.csv",
                             header = FALSE)

# load the year-by-year predicted global CO2 sediment remineralization rates from Sala et al.
# (this is a static .csv version of the object "results," generated from the simple
# model beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd)

Sala_et_al_trawlTiming_results.raw <- read.csv(file = "data/derived/sala_et_al_2021_model/Sala_et_al_trawlTiming_results.csv",
                                           header = TRUE) 

colnames(Sala_et_al_trawlTiming_results.raw)[1] <- c("Year")

# since the underlying Sala et al. dataset is very sparse, can create a subset
# and run our functions on just that subset

ind.nonZeroCO2 <- which(!is.na(Sala_CO2_efflux.df$co2_efflux))
length(ind.nonZeroCO2)/length(Sala_CO2_efflux.df$co2_efflux) # values that aren't NA represent < 1% of the total dataset

# need to define cell area per email from jmayorga@bren.ucsb.edu:
# "the fluxes in the Geotiff are per km2, so please make sure to multiply times
# the pixel’s area before summing up."

SalaModel_cell_area <- 934.4789^2/1000000 # from https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/data_prep/update_bottom_trawling_impact.Rmd

# generalized version of functions

genEffluxFracs <- function(year){
  yearInd <- which(seqFracYears.raw==year)
  fseq_bottom.thisyear <- fseq_bottom.multyears[,,yearInd]
  effluxFrac_bottom.thisyear <- 1-as.numeric(unlist(fseq_bottom.thisyear))
  return(effluxFrac_bottom.thisyear)
}

constrainFlux <- function(dataIndex, effluxFracs){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[dataIndex]*effluxFracs[coord.matches$Closest[dataIndex]]
  return(adjFlux)
}

# run the calculations

# *** these still blew up memory running with the entire dataset; even worse
# when using mclapply ... but much better with a reduced dataset

# set up structure to hold results

predicted.PgCO2_per_year_to_atmos <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 2))
colnames(predicted.PgCO2_per_year_to_atmos) = c("Year",
                                                "PgCO2_per_year_to_atmos")
predicted.PgCO2_per_year_to_atmos[,1] <- unlist(seqFracYears.raw)

# iterate

for (i in 200:nrow(predicted.PgCO2_per_year_to_atmos)) {

  print(predicted.PgCO2_per_year_to_atmos[i,1])
  
  time0 <- Sys.time()
  
  thisYear <- predicted.PgCO2_per_year_to_atmos[i,1]
  EffluxFracs.thisyear <- genEffluxFracs(thisYear)
  adjCO2efflux.thisyear <- unlist(lapply(ind.nonZeroCO2, constrainFlux, EffluxFracs.thisyear))
  predicted.PgCO2_per_year_to_atmos[i,2] <- sum(adjCO2efflux.thisyear*SalaModel_cell_area, na.rm=T)*(1/10^9)
  
  time1 <- Sys.time()
  print(time1 - time0)
  
}

# save

write.csv(predicted.PgCO2_per_year_to_atmos, file = "data/derived/output/adjCO2efflux_PgCO2_yr.csv",
          row.names = FALSE)

# now we can make time-integrated estimates of total emissions to the atmosphere
# can pick our base year, using 100 y as example

# assumes:
# 1. annual CO2 flux from sediments follows the trend posited by Sala et al. 2021,
# contained in the object "results" that is generated beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd
# *** this is where the assertion in the Sala et al. paper that sediment emissions
# in the first year are 1.47 Pg CO2, declining after 10 y to 0.58 Pg CO2, comes from 
# 2. the fraction of the sediment CO2 emissions in year n that will have reached the
# atmosphere 100-n years later is described by the benthic emissions fraction
# for the nth year based on Siegel et al. 2021 (1 - fseq_bottom)

# first, over a 100 year time period, using the same trawling scenario
# as in the original Sala et al paper

timeInt.PgCO2_to_atmos.100y <- sum(rev(predicted.PgCO2_per_year_to_atmos[1:100,2])*(Sala_et_al_trawlTiming_results.raw$C_remin[1:100]/
                                                  Sala_et_al_trawlTiming_results.raw$C_remin[1]))

# [1] 1.265947

# after 100 y of continuous trawling, a cumulative 1.27 Pg CO2 will have reached the atmosphere 

# 200 year time period

timeInt.PgCO2_to_atmos.200y <- sum(rev(predicted.PgCO2_per_year_to_atmos[1:200,2])*(Sala_et_al_trawlTiming_results.raw$C_remin[1:200]/
                                                                                 Sala_et_al_trawlTiming_results.raw$C_remin[1]))

# [1] 3.945729

# after 100 y of continuous trawling, a cumulative 3.95 Pg CO2 will have reached the atmosphere 

# # send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.
# 
# system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))