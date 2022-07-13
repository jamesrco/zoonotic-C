# 03_ConstrainCO2Flux_adjFlux.R
# Created July 5, 2022
# Purpose: Third in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This third script performs the actual adjustment of the Sala et al. benthic
# CO2 flux data and includes a time-integrated calculation for total emissions
# to the atmosphere over a 100 year time period

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

# we are ready at this point to adjust the fluxes in the Sala et al dataset
# using the appropriate (nearest match) benthic sequestration fractions in the
# Siegel et al dataset

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# zoonotic-c/data/derived/output/

load("data/derived/output/coord.matches.RData")

# load the summed emissions fractions for 1 to 100 y (previously generated using 
# gen_fracs_to_constrain_trawlCO2.m)

sum_of_fseq_bottom_1to100years.raw <- read.csv(file = "data/derived/benthic_seqfractions/sum_of_fseq_bottom_1to100years.csv",
                                               header = FALSE)

# load the year-by-year predicted global CO2 sediment remineralization rates from Sala et al.
# (this is a static .csv version of the object "results," generated from the simple
# model beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd)

Sala_et_al_trawlTiming_results.raw <- read.csv(file = "data/derived/sala_et_al_2021_model/Sala_et_al_trawlTiming_results.csv",
                                           header = TRUE) 

colnames(Sala_et_al_trawlTiming_results.raw)[1] <- c("Year")

# define some efflux fractions for each of the timescales (should be 1 - the
# sequestration fraction) and some functions to perform the calculations

effluxFrac.1yr <- 1-as.numeric(unlist(fseq_bottom_1yr.raw))
effluxFrac.5yr <- 1-as.numeric(unlist(fseq_bottom_5yr.raw))
effluxFrac.10yr <- 1-as.numeric(unlist(fseq_bottom_10yr.raw))
effluxFrac.25yr <- 1-as.numeric(unlist(fseq_bottom_25yr.raw))
effluxFrac.50yr <- 1-as.numeric(unlist(fseq_bottom_50yr.raw))
effluxFrac.75yr <- 1-as.numeric(unlist(fseq_bottom_75yr.raw))
effluxFrac.100yr <- 1-as.numeric(unlist(fseq_bottom_100yr.raw))
effluxFrac.1000yr <- 1-as.numeric(unlist(fseq_bottom_1000yr.raw))

# since the underlying Sala et al. dataset is very sparse, can create a subset
# and run our functions on just that subset

ind.nonZeroCO2 <- which(!is.na(Sala_CO2_efflux.df$co2_efflux))
length(ind.nonZeroCO2)/length(Sala_CO2_efflux.df$co2_efflux) # values that aren't NA represent < 1% of the total dataset

constrainFlux.1yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.1yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.5yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.5yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.10yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.10yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.25yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.25yr[coord.matches$Closest[a]]
  # # some code for a progress indicator, if desired
  # thisRecNum <- which(ind.nonZeroCO2 == a)
  # if (thisRecNum %% 1000==0) {
  #   progress <- (thisRecNum/length(ind.nonZeroCO2))*100
  #   print(paste0("Progress: ",progress," %"))
  # }
    return(adjFlux)
}

constrainFlux.50yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.50yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.75yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.75yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.100yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.100yr[coord.matches$Closest[a]]
  return(adjFlux)
}

constrainFlux.1000yr <- function(a){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[a]*effluxFrac.1000yr[coord.matches$Closest[a]]
  return(adjFlux)
}

# run the calculations

# *** these still blew up memory running with the entire dataset; even worse
# when using mclapply ... but much better with a reduced dataset

# 1 year scenario
time0 <- Sys.time()
adjCO2efflux.1yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.1yr))
time1 <- Sys.time()
print(time1 - time0)

# clean up, reformat, save output; clear memory
adjCO2efflux.1yr <- as.data.frame(matrix(data = NA, nrow = length(Sala_CO2_efflux.df$co2_efflux), ncol = 1))
colnames(adjCO2efflux.1yr) <- c("adjCO2_efflux.1yr")
adjCO2efflux.1yr$adjCO2_efflux.1yr[ind.nonZeroCO2] <- adjCO2efflux.1yr_nonZero
save(adjCO2efflux.1yr, file = "data/derived/output/adjCO2efflux.1yr.RData")
gc()

# 25 year scenario
time0 <- Sys.time()
adjCO2efflux.25yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.25yr))
time1 <- Sys.time()
print(time1 - time0)

# clean up, reformat, save output; clear memory
adjCO2efflux.25yr <- as.data.frame(matrix(data = NA, nrow = length(Sala_CO2_efflux.df$co2_efflux), ncol = 1))
colnames(adjCO2efflux.25yr) <- c("adjCO2_efflux.25yr")
adjCO2efflux.25yr$adjCO2_efflux.25yr[ind.nonZeroCO2] <- adjCO2efflux.25yr_nonZero
save(adjCO2efflux.25yr, file = "data/derived/output/adjCO2efflux.25yr.RData")
gc()

# 50 year scenario
time0 <- Sys.time()
adjCO2efflux.50yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.50yr))
time1 <- Sys.time()
print(time1 - time0)

# clean up, reformat, save output; clear memory
adjCO2efflux.50yr <- as.data.frame(matrix(data = NA, nrow = length(Sala_CO2_efflux.df$co2_efflux), ncol = 1))
colnames(adjCO2efflux.50yr) <- c("adjCO2_efflux.50yr")
adjCO2efflux.50yr$adjCO2_efflux.50yr[ind.nonZeroCO2] <- adjCO2efflux.50yr_nonZero
save(adjCO2efflux.50yr, file = "data/derived/output/adjCO2efflux.50yr.RData")
gc()

# 100 year scenario
time0 <- Sys.time()
adjCO2efflux.100yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.100yr))
time1 <- Sys.time()
print(time1 - time0)

# clean up, reformat, save output; clear memory
adjCO2efflux.100yr <- as.data.frame(matrix(data = NA, nrow = length(Sala_CO2_efflux.df$co2_efflux), ncol = 1))
colnames(adjCO2efflux.100yr) <- c("adjCO2_efflux.100yr")
adjCO2efflux.100yr$adjCO2_efflux.100yr[ind.nonZeroCO2] <- adjCO2efflux.100yr_nonZero
save(adjCO2efflux.100yr, file = "data/derived/output/adjCO2efflux.100yr.RData")
gc()

# 1000 year scenario
time0 <- Sys.time()
adjCO2efflux.1000yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.1000yr))
time1 <- Sys.time()
print(time1 - time0)

# clean up, reformat, save output; clear memory
adjCO2efflux.1000yr <- as.data.frame(matrix(data = NA, nrow = length(Sala_CO2_efflux.df$co2_efflux), ncol = 1))
colnames(adjCO2efflux.1000yr) <- c("adjCO2_efflux.1000yr")
adjCO2efflux.1000yr$adjCO2_efflux.1000yr[ind.nonZeroCO2] <- adjCO2efflux.1000yr_nonZero
save(adjCO2efflux.1000yr, file = "data/derived/output/adjCO2efflux.1000yr.RData")
gc()

# 5, 10, 75 year scenarios; for prediction curve
time0 <- Sys.time()
adjCO2efflux.5yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.5yr))
time1 <- Sys.time()
print(time1 - time0)
sum(adjCO2efflux.5yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9)

time0 <- Sys.time()
adjCO2efflux.10yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.10yr))
time1 <- Sys.time()
print(time1 - time0)
sum(adjCO2efflux.10yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9)

time0 <- Sys.time()
adjCO2efflux.75yr_nonZero <- unlist(lapply(ind.nonZeroCO2, constrainFlux.75yr))
time1 <- Sys.time()
print(time1 - time0)
sum(adjCO2efflux.75yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9)

# get some basic statistics; with unadjusted Sala et al data for comparison

sum(adjCO2efflux.1yr, na.rm=T)
# [1] 1859877
sum(adjCO2efflux.25yr, na.rm=T)
# [1] 26483089
sum(adjCO2efflux.50yr, na.rm=T)
# [1] 36595847
sum(adjCO2efflux.100yr, na.rm=T)
# [1] 52273243
sum(adjCO2efflux.1000yr, na.rm=T)
# [1] 559481566
sum(Sala_CO2_efflux.df$co2_efflux, na.rm=T) # in Mg per km^2
# [1] 1691236544 # this is the same as cellStats(Sala_CO2_efflux.raw, stat="sum")

# but this value (1.69 Pg CO2) is not the same as the 1.47 Pg CO2 reported in the 
# Sala et al. paper

# email to jmayorga@bren.ucsb.edu (replied within an hour -- thanks Juan) yielded
# following advice: "the fluxes in the Geotiff are per km2, so please make sure
# to multiply times the pixel’s area before summing up."

# so we can adjust accordingly

SalaModel_cell_area <- 934.4789^2/1000000 # from https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/data_prep/update_bottom_trawling_impact.Rmd

sum(adjCO2efflux.1yr*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 0.001624139
sum(adjCO2efflux.25yr*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 0.02312638
sum(adjCO2efflux.50yr*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 0.03195735
sum(adjCO2efflux.100yr*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 0.04564765
sum(adjCO2efflux.1000yr*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 0.4885677
sum(Sala_CO2_efflux.df$co2_efflux*SalaModel_cell_area, na.rm=T)*(1/10^9) # in Pg CO2
# [1] 1.476874 # *** now this matches the value in the Sala et al paper
# also as check: cellStats(Sala_CO2_efflux.raw*SalaModel_cell_area, stat="sum")*(1/10^9)
# [1] 1.476874

# now, save these values to file

adjCO2efflux_PgCO2_yr <- c(sum(adjCO2efflux.1yr*SalaModel_cell_area, na.rm=T)*(1/10^9),
                           sum(adjCO2efflux.5yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9),
                           sum(adjCO2efflux.10yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9),
                           sum(adjCO2efflux.25yr*SalaModel_cell_area, na.rm=T)*(1/10^9),
                           sum(adjCO2efflux.50yr*SalaModel_cell_area, na.rm=T)*(1/10^9),
                           sum(adjCO2efflux.75yr_nonZero, na.rm=T)*SalaModel_cell_area*(1/10^9),
                           sum(adjCO2efflux.100yr*SalaModel_cell_area, na.rm=T)*(1/10^9),
                           sum(adjCO2efflux.1000yr*SalaModel_cell_area, na.rm=T)*(1/10^9),
                           sum(Sala_CO2_efflux.df$co2_efflux*SalaModel_cell_area, na.rm=T)*(1/10^9))

names(adjCO2efflux_PgCO2_yr) <- c("adjCO2efflux.1yr_PgCO2_per_y",
                                  "adjCO2efflux.5yr_PgCO2_per_y",
                                  "adjCO2efflux.10yr_PgCO2_per_y",
                                  "adjCO2efflux.25yr_PgCO2_per_y",
                                  "adjCO2efflux.50yr_PgCO2_per_y",
                                  "adjCO2efflux.75yr_PgCO2_per_y",
                                  "adjCO2efflux.100yr_PgCO2_per_y",
                                  "adjCO2efflux.1000yr_PgCO2_per_y",
                                  "CO2efflux_PgCO2_per_y.unadjusted")

write.csv(adjCO2efflux_PgCO2_yr, file = "data/derived/output/adjCO2efflux_PgCO2_yr.csv",
          row.names = TRUE)

# now we can make a time-integrated estimate of total emissions to the
# atmosphere over a 100 year time period, using the same trawling scenario
# as in the original Sala et al paper

# first, need to make a prediction curve from the model results we have

adjCO2efflux_PgCO2_yr.for_prediction <- as.data.frame(matrix(data = NA, 
                                                             nrow = 7,
                                                             ncol = 3))
colnames(adjCO2efflux_PgCO2_yr.for_prediction) = c("Year",
             "PgCO2_per_year_to_atmos",
             "Sum_emitFrac")

adjCO2efflux_PgCO2_yr.for_prediction[,2] <- adjCO2efflux_PgCO2_yr[1:7]
adjCO2efflux_PgCO2_yr.for_prediction[,1] <- c(1,5,10,25,50,75,100)
adjCO2efflux_PgCO2_yr.for_prediction[,3] <- sum_of_fseq_bottom_1to100years.raw[sum_of_fseq_bottom_1to100years.raw[,1] %in% adjCO2efflux_PgCO2_yr.for_prediction[,1],2]

model <- lm(PgCO2_per_year_to_atmos ~ poly(Sum_emitFrac,2, raw = TRUE), data = adjCO2efflux_PgCO2_yr.for_prediction)

# > summary(model)
# 
# Call:
#   lm(formula = PgCO2_per_year_to_atmos ~ poly(Sum_emitFrac, 2, 
#                                               raw = TRUE), data = adjCO2efflux_PgCO2_yr.for_prediction)
# 
# Residuals:
#   1          2          3          4          5          6          7 
# -2.372e-04  1.114e-04  6.728e-05  3.818e-04  7.092e-06 -7.310e-04  4.006e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                        -2.133e-03  5.112e-04  -4.172  0.01400 *  
#   poly(Sum_emitFrac, 2, raw = TRUE)1  6.214e-05  2.251e-06  27.603 1.02e-05 ***
#   poly(Sum_emitFrac, 2, raw = TRUE)2 -1.662e-08  1.935e-09  -8.590  0.00101 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.000478 on 4 degrees of freedom
# Multiple R-squared:  0.9994,	Adjusted R-squared:  0.9991 
# F-statistic:  3347 on 2 and 4 DF,  p-value: 3.567e-07

# make some predictions for 1 to 100 y using the model we just fit

predicted.PgCO2_per_year_to_atmos <- as.data.frame(matrix(data = NA, 
                                                          nrow = 100,
                                                          ncol = 2))
colnames(predicted.PgCO2_per_year_to_atmos) = c("Year",
                                                "PgCO2_per_year_to_atmos.predicted")
predicted.PgCO2_per_year_to_atmos[,1] <- c(1:100)
predicted.PgCO2_per_year_to_atmos[,2] <- predict(model2,newdata=data.frame(Sum_emitFrac = sum_of_fseq_bottom_1to100years.raw[,2]))

# now, can make a time-integrated calculation over 100 years

# assumes:
# 1. annual CO2 flux from sediments follows the trend posited by Sala et al. 2021,
# contained in the object "results" that is generated beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd
# *** this is where the assertion in the Sala et al. paper that sediment emissions
# in the first year are 1.47 Pg CO2, declining after 10 y to 0.58 Pg CO2, comes from 
# 2. the fraction of the sediment CO2 emissions in year n that will have reached the
# atmosphere 100-n years later is described by the benthic emissions fraction
# for the nth year based on Siegel et al. 2021 (1 - fseq_bottom)

sum(rev(predicted.PgCO2_per_year_to_atmos[,2])*(Sala_et_al_trawlTiming_results.raw$C_remin[1:100]/
                                                  Sala_et_al_trawlTiming_results.raw$C_remin[1]))

# after 100 y of continuous trawling, a cumulative 1.27 Pg CO2 will have reached the atmosphere 

# # send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.
# 
# system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))