# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# libraries

library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(rgdal)
library(sp)
library(data.table)

# library(ncdf4) # not needed anymore

# load datasets

# # Attempt to load Siegel et al. "sequestration fractions" from NetCDF
# 
# OCIM2_fseq_48L <- nc_open("data/siegel_et_al_2021_v2/fseq_OCIM2_48L.nc", verbose=TRUE) # currently returns an error

# instead we'll use a modified version of the MATLAB script provided by Siegel et al. to get what we need: see "gen_fracs_to_constrain_trawlCO2.m" which should be in this same directory

# take the necessary detour into MATLAB at this point, if the output hasn't been generated; then return to R

# the output from gen_fracs_to_constrain_trawlCO2.m -- 25, 50, and 100 year sequestration fractions for ocean bottom depths, plus the necessary metadata -- should now be in several .csv files found in data/derived/benthic_seqfractions

# load these output files

fseq_bottom_25yr.raw <- read.csv("data/derived/benthic_seqfractions/fseq_bottom_25yr.csv",
                            header = FALSE)
fseq_bottom_50yr.raw <- read.csv("data/derived/benthic_seqfractions/fseq_bottom_50yr.csv",
                                 header = FALSE)
fseq_bottom_100yr.raw <- read.csv("data/derived/benthic_seqfractions/fseq_bottom_100yr.csv",
                                 header = FALSE)
fseq_bottom_depth_m.raw <- read.csv("data/derived/benthic_seqfractions/bottom_depth_m.csv",
                                 header = FALSE)
fseq_bottom_lat_degN.raw <- read.csv("data/derived/benthic_seqfractions/lat_degN.csv",
                                 header = FALSE)
fseq_bottom_long_degE.raw <- read.csv("data/derived/benthic_seqfractions/long_degE.csv",
                                 header = FALSE)

# load Sala et al. pCO2 data, as GeoTIFF, then convert to data frames
# with some help from https://datacarpentry.org/r-raster-vector-geospatial/01-raster-structure/ and
# https://www.neonscience.org/resources/learning-hub/tutorials/raster-data-r

# good instructions on enabling multithreading on Mac (for data.table) here:
# https://github.com/Rdatatable/data.table/wiki/Installation

Sala_bottomtrawl_Ia.raw <- 
  raster("data/raw/sala_et_al_2021/bottom_trawling_Ia.tif")
Sala_bottomtrawl_Ia.df <- as.data.frame(Sala_bottomtrawl_Ia.raw, xy = TRUE)

Sala_carbon_ranking.raw <- 
  raster("data/raw/sala_et_al_2021/carbon_ranking.tif")
Sala_carbon_ranking.df <- as.data.frame(Sala_carbon_ranking.raw, xy = TRUE)

Sala_CO2_efflux.raw <- 
  raster("data/raw/sala_et_al_2021/co2_efflux.tif")
# Sala_CO2_efflux.df <- as.data.frame(Sala_CO2_efflux.raw, xy = TRUE) 
# # note that creating a data frame from this last object will require a lot of memory ... I had to up 
# # the R_MAX_VSIZE variable in .Renviron according to the directions here:
# # https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

# now, can find most appropriate (nearest) sequestration fraction for each point in the CO2 flux dataset

# create data frames of the coordinates of the points in the two different datasets

# CO2 flux data
Sala_CO2_efflux.coords <- xyFromCell(Sala_CO2_efflux.raw,c(1:length(Sala_CO2_efflux.raw))) # easy, can just use the xyFromCell function in raster package

# sequestration fractions
Siegel_fseq.coords <- data.frame(as.numeric(rep(fseq_bottom_long_degE.raw[1,],91)),rep(fseq_bottom_lat_degN.raw[,1],180))
colnames(Siegel_fseq.coords) <- c("x","y")
# convert to eastings & westings rather than just eastings, for compatibility 
Siegel_fseq.coords$x[Siegel_fseq.coords$x>180] <- Siegel_fseq.coords$x[Siegel_fseq.coords$x>180]-360

# now have (more or less) apples to oranges; find best match for each point in the CO2 flux dataset

# convert to data tables
Sala_CO2_efflux.coords.dt <- data.table(Sala_CO2_efflux.coords)
Siegel_fseq.coords.dt <- data.table(Siegel_fseq.coords)

# define function to find nearest match 

dist <- function(a, b){
  dt <- data.table((Siegel_fseq.coords.dt$x-a)^2+(Siegel_fseq.coords.dt$y-b)^2)
  return(which.min(dt$V1))}

coord.matches <- Sala_CO2_efflux.coords.dt[, j = list(Closest =  dist(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt)]

dist.alt <- function(a, b){
  dt <- data.table((Sala_CO2_efflux.coords.dt$x-a)^2+(Sala_CO2_efflux.coords.dt$y-b)^2)
  return(which.min(dt$V1))}

coord.matches.alt <- Siegel_fseq.coords.dt[, j = list(Closest =  dist.alt(x, y)), by = 1:nrow(Siegel_fseq.coords.dt)]

