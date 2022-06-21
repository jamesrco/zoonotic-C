# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# libraries

library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(rgdal)

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
fseq_bottom_lat_degNr.raw <- read.csv("data/derived/benthic_seqfractions/lat_degN.csv",
                                 header = FALSE)
fseq_bottom_long_degE.raw <- read.csv("data/derived/benthic_seqfractions/long_degE.csv",
                                 header = FALSE)

# load Sala et al. pCO2 data, as GeoTIFF
# with some help from https://datacarpentry.org/r-raster-vector-geospatial/01-raster-structure/

Sala_bottomtrawl_Ia.raw <- 
  raster("data/raw/sala_et_al_2021/bottom_trawling_Ia.tif")

Sala_carbon_ranking.raw <- 
  raster("data/raw/sala_et_al_2021/carbon_ranking.tif")

Sala_CO2_efflux.raw <- 
  raster("data/raw/sala_et_al_2021/co2_efflux.tif")