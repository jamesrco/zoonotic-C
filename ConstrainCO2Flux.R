# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# set the working directory; create directory for output

setwd("~/zoonotic-c")
dir.create("output")

# libraries

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed

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
# https://github.com/Rdatatable/data.table/wiki/Installation and (even more helpful)
# here: https://firas.io/post/data.table_openmp/
# *** removing & recompiling data.table from source is critical if you already have it 
# installed

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
Siegel_fseq.coords.df <- data.frame(as.numeric(rep(fseq_bottom_long_degE.raw[1,],91)),rep(fseq_bottom_lat_degN.raw[,1],180))
colnames(Siegel_fseq.coords.df) <- c("x","y")
# convert to eastvings & westings rather than just eastings, for compatibility 
Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180] <- Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180]-360

# now have (more or less) apples to oranges; find best match for each point in the CO2 flux dataset

# ******************************************************************************
# one approach: using just data.table
# ******************************************************************************

# *** seemed to take too long on my laptop and didn't lend itself to obvious parallelization
# but, given the error with approach #2, below, trying this with a 32 vCPU EC2
# instance on AWS (64 GB RAM plus a 35 GB disk swap)

# first, set number of threads (need to change depending on cores or vCPUS)

setDTthreads(32) # shouldn't allow you to exceed actual # of cores or vCPUS,
                # at least on linux

# convert to data tables; add decimal where necessary
Sala_CO2_efflux.coords.dt <- data.table(Sala_CO2_efflux.coords/10^5)
Siegel_fseq.coords.dt <- data.table(Siegel_fseq.coords.df)

# define function to find nearest match
# adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame

dist1 <- function(a, b){
  dt <- data.table((Siegel_fseq.coords.dt$x-a)^2+(Siegel_fseq.coords.dt$y-b)^2)
  return(which.min(dt$V1))}

# find matches

# # test with a subset first; with benchmarking
# 
# Sala_CO2_efflux.coords.dt.sub <- Sala_CO2_efflux.coords.dt[1:10000,]
# 
# # find matches
# 
# time0 <- Sys.time()
# 
# coord.matches.test1 <- Sala_CO2_efflux.coords.dt.sub[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt.sub)]
# 
# time1 <- Sys.time()
# print(time1 - time0)

# now with the whole enchilada
# *** with whole dataset, just ran and ran and ran, at least on my 2015 quad-core Intel i7 MBP

time0 <- Sys.time()

coord.matches <- Sala_CO2_efflux.coords.dt[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt)]

time1 <- Sys.time()
print(time1 - time0)

# # ******************************************************************************
# # second approach: using data.table and apply to take advantage of parallelization
# # ******************************************************************************
# 
# # *** unfortunately, this approach didn't seem to work so well when I sent it to
# # a c5a.8xlarge 32-vCPU AWS EC2 machine with 64 GB of RAM (and a 35 GB disk swap):
# # ran for about half a day (benchmark time difference of 13.84486 hours),
# # which seemed right, but then returned this error:
# #
# # Warning message:                                                        
# # In parallel::mclapply(X = X, FUN = FUN, ...) :
# #  scheduled cores 1, 3, 5, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32 did not deliver results, all 
# #  values of the jobs will be affected
# #
# # the output object did save, with the correct number of rows, but I don't think
# # it worked correctly ... will try approach #1, above, on a similar EC2 instance
# 
# # convert Sala data to a data.frame, move decimal
# Sala_CO2_efflux.coords.df <- as.data.frame(Sala_CO2_efflux.coords/10^5)
# 
# # create a parallel version of sapply
# 
# mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
#   FUN <- match.fun(FUN)
#   answer <- parallel::mclapply(X = X, FUN = FUN, ...)
#   if (USE.NAMES && is.character(X) && is.null(names(answer))) 
#     names(answer) <- X
#   if (!isFALSE(simplify) && length(answer)) 
#     simplify2array(answer, higher = (simplify == "array"))
#   else answer
# }
# 
# # define function to find nearest match 
# # adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame
# 
# dist2 <- function(df){
#   dt <- data.table((Siegel_fseq.coords.df$x-df$x)^2+(Siegel_fseq.coords.df$y-df$y)^2)
#   return(which.min(dt$V1))}
# 
# # test with a subset first; with benchmarking
# 
# Sala_CO2_efflux.coords.df.sub <- Sala_CO2_efflux.coords.df[1:10000,]
# 
# # find matches; may need to adjust no. of cores
# 
# time0 <- Sys.time()
# 
# coord.matches.test2 <- mcsapply(1:nrow(Sala_CO2_efflux.coords.df.sub), function(x) return(dist2(Sala_CO2_efflux.coords.df.sub[x,])), mc.cores = 4)
# 
# time1 <- Sys.time()
# print(time1 - time0)
# 
# # *** definitely faster than the approach #1 above, at least on my 2015 quad-core Intel i7 MBP
# 
# # save output
# 
# save(coord.matches.test2, file = "output/coord.matches.test2.RData")
# 
# # now with the whole enchilada
# 
# time0 <- Sys.time()
# 
# coord.matches <- mcsapply(1:nrow(Sala_CO2_efflux.coords.df), function(x) return(dist2(Sala_CO2_efflux.coords.df[x,])), mc.cores = 32)
# 
# time1 <- Sys.time()
# print(time1 - time0)

# save output

save(coord.matches, file = "output/coord.matches.RData")

# send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.

system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))


