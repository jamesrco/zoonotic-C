# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# libraries

library(ncdf4)

# load datasets

# Attempt to load Siegel et al. "sequestration fractions" from NetCDF

OCIM2_fseq_48L <- nc_open("data/siegel_et_al_2021_v2/fseq_OCIM2_48L.nc", verbose=TRUE) # currently returns an error

# instead we'll use a modified version of the MATLAB script provided by Siegel et al. to get what we need: see "gen_fracs_to_constrain_trawlCO2.m" which should be in this same directory

# take the necessary detour into MATLAB at this point, if the output hasn't been generated; then return to R ... ok, done.

# the output from gen_fracs_to_constrain_trawlCO2.m -- 25, 50, and 100 year sequestration fractions for ocean bottom depths, plus the necessary metadata -- should now be in several .csv files found in data/derived/benthic_seqfractions

# load these output files


# Sala et al. pCO2 data, as GeoTIFF

