# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# libraries

library(ncdf4)

# load datasets

# Attempt to load Siegel et al. "sequestration fractions" from NetCDF

OCIM2_fseq_48L <- nc_open("data/siegel_et_al_2021_v2/fseq_OCIM2_48L.nc", verbose=TRUE)

# currently returns an error ... can use the MATLAB script provided by Siegel et al. to get what we need a bit quicker 

# Sala et al. pCO2 data, as GeoTIFF

