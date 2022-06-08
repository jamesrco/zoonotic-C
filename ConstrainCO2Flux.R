# ConstrainCO2Flux.R
# Created June 7, 2022
# Purpose: Attempt to constrain the crude estimate of benthic CO2 flux in Sala  et al. 2021
# using the sequestration fractions in Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# libraries

library(ncdf4)

# load datasets

# Siegel et al. "sequestration fractions," from NetCDF

OCIM2_fseq_48L <- nc_open("data/siegel_et_al_2021/fseq_OCIM2_48L.nc", verbose=TRUE)
ncvar_get("data/siegel_et_al_2021/fseq_OCIM2_48L.nc")

# Sala et al. pCO2 data, as GeoTIFF

