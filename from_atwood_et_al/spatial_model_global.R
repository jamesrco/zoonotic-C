#setwd("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data")
#install.packages('ranger')
#install.packages('randomForest')
#install.packages('plyr')
#install.packages('fields')

library(dismo)
library(raster)
library(rgdal)
library(randomForest)
library(plyr)


mg <- read.csv("final_mg.csv", header = TRUE)
sg <- read.csv("final_sg.csv", header = TRUE)
sm <- read.csv("final_sm.csv", header = TRUE)
blue_carbon <-rbind.fill(mg, sg, sm)



coordinates(blue_carbon) =~ Longitude.dd + Latitude.dd
proj4string(blue_carbon)<- CRS("+proj=longlat +datum=WGS84")
require(raster)
shapefile(blue_carbon, "blue_carbon_fnl_1.shp", overwrite = TRUE)

###############################################################################################################
# In ArcGIS I was able to clip the points to within 50km of the coastline using a select by location


###############################################################################################################
#Now I'm bringing back in the shapefile created in ArcGIS

#load shapefile of blue_carbon
bc <-shapefile('blue_carbon_fnl_1.shp')


#load predictors from World Clim
#climate <- getData('worldclim', var='bio', res=2.5)

ext <- extent(-180.0, 180.0, -80.0, 80.0)


mean_temp <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/mean_temp_re9.tif")
mean_temp <- setExtent(mean_temp, ext, keepres=TRUE)
precip <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/precip_re3.tif")
precip <- setExtent(precip, ext, keepres=TRUE)
sst <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/sst_re2.tif")
sst <- setExtent(sst, ext, keepres=TRUE)
chla <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/chla_re3.tif")
chla <- setExtent(chla, ext, keepres=TRUE)
dem <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/dem_re2.tif")
dem <- setExtent(dem, ext, keepres=TRUE)
ssh <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/sea_ht_re2.tif")
ssh <- setExtent(ssh, ext, keepres=TRUE)

#BioClim mean temp and SST combined into temp_re2.tif
temp <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/temp_re10.tif")
temp <- setExtent(temp, ext, keepres=TRUE)

###########################################################################

preds <- stack(temp, precip, chla, dem, ssh)

##############################################################################################

names(preds)[1] <- "temp_re10"
names(preds)[2] <- "precip_re3"
names(preds)[3] <- "chla_re3"
names(preds)[4] <- "dem_re2"
names(preds)[5] <- "sea_ht_re2"

#####################################################################################

#add bc and covariates to one data frame (buffer in meters)
bc_covar <- data.frame(coordinates(bc), bc, extract(temp, bc, method = 'bilinear'), extract(precip, bc, method = 'bilinear'), 
                       extract(chla, bc, method = 'bilinear'), extract(dem, bc, method = 'bilinear'), 
                       extract(ssh, bc, method = 'bilinear'))



#here i renamed the covariate column names to match the raster layer in the preds raster stack 
#(needed to make the prediction work, it pretty much fails otherwise), thanks to all the previous steps for fucking up the names over and over and over
names(bc_covar)[8] <- "stock.1m"
names(bc_covar)[12] <- "temp_re10"
names(bc_covar)[13] <- "precip_re3"
names(bc_covar)[14] <- "chla_re3"
names(bc_covar)[15] <- "dem_re2"
names(bc_covar)[16] <- "sea_ht_re2"


# random Forest model on blue carbon data
print(bc_rf <- randomForest(stock.1m ~ temp_re10 + precip_re3 + chla_re3 + dem_re2 + sea_ht_re2
                      , data= bc_covar, mtry=1, importance=TRUE, na.action=na.omit))

#Plot of error improvements
plot(bc_rf)
#Variable Importance Plot
varImpPlot(bc_rf)

#Predict new conditions with random Forest model
rp_raster <- predict(preds, bc_rf, na.rm = TRUE, filename="bcstock_global_rf.tif", overwrite = TRUE)
#plot(rp_raster) #just a check


#now resample to the correct cell size
#loading raster guide for resampling procedure, cell size 1km X 1km = 1 HA
resample_grid <- raster("C:/Users/Andy/Desktop/Carbon Mapping/Spatial Data/resample_grid.tif")
rp_raster <- resample(rp_raster, resample_grid, method = 'ngb', filename="bcstock_rf_fnl_output_global.tif", overwrite = TRUE)

#now we can plot the result, which has a grid area of 1HA (This was needed to match the carbon stock calculations)
spplot(rp_raster)

# HOORAY! 

#bc_rf_results <- data.frame(predict(bc_rf))
#bc_cov_rf_results <- cbind(bc_covar,bc_rf_results)














###################################################################################

#I want to try now with the sdm package


