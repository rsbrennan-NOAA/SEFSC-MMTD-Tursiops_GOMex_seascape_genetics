
library(terra)

setwd("C:/Users/Reid.Brennan/Downloads")
dat <- terra::rast("oisst-avhrr-v02r01.20240901.nc")

https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
as.date(dat$collection.dat)
years <- format(as.Date(dat$collection.dat, format="%m/%d/%Y"),"%Y")

uniq_years <- unique(years)
uniq_years[order(uniq_years)]

#------------------------------------------------------------------------------
# calc weekly sst around collection date:

# "1994" "1996" "1997" "1999" "2000" "2001" "2002" "2003" "2004" "2006" "2007" "2008"

wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.1994.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.1996.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.1997.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.1999.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2000.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2001.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2002.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2003.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2004.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2006.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2007.nc
wget https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.2008.nc


# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)

clim <- rast("analysis/environmental_variables/temperature_oisst/sst.week.mean.nc")
plot(clim)

time(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)
# need to adjust these to match those used in the raster file
coords$lon <- 360 +(coords$lon)


library(maps)
head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, get the dates, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

# find closest week to each collection date
# if there's tie, should return the earlier week, which seems correct
date_index <- outer(location$date_correct, time(clim), `-`) |> abs() |> apply(1, which.min)

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "weekly_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct
  
for(i in 1:length(date_index)){
  tmp_clim <- clim[[date_index[i]]]
  val<-terra::extract(x=tmp_clim, y=coords[i,])
  dfout$weekly_mean_temp[i] <- val[1,2]
}


# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$weekly_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(277,278), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=1)

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA


for(i in 1:length(missing_index)){

  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="+proj=longlat +datum=WGS84")
  sample_ext <- distance(x=clim[[date_index[1]]], y=sample_vect)
  df_ext <- values(sample_ext)

  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  # take this index, and pull it from the original 
  # use the same closest date as above
  tmp_clim <- clim[[date_index[missing_index[i]]]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords,
                      method="bilinear")
  val<-terra::extract(x=clim_wgs84, y=close_coords)
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
climplt <- clim[[1]]
plot(climplt, xlim=c(277,278), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/weekly_mean_temp.csv", 
              quote=F, row.names=F)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# weekly temperature anomaly
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# need to handle this differently, bc daily data. 
# find the date, select that day and 7 days prior. take mean anomaly

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")
location$year <- format(as.Date(location$date_correct, format="%d/%m/%Y"),"%Y")


coords<-data.frame(lon=location$long, lat=location$lat)
# need to adjust these to match those used in the raster file
coords$lon <- 360 +(coords$lon)


dfout <- as.data.frame(matrix(ncol=7, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "anom_temp_mean", "anom_temp_median", "n_days")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

# note that 6 samples still fall just outside of the grid. move them to the same grid as above. 
#
for(i in 1:length(missing_index)){
  coords$lon[missing_index[i]] <-missingdat$lon_new[i]
  coords$lat[missing_index[i]] <- missingdat$lat_new[i]
}

# loop over each sample

for(i in 1:nrow(location)){
  # get the year index
  year_tmp <- location$year[i]
  # read in the data
  anom_nm <- paste0("analysis/environmental_variables/temperature_anom/sst.day.anom.", year_tmp, ".nc")
  anom_in <- rast(anom_nm)
  # get date
  date_tmp <- location$date_correct[i]
  date_index <- which(time(anom_in) == date_tmp)
  date_week <- c(date_index, date_index-seq(1, 6, 1))
  # get the anom vals for each day:
  val_tmp <- rep(NA, 7)  
  for(day in 1:length(date_week)){
    tmp_clim <- anom_in[[date_week[day]]]
    anom_tmp <- terra::extract(x=tmp_clim, y=coords[i,])
    val_tmp[day] <- anom_tmp[1,2]
  }
  # get median for the week and save to output
  dfout$anom_temp_mean[i] <- median(val_tmp, na.rm=T)
  dfout$anom_temp_median[i] <-mean(val_tmp, na.rm=T)
  dfout$n_days[i] <- sum(!is.na(val_tmp))
    
}

# check for missing data:
dfout[which(dfout$n_days<7),]

plot(dfout$anom_temp_mean, dfout$anom_temp_median)

# write output:
write.csv(dfout[,1:5], file="analysis/environmental_variables/weekly_anomaly_temp.csv", 
          quote=F, row.names=F)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# gather all data so far:
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
dat_depth <- read.csv("analysis/environmental_variables/depth_distance.csv")
dat_anom <- read.csv("analysis/environmental_variables/weekly_anomaly_temp.csv")
dat_temp <- read.csv("analysis/environmental_variables/weekly_mean_temp.csv")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# long-term data, annual avg
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# sst- annual mean from 1995_2004 1/10 degree grids
#mkdir annual_mean_sst_1995_2004
#cd annual_mean_sst_1995_2004
#wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/temperature/netcdf/95A4/0.10/gom_95A4_t00_10.nc
#cd ..
# salinity - annual mean from 1995_2004 1/10 degree grids
#mkdir annual_mean_salinity_1995_2004
#cd annual_mean_salinity_1995_2004
#wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/salinity/netcdf/95A4/0.10/gom_95A4_s00_10.nc#cd ..
# nitrate 1 degree grids  1955–2017
#mkdir annual_mean_nitrate
#cd annual_mean_nitrate
#wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/nitrate/netcdf/all/1.00/gom_all_n00_01.nc
#cd ..
# oxygen sat 1 degree grids  1955–2017
#mkdir annual_mean_oxygen
#cd annual_mean_oxygen
# wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/oxygen/netcdf/all/1.00/gom_all_o00_01.nc
#cd ..
# silicate - Statistical mean on 1° grid  1955–2017
#wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/silicate/netcdf/all/1.00/gom_all_i00_01.nc
# phosphate - Statistical mean on 1° grid for  1955–2017
#mkdir annual_mean_phosphate
#cd annual_mean_phosphate
#wget https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/phosphate/netcdf/all/1.00/gom_all_p00_01.nc
#cd ..

# ------------------------------
# temperature
# ------------------------------



# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)
library(maps)

clim <- rast("analysis/environmental_variables/annual_mean_sst_1995_2004/gom_95A4_t00_10.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  # take this index, and pull it from the original 
  # use the same closest date as above
  #tmp_clim <- clim[[1]][missing_index[i]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  if(is.na(val[1,2]) == TRUE){
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[3],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_temp.csv", 
          quote=F, row.names=F)






# ------------------------------
# salinity
# ------------------------------

# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)
library(maps)

clim <- rast("analysis/environmental_variables/annual_mean_salinity_1995_2004/gom_95A4_s00_10.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  # take this index, and pull it from the original 
  # use the same closest date as above
  #tmp_clim <- clim[[1]][missing_index[i]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  if(is.na(val[1,2]) == TRUE){
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[3],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_temp.csv", 
          quote=F, row.names=F)







