
library(terra)

setwd("C:/Users/Reid.Brennan/Downloads")
dat <- terra::rast("oisst-avhrr-v02r01.20240901.nc")

https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
as.date(dat$collection.dat)
years <- format(as.Date(dat$collection.dat, format="%m/%d/%Y"),"%Y")

uniq_years <- unique(years)
uniq_years[order(uniq_years)]

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
  sample_ext <- distance(x=tmp_clim, y=sample_vect)
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
  val<-terra::extract(x=tmp_clim, y=close_coords)
  dfout$weekly_mean_temp[missing_index[i]] <- val[1,2]
  
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
# temperature anomaly
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





