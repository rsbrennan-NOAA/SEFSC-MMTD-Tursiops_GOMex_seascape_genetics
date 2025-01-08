
library(terra)

setwd("C:/Users/Reid.Brennan/Downloads")
dat <- terra::rast("oisst-avhrr-v02r01.20240901.nc")

https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
as.date(dat$collection.dat)
years <- format(as.Date(dat$collection.dat, format="%m/%d/%Y"),"%Y")

uniq_years <- unique(years)




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

# plot and look at them, to make sure there are no mistakes
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


plot(clim[[1]], ylim=c(20,32), xlim=c(-99, -79))
points(coords, pch=21, col="black", bg="red2", cex=1.3)


library(maps)
library(sf)

# Get states data
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

# Convert raster to dataframe for ggplot
clim_df <- as.data.frame(clim, xy = TRUE)
names(clim_df)[3] <- "value"  # rename the value column

# Create the plot
p <- ggplot() +
  # Plot the raster data
  geom_tile(data = clim_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  # Add the points
  
  geom_point(data = coords, aes(x = lon, y = lat),
             fill="red2",
             shape= 21, color="black", size = 2.5,
             alpha=0.5) +
  # Set the coordinate limits
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  # Customize the theme
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  # Customize the color scale for the raster
  scale_fill_viridis_c() +
  # Labels
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  annotation_scale()

ggsave("figures/map_annualSST_pops.pdf", p, h=4, w=5)

p <- ggplot() +
  # Plot the raster data
  geom_tile(data = clim_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  # Set the coordinate limits
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  # Customize the theme
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  # Customize the color scale for the raster
  scale_fill_viridis_c() +
  # Labels
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  annotation_scale()

ggsave("figures/map_annualSST.pdf", p, h=4, w=5)
ggsave("figures/map_annualSST.png", p, h=4, w=5)


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
write.csv(dfout, file="analysis/environmental_variables/annual_mean_salinity.csv", 
          quote=F, row.names=F)









# ------------------------------
# oxygen
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_oxygen/gom_all_o00_01.nc")
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
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
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

  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_oxygen.csv", 
          quote=F, row.names=F)









# ------------------------------
# annual_mean_nitrate 
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_nitrate/gom_all_n00_01.nc")
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
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
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
  
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_nitrate.csv", 
          quote=F, row.names=F)





# ------------------------------
# annual_mean_phosphate 
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_phosphate/gom_all_p00_01.nc")
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
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
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
  
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_phosphate.csv", 
          quote=F, row.names=F)






