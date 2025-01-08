#!/bin/bash

# get environmental data for Gulf of Mexico

# weekly tempature mean
## from 1981-2024.09.15
## 0.25 degree grids
mkdir -p analysis/environmental_variables/temperature_oisst/

wget -P analysis/environmental_variables/temperature_oisst/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.week.mean.nc

# get daily temperature anomaly to calculate weekly temp anomaly
## for each year
## 0.25 degree grids
mkdir -p analysis/environmental_variables/temperature_anom/

wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.1994.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.1996.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.1997.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.1999.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2000.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2001.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2002.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2003.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2004.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2006.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2007.nc
wget -P analysis/environmental_variables/temperature_anom/ https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.anom.2008.nc

######
# annual means
######
#  All from https://www.ncei.noaa.gov/

# temperature:
## 1/10 degree grid
## mean from 1995-2024

mkdir -p analysis/environmental_variables/annual_mean_sst_1995_2004/
wget -P analysis/environmental_variables/annual_mean_sst_1995_2004/ https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/temperature/netcdf/95A4/0.10/gom_95A4_t00_10.nc

# salinity 
## 1/10 degree grids
## annual mean from 1995_2004 

mkdir -p analysis/environmental_variables/annual_mean_salinity_1995_2004
wget -P analysis/environmental_variables/annual_mean_salinity_1995_2004 https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/salinity/netcdf/95A4/0.10/gom_95A4_s00_10.nc#cd ..

# nitrate 
## 1 degree grids  
## annual means 1955–2017
mkdir -p analysis/environmental_variables/annual_mean_nitrate
wget -P analysis/environmental_variables/annual_mean_nitrate https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/nitrate/netcdf/all/1.00/gom_all_n00_01.nc

# oxygen 
## 1 degree grids
## annual means 1955–2017
mkdir -p analysis/environmental_variables/annual_mean_oxygen 
wget -P  analysis/environmental_variables/annual_mean_oxygen  https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/oxygen/netcdf/all/1.00/gom_all_o00_01.nc

# phosphate - 
# 1° grid
# assnual mean 1955–2017
mkdir -p analysis/environmental_variables/annual_mean_phosphate
wget -P analysis/environmental_variables/annual_mean_phosphate https://www.ncei.noaa.gov/thredds-ocean/fileServer/regional-climatologies/gulf-of-mexico/DATA/phosphate/netcdf/all/1.00/gom_all_p00_01.nc

