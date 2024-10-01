#map

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)


# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

head(dat)

usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))


# 
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = dat, aes(x = long, y = lat, fill="burntorange"), 
             shape= 21, color="black", size = 2.5,
             alpha=0.5) +
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("latitude")+
  ylab("longitude") +
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  annotation_scale()
  
p

ggsave("figures/map.pdf", p, h=4, w=5)
ggsave("figures/map.png", p, h=4, w=5)


library(marmap)

#--------------#
#
# Calculate least-cost distances to shore
#
#--------------#

# from https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/2.Prepare_spatial_data/2.prepare_spatial_data.R

# Import coordinates of sites

coords.gps = dplyr::select(dat, long , lat)

str(coords.gps)

# Get bathymetry data from NOAA using marmap 
bathydata = getNOAA.bathy(lon1 = min(coords.gps$long) -1,
                          lon2 = max(coords.gps$long)+1,
                          lat1 = min(coords.gps$lat)-1,
                          lat2 = max(coords.gps$lat) +1,
                          resolution = 2)

# Get depth of coordinates

depths = cbind(site = dat$Sample,
               get.depth(bathydata, coords.gps, locator = FALSE))
depths
# one is > 0, so move it slightly. This is a bug bc of the grid
# this also causes problems with distance calcs
coords.gps$long[depths$depth>0] <- -82.58 # from -82.56

depths = cbind(site = dat$Sample,
               get.depth(bathydata, coords.gps, locator = FALSE))
depths
# depths in meters

# Plot bathymetry data and coordinates
#coords <- coords %>%
#  mutate(colors = case_when(
#    Pop.Structure.Location == "Atlantic" ~ "#eac435",
#    Pop.Structure.Location == "Dry Tortuga" ~ "#557fc3",
#    Pop.Structure.Location == "NGOMex" ~ "#03cea4",
#    Pop.Structure.Location == "WGOMex" ~ "#fb4d3d",
#    TRUE ~ "#999999"  # Default color if none of the above match
#  ))

plot(bathydata, n=20, deepest.isobath = min(depths$depth), 
     shallow=0, step=10)
points(dat$long, dat$lat, pch = 21, 
       bg ="orange", col = "black", cex = 2)


# calculate the distance to nearest shoreline:

d <- dist2isobath(bathydata, x=dat$long, y=dat$lat, isobath = 0)
# isobath=0 means find coastline.
# output distance is in meters


# Visualize the great circle distances
blues <- c("lightsteelblue4","lightsteelblue3","lightsteelblue2","lightsteelblue1")
plot(bathydata, image=TRUE, lwd=0.1, land=TRUE, bpal = list(c(0,max(bathydata),"grey"), c(min(bathydata),0,blues)))
plot(bathydata, deep=-200, shallow=-200, step=0, lwd=0.6, add=TRUE)
points(dat$long,dat$lat, pch=21, col="orange4", bg="orange2", cex=.8)
linesGC(d[2:3],d[4:5])

head(d)

# make df of depth and distance to shore and write

dfout <- data.frame(
                    id = dat$Sample,
                    depth = depths$depth,
                    distance_to_shore = d$distance
)

write.csv(dfout, file="analysis/environmental_variables/depth_distance.csv", 
          row.names=F, quote=F)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# get pairwise geographic distances between indivs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# I thought I did this... 
# min depth 1 stops them from crossing land
# these are slow to run, especially lc.dist
trans1 <- trans.mat(bathydata, min.depth = -1, max.depth = NULL)
lc_paths <- lc.dist(trans1, coords.gps, res = "path")
# Compute least-cost distances (km) matrix
lc_dist <- lc.dist(trans1, coords.gps, res = "dist")
save(lc_paths, file = "least_cost_paths.RData")
save(lc_dist, file = "least_cost_paths_dist.RData")


plot.bathy(bathydata, image= TRUE, land = TRUE, n = 0,
           bpal = list(c(0, max(bathydata), "grey"),
                       c(min(bathydata), 0, "royalblue")))
lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)

# Compute least-cost distances (km) matrix

# Convert to matrix, rename columns and rows, and export as csv file
lc_mat = as.matrix(lc_dist)
colnames(lc_mat) = as.vector(location$Sample)
rownames(lc_mat) = as.vector(location$Sample)
lc_mat

# calculate geographic straight line distances directly. to resolve the close together samples that are called 0.

library(geodist)
distmat <- geodist (coords, measure = "geodesic")
distmat_km <- distmat/1000
                
lc_mat[lc_mat == 0] <- distmat_km[lc_mat == 0]

lc_mat[1:10,1:10]

str(lc_mat)                                        

write.csv(lc_mat, file="analysis/environmental_variables/lc_distances_km.csv")






