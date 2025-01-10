# fig1: Combined map + PCA and DAPC. map on top, dapc and pca on bottom. 

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
library(RColorBrewer)

# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

head(dat)

usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

#---------------
# add pop colors:

# read in the pop labels.
my_colors <- data.frame(color = c("firebrick3", "#1B9E77", "#8B4513", "#E6AB02", 
                                  "black", "#7570B3"),
                        population = c("Eastern Coastal", "Eastern Oceanic", 
                                       "Northeastern Oceanic", "Western Coastal",
                                       "Western Oceanic", "Shelf"))

# read in populations
population_ids <- read.csv("analysis/population_assignments.csv", header=T)
pop_ids_new <- full_join(dat, population_ids, by=c("Sample" = "indiv"))

table(population_ids$population_name)

# 
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = pop_ids_new, aes(x = long, y = lat, fill=population_name),
             shape= 21, color="black", size = 2.5,
             alpha=1) +
  coord_sf() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),

    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") +   
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  annotation_scale() +
  scale_fill_manual(values = setNames(my_colors$color, my_colors$population),
                    guide = guide_legend(override.aes = list(alpha = 1, size = 2.5)))

p


ggsave("figures/map_pop_assignments.pdf", p, h=4, w=5)
ggsave("figures/map_pop_assignments.png", p, h=4, w=5)

write.csv(my_colors, file="analysis/popColors.csv", row.names=F)

##----------------------------------------------------
# make PCA
##----------------------------------------------------

p1 <- ggplot(population_ids, aes(x=PC1, y=PC2, fill=population_name))+
  geom_point(size=3, pch=21) +
  ggtitle("PCA")  +
  xlab("PC1: 25% variation")+
  ylab("PC2: 3% variation") +
  theme_classic(base_size = 14) +
    scale_fill_manual(values = setNames(my_colors$color, my_colors$population),
                    guide = guide_legend(override.aes = list(alpha = 1, size = 2.5)))

p1


p2 <- ggplot(pltdat, aes(x=LD1, y=LD2*-1, fill=population_name))+
  geom_point(size=3, pch=21) +
  ggtitle("DAPC")  +
  ylab("LD2") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = setNames(my_colors$color, my_colors$population),
                    guide = guide_legend(override.aes = list(alpha = 1, size = 2.5)))
p2

ggsave(ggpubr::ggarrange(p1, p2, common.legend = T), file="figures/pca_vs_dapc.png",
       h=4, w=7)

# combine map and pca:










