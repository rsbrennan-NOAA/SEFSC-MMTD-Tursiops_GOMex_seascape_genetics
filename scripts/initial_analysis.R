#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate environmental distances between samples
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library(algatr)
# https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
# https://thewanglab.github.io/algatr/articles/enviro_data_vignette.html
# check for collinearity

# read in all env. data

dat_depth_dist <- read.csv("analysis/environmental_variables/depth_distance.csv") 
dat_temp_mean <- read.csv("analysis/environmental_variables/weekly_mean_temp.csv")
dat_temp_anom <- read.csv("analysis/environmental_variables/weekly_anomaly_temp.csv")
dat_temp_annual <- read.csv("analysis/environmental_variables/annual_mean_temp.csv")
dat_salinity_annual <- read.csv("analysis/environmental_variables/annual_mean_salinity.csv")
dat_oxygen_annual <- read.csv("analysis/environmental_variables/annual_mean_oxygen.csv")
dat_nitrate_annual <- read.csv("analysis/environmental_variables/annual_mean_nitrate.csv")
dat_phosphate_annual <- read.csv("analysis/environmental_variables/annual_mean_phosphate.csv")

# pairwise distance between samples. 
dat_lcdist <- read.csv("analysis/environmental_variables/lc_distances_km.csv", row.names=1)



# merge into single df
dat <- (cbind(cbind(dat_depth_dist, dat_temp_mean$weekly_mean_temp), dat_temp_anom$anom_temp_mean))
colnames(dat) <- c("id", "depth", "distance_to_shore", "weekly_mean_temp", "anom_temp")
# are there differences between the genetic pops for env factors?
  # some boxplots here.

dat <- cbind(
  dat_depth_dist[, c("id", "depth", "distance_to_shore")],
  weekly_mean_temp = dat_temp_mean[, 5],
  anom_temp = dat_temp_anom[, 5],
  annual_mean_temp = dat_temp_annual[, 5],
  annual_mean_salinity = dat_salinity_annual[, 5],
  annual_mean_oxygen = dat_oxygen_annual[, 5],
  annual_mean_nitrate = dat_nitrate_annual[, 5],
  annual_mean_phosphate = dat_phosphate_annual[, 5]
)

pops <- read.csv("population_assignments.csv")

head(dat)

all <- merge(dat, pops, by.x="id", by.y="indiv")

boxplot(all$depth ~ all$newpop)
boxplot(all$distance_to_shore  ~ all$newpop)
boxplot(all$weekly_mean_temp   ~ all$newpop)
boxplot(all$annual_mean_oxygen     ~ all$newpop)
boxplot(all$annual_mean_salinity     ~ all$newpop)
boxplot(all$annual_mean_temp    ~ all$newpop)
boxplot(all$annual_mean_phosphate     ~ all$newpop)


#--------------------------------------------------------------------------------
# process environmental data

# make correlation plot to check for collinearity
# Calculate correlation matrix
corr_matrix <- cor(dat[2:ncol(dat)], method = "pearson")

library(corrplot)
corrplot.mixed(corr_matrix,
         upper = "ellipse",  # Ellipses on upper triangle
         lower = "number",   # Numbers on lower triangle
         tl.col = "black",   # Text label color
         tl.cex = 0.7,       # Text label size
         number.cex = 0.7)   # Number size

corrplot(corr_matrix,
         method='ellipse',
              type = 'lower', diag = TRUE,   
               tl.col = "black",   # Text label color
               tl.cex = 0.7,       # Text label size
               number.cex = 0.7)   # Number size

png("figures/correlation_plot.png", width=600, height=600, res=120)
corrplot(corr_matrix,
         method='number',
         type = 'lower', diag = TRUE,   
         tl.col = "black",   # Text label color
         tl.cex = 0.7,       # Text label size
         number.cex = 0.7)   # Number size
dev.off()

# drop correlated env. variables
# phosphate and nitrate 0.82. 
# salinity and temp (annual) are 0.71, borderline. keep for now.

# drop phosphate:
library(dplyr)
dat <- dat %>% select(-annual_mean_phosphate)
 

#--------------------------------------------------------------------------
# calculate environmental distances
# A list with one matrix for each parameter
numeric_cols <- names(dat)[sapply(dat, is.numeric)]
numeric_cols <- numeric_cols[numeric_cols != "id"]

# Initialize empty list
distmat <- list()

# Create distance matrix for each numeric column
for(col in numeric_cols) {
  # Create distance matrix
  dist_tmp <- as.matrix(dist(dat[[col]], diag = TRUE, upper = TRUE))
  
  # Add row and column names
  row.names(dist_tmp) <- dat$id
  colnames(dist_tmp) <- dat$id
  
  # Add to list
  distmat[[col]] <- dist_tmp
}


# add geographic distance:
distmat[['geographic_dist']] <- as.matrix(geo_dist)

# collinearity between geographic and environmental distances:

library(viridis)
# Make a fun heat map with the pairwise distances
geo_dist <- as.data.frame(dat_lcdist)
colnames(geo_dist) <- rownames(geo_dist)


geo_dist %>%
  rownames_to_column("sample") %>%
  gather("sample_comp", "dist", -"sample") %>%
  ggplot(aes(x = (sample), y = (sample_comp), fill = dist)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis() +
  xlab("Sample") +
  ylab("Sample")










#-----------------------------------------------------------------
# nj tree. genetic distances
#-----------------------------------------------------------------
library(ape)
library(Matrix)
library(tidyverse)
library(RColorBrewer)

gendist <- as.matrix(read.csv("genetic_distances_pairwise_matrix.csv"))
gendist2 <- as.matrix(Matrix::forceSymmetric(gendist,uplo="L"))

all$newpop <- as.factor(all$newpop)

Colorsdf <-
  with(all,
       data.frame(population = levels(all$newpop),
                  color = I(brewer.pal(nlevels(newpop), name = 'Dark2'))))
cols <- Colorsdf$color[match(all$newpop, Colorsdf$population)]

nj(gendist2) %>% plot(.,"unrooted", tip.color = cols)
nj(gendist2) %>% plot(., tip.color = cols)

library(ggtree)
out_tree <- nj(gendist2)
out <- as_tibble(nj(gendist2))
out$label <- gsub("X", "", out$label)
dm <- left_join(out, all, by=c("label" ="id"))

ggtree(out_tree) + 
  theme_tree()

#p <- ggtree(out_tree,layout="daylight") + theme_tree()
p <- ggtree(out_tree, layout="unrooted") + theme_tree()
p <- ggtree(out_tree, layout="ape") + theme_tree()


pout <- p %<+% dm + 
  #geom_tiplab(aes(color=newpop), size=0.9) +
  theme(legend.position="right")+ 
  #geom_text( show.legend  = F ) +
  geom_tippoint(aes(color=newpop), size=4, alpha=0.7) 

ggsave(pout, filename="figures/nj_tree_unrooted.png", h=5, w=6)



#-----------------------------------------------------------------
# prep genetic distances for gea
#-----------------------------------------------------------------


# genetic distances matrix needs to be symmetrical
gendist2 <- as.matrix(Matrix::forceSymmetric(gendist,uplo="L"))

#-----------------------------------------------------------------
# mmr: multiple matrix regression: 
#-----------------------------------------------------------------
#https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
# disentangling the contribution of geographic and environment

# significance found via randomization permutaions
# get ndividual regression coefficients and p-values for each dependent variable and a 
  # “coefficient ratio,” which is the ratio between regression coefficients, 
   # which thus provides an idea of the relative contributions of IBD and IBE in explaining variation in the genetic distances in your data.

# assumptions:
    # (1) the coordinates and genetic distance files MUST have the same ordering of individuals;
    # (2) this function assumes that each individual has its own sampling coordinates (even if population-based sampling was performed).

set.seed(01)
results_full <- mmrr_run(gendist2, distmat, nperm = 99, stdz = TRUE, model = "full")
# The results from running the “full” MMRR model contains four elements:
  # 1. coeff_df: a dataframe with statistics relating to each variable’s distance related to genetic distance, including coefficient values for each environmental variable and geographic distance
  # 2 mod: a dataframe containing statistics for the results of the model, including an R^2 value, and F statistics
  # 3/4: X and Y: the input data

#results_full

table_full <- mmrr_table(results_full, digits = 2, summary_stats = TRUE)
gt::gtsave(table_full, filename = "figures/mmrr_full_results_table.png")   # Save as PNG


# Single variable plot
variables_full <- mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "vars", stdz = TRUE)
ggsave(filename="figures/mmrr_variables_full.png", h=7, w=7)

# Fitted variable plot
fitted_full <- mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "fitted", stdz = TRUE)
ggsave(filename="figures/mmrr_fitted_full.png", h=7, w=7)

# covariance plot
covariance_full <- mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "cov", stdz = TRUE)
ggsave(filename="figures/mmrr_covariance_full.png", h=7, w=7)


#### best model only #######
# Run MMRR with all variables and select best model
set.seed(01)
results_best <- mmrr_run(gendist2, distmat, nperm = 99, stdz = TRUE, model = "best")

table_best <- mmrr_table(results_best, digits = 2, summary_stats = TRUE)
gt::gtsave(table_best, filename = "figures/mmrr_best_results_table.png")   # Save as PNG




mmrr_plot(gendist2, distmat, mod = results_best$mod, plot_type = "all", stdz = TRUE)








#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# gdm
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# Generalized dissimilarity modeling (GDM) is a matrix regression method in which 
# explanatory variables (in our case, genetic data, in the form of a distance matrix) 
# is regressed against a response matrix (e.g., environmental variables for locations 
# from which samples were obtained and geographic distances between those locations).
# GDM calculates the compositional dissimilarity between pairs of sites, 
# and importantly allows for nonlinear relationships to be modeled.


# assumptions:
  # (1) the coords and gendist files must have the same ordering of individuals; 
  # (2) each individual has its own sampling coordinates
location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)
envs <- dat[,2:ncol(dat)]

gdm_full <- gdm_run(
  gendist = gendist,
  coords = coords,
  env = envs,
  model = "full",
  scale_gendist = TRUE
)

summary(gdm_full$model)

gdm_plot_diss(gdm_full$model)
# Plot the I-splines with free x and y-axes
gdm_plot_isplines(gdm_full$model, scales = "free")

# Plot the I-splines with a free x-axis and a fixed y-axis
# This allows for visualization of relative importance (i.e., the height of the I-splines)
gdm_plot_isplines(gdm_full$model, scales = "free_x")
gdm_table(gdm_full)

# calculate variable importance and significance:
# To run gdm.varImp() you need a gdmData object, which you can create using gdm_format()
gdmData <- gdm_format(gendist, coords, envs, scale_gendist = TRUE)

# Then you can run gdm.varImp(), 
# specifying whether you want to use geographic distance as a variable as well 
varimp <- gdm::gdm.varImp(gdmData,
                          predSelect = FALSE ,
                          geo = TRUE, nPerm = 50,parallel=TRUE, cores = 9)

# visualize the results 
gdm_varimp_table(varimp)

# depth by far the most important. then annual mean temp, 
    # then annual mean nitrate, then dist to shore. then annual mean oxygen

# I think some of these are confounded, so depth and distance to shore and temp maybe?

# make a map, but I need a raster to do this. doesn't work currenty. 
map <- gdm_map(gdm_full$model, envs, coords)


gdm_best <- gdm_run(gendist = gendist, 
                    coords = coords, 
                    env = envs, 
                    model = "best", 
                    scale_gendist = TRUE,
                    nperm = 10, # should be 1000, just 10 to run quickly to see if working
                    sig = 0.05)

# Look at p-values
gdm_best$pvalues
gdm_best$varimp




