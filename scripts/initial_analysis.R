#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate environmental distances between samples
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
# https://thewanglab.github.io/algatr/articles/enviro_data_vignette.html
# check for collinearity

# read in all env. data

dat_depth_dist <- read.csv("analysis/environmental_variables/depth_distance.csv") 
dat_lcdist <- read.csv("analysis/environmental_variables/lc_distances_km.csv")
dat_temp_mean <- read.csv("analysis/environmental_variables/weekly_mean_temp.csv")
dat_temp_anom <- read.csv("analysis/environmental_variables/weekly_anomaly_temp.csv")
dat_temp_annual <- read.csv("analysis/environmental_variables/annual_mean_temp.csv")
dat_salinity_annual <- read.csv("analysis/environmental_variables/annual_mean_salinity.csv")
dat_oxygen_annual <- read.csv("analysis/environmental_variables/annual_mean_oxygen.csv")
dat_nitrate_annual <- read.csv("analysis/environmental_variables/annual_mean_nitrate.csv")
dat_phosphate_annual <- read.csv("analysis/environmental_variables/annual_mean_phosphate.csv")



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
  annual_mean_phosphate = dat_phosphate_annual[, 5],
  lc_dist_km = dat_lcdist[, 5]
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




# check below here--


# nj tree. 
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
p <- ggtree(out_tree) + theme_tree()

  
pout <- p %<+% dm + 
    #geom_tiplab(aes(color=newpop), size=0.9) +
    theme(legend.position="right")+ 
    #geom_text( show.legend  = F ) +
    geom_tippoint(aes(color=newpop), size=4, alpha=0.7) 

ggsave(pout, filename="figures/nj_tree.png", h=5, w=6)


# make correlation plot


# drop correlated env. variables


# calculate environmental distances
# make matrix of geographic distances + environmental distances
# a little unclear exactly what format the env. distances should be. 
# A list with one matrix for each parameter, I think. but double check.
env <- dat_depth_dist$depth 
distmat <- as.matrix(dist(env, diag = TRUE, upper = TRUE))
row.names(distmat) <-dat_depth_dist$id
colnames(distmat)<-dat_depth_dist$id
distmat <- list(distmat)
names(distmat) <- c("depth")

env <- dat_depth_dist$distance_to_shore 
distmat_tmp <- as.matrix(dist(env, diag = TRUE, upper = TRUE))
row.names(distmat_tmp) <-dat_depth_dist$id
colnames(distmat_tmp)<-dat_depth_dist$id
distmat[[2]] <- distmat_tmp
names(distmat) <- c("depth", "distance_From_shore")

# genetic distances matrix needs to be symmetrical
gendist2 <- as.matrix(Matrix::forceSymmetric(gendist,uplo="L"))

results_full <- mmrr_run(gendist2, distmat, nperm = 99, stdz = TRUE, model = "full")

mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "vars", stdz = TRUE)
mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "fitted", stdz = TRUE)
mmrr_plot(gendist2, distmat, mod = results_full$mod, plot_type = "cov", stdz = TRUE)
mmrr_table(results_full, digits = 2, summary_stats = TRUE)

######
# gdm

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)
envs <- dat_depth_dist[,2:3]

gdm_full <- gdm_run(
  gendist = gendist2,
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

gdm_best <- gdm_run(gendist = gendist2, 
                    coords = coords, 
                    env = envs, 
                    model = "best", 
                    scale_gendist = TRUE,
                    nperm = 1000, 
                    sig = 0.05)


# Look at p-values
gdm_best$pvalues
gdm_best$varimp

summary(gdm_full$model)


gdm_table(gdm_full)
