#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate environmental distances between samples
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# useful papers:
# uses both GDM and mmrr: https://doi.org/10.1111/mec.15301
# mmrr with backwards elimination: https://www.nature.com/articles/s41437-021-00405-0
# can test the hypothesis that IBD vs IBE

library(algatr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(dplyr)


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

dat_depth_dist$depth_log10 <- log10(dat_depth_dist$depth * -1)*-1
dat_depth_dist$distance_to_shore_log10 <- log10(dat_depth_dist$distance_to_shore)

# merge into single df

dat <- cbind(
  dat_depth_dist[, c("id", "depth", "depth_log10",
                     "distance_to_shore_log10","distance_to_shore")],
  weekly_mean_temp      = dat_temp_mean[, 5],
  weekly_anom_temp      = dat_temp_anom[, 5],
  annual_mean_temp      = dat_temp_annual[, 5],
  annual_mean_salinity  = dat_salinity_annual[, 5],
  annual_mean_oxygen    = dat_oxygen_annual[, 5],
  annual_mean_nitrate   = dat_nitrate_annual[, 5],
  annual_mean_phosphate = dat_phosphate_annual[, 5]
)

pops <- read.csv("analysis/population_assignments.csv")
colnames(pops) <- c("Sample", "pop", "pred.post", "geo_pop", "comp", "Pop")
# pops <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

head(dat)

all <- merge(dat, pops[,c('Sample','Pop')], by.x="id", by.y="Sample")

popColors <- data.frame(population= c("East_Oceanic", "East_shelf", "NE_Oceanic", "NW_InnerShelf", "NW_Oceanic", "NW_Outer_shelf"),
                          color= c("#1B9E77","firebrick3","#7570B3","#E6AB02","#8B4513","black"))

table(all$Pop)

long_data <- all %>%
  select(Pop, 2:10) %>%  # Select columns 5 through 10
  pivot_longer(cols = -Pop,
               names_to = "variable",
               values_to = "value")


new_order <- c("East_shelf", "NW_InnerShelf", "NW_Outer_shelf", "NW_Oceanic", "NE_Oceanic", "East_Oceanic")
ordered_colors <- popColors$color[match(new_order, popColors$population)]  

ggplot(long_data, aes(x = factor(Pop, levels = new_order), 
                y = value, 
                fill = factor(Pop, levels = new_order))) +
  facet_wrap(~variable, 
             scales = "free_y") +
  geom_boxplot() +
  scale_fill_manual(values = ordered_colors) +
  scale_color_manual(values = ordered_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "")

ggsave("figures/all_environments_pops.png", h=10, w=10)

# break out for presentation
p_depth <- ggplot(all, aes(x = factor(Pop, levels = new_order), 
                           y = depth, 
                           fill = factor(Pop, levels = new_order), 
                           color = factor(Pop, levels = new_order))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = ordered_colors) +
  scale_color_manual(values = ordered_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "") +
  coord_cartesian(ylim=c(-900,0))+
  ggtitle("depth")

ggsave(file="figures/popEnv_depth.png", p_depth, h=3, w=5)


p_temp <- ggplot(all, aes(x = factor(Pop, levels = new_order), 
                           y = annual_mean_temp, 
                           fill = factor(Pop, levels = new_order), 
                           color = factor(Pop, levels = new_order))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = ordered_colors) +
  scale_color_manual(values = ordered_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "", y="annual mean temperature") +
  ggtitle("annual mean temperature")

ggsave(file="figures/popEnv_annualTemp.png", p_temp, h=3, w=5)



temp_depth <- ggplot(all, aes(x = log10(depth*-1), 
                y = annual_mean_temp, 
                fill = factor(Pop, levels = new_order))) +
  geom_point(shape=21, size=2.5) +
  scale_fill_manual(values = ordered_colors) +
  scale_color_manual(values = ordered_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "log10(depth)", y="annual mean temperature") +
  ggtitle("temperature vs. depth")
  

ggsave(file="figures/temp_depth.png", temp_depth, h=3.5, w=3.5)


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# process environmental data
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# make correlation plot to check for collinearity
# Calculate correlation matrix
corr_matrix <- cor(dat[2:ncol(dat)], method = "pearson")

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
# salinity and temp (annual) are 0.71,
high_corr <- which(abs(corr_matrix) > 0.7 & corr_matrix != 1 & upper.tri(corr_matrix), arr.ind = TRUE)
data.frame(
  row = rownames(corr_matrix)[high_corr[,1]],
  col = colnames(corr_matrix)[high_corr[,2]], 
  correlation = corr_matrix[high_corr]
)

#row                     col correlation
#1                   depth             depth_log10   0.7540154
#2             depth_log10 distance_to_shore_log10  -0.7645013
#3             depth_log10       distance_to_shore  -0.7217551
#4 distance_to_shore_log10       distance_to_shore   0.8774960
#5             depth_log10        annual_mean_temp  -0.7305640
#6        annual_mean_temp    annual_mean_salinity   0.7063588
#7     annual_mean_nitrate   annual_mean_phosphate   0.8160668

# to drop:
dat <- dat %>%  select(-depth, -distance_to_shore, -distance_to_shore_log10, 
-annual_mean_temp, -annual_mean_phosphate)

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
  #tmp_scaled <- scale(dat[[col]], center = TRUE, scale = TRUE)
  # Create distance matrix
  #dist_tmp <- as.matrix(dist(tmp_scaled, diag = TRUE, upper = TRUE))
  dist_tmp <- as.matrix(dist(dat[[col]], diag = TRUE, upper = TRUE))
  #dist_tmp <- as.matrix(tmp_scaled, diag = TRUE, upper = TRUE, method="euclidean"))
  # Add row and column names
  row.names(dist_tmp) <- dat$id
  colnames(dist_tmp) <- dat$id
  
  # Add to list
  distmat[[col]] <- dist_tmp
}


# add geographic distance:
#distmat[['geographic_dist']] <- as.matrix(geo_dist)
#distmat[['geographic_dist']] <- as.matrix(log10(dat_lcdist))
geo_dist_vector <- (dat_lcdist)

geo_dist_log <- as.matrix(log10(dat_lcdist))
geo_dist_log[is.infinite(geo_dist_log)] <- 0


distmat[['geographic_dist']] <- (geo_dist_log)
#dat$geo_distance <- geo_dist_log


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



# pca of environmental factors

library(vegan)
env_pca <- rda(dat[,2:ncol(dat)], scale=TRUE)
summary(env_pca)




# Get site scores
site_scores <- as.data.frame(scores(env_pca, display = "sites"))
# Get variable loadings (for arrows)
var_scores <- as.data.frame(scores(env_pca, display = "species"))

site_scores$Pop <- all$Pop  # Assuming ids line up, which they do.

# Create data frame for arrows
var_arrows <- data.frame(
  x1 = 0,
  y1 = 0,
  x2 = var_scores$PC1,
  y2 = var_scores$PC2,
  labels = rownames(var_scores)
)
var_arrows$labels <- c("depth", 
                       "weekly temp\nmean", 
                       "weekly temp\nanomaly",
                       "annual salinity\nmean", 
                       "annual oxygen\nmean", 
                       "annual nitrate\nmean")

library(ggrepel)

# Create the plot
pca_plot <- ggplot() +
  # Add points for sites
  geom_point(data = site_scores, 
             aes(x = PC1, y = PC2, fill = factor(Pop, levels = new_order)),
             shape = 21, size = 2.5) +
  # Add arrows
  geom_segment(data = var_arrows,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red") +
  # Add labels for arrows
  geom_text_repel(
    data = var_arrows,
    aes(x = x2, y = y2, label = labels),
    color = "black",
    size = 3.5,
    direction = "both",
    box.padding = 1,
    point.padding = 0.5,
    min.segment.length = 0.2,
    max.overlaps = Inf,
    segment.color = "grey0",
    segment.linetype=3,
    force = 7,
    force_pull = 0.1,
    nudge_x = ifelse(var_arrows$x2 > -0.5, 0.2, -0.1)+ 
      ifelse(var_arrows$y2 == min(var_arrows$y2), 0.5, 0),
    nudge_y = ifelse(var_arrows$y2 > 1, 0.3, -0.9) + 
      ifelse(var_arrows$y2 == min(var_arrows$y2), 1, 0)
  ) +
  scale_fill_manual(values = ordered_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.title = element_blank()
  ) +
  labs(x = "PC1", y = "PC2") +
  ggtitle("PCA of Environmental Variables")
pca_plot

ggsave(pca_plot, file="figures/pca_environment.pdf", h=5, w=6.5)
ggsave(pca_plot, file="figures/pca_environment.png", h=5, w=6.5)

# what are the individuals on the bottom left of the plot?

dat[which(site_scores$PC1 < -1),]

dat[255:262,]
dat_depth_dist[255:262,]

# those are shallow ones.
hist(dat_depth_dist$depth, breaks=500)

dat_depth_dist[dat_depth_dist$depth< -2000,]

# these deep ones are real and legit

#-----------------------------------------------------------------
# prep genetic distances for gea
#-----------------------------------------------------------------

# genetic distances matrix needs to be symmetrical
gendist <- as.matrix(read.csv("genetic_distances_pairwise_matrix.csv"))
gendist2 <- as.matrix(Matrix::forceSymmetric(gendist,uplo="L"))

#-----------------------------------------------------------------
# mmr: multiple matrix regression: 
#-----------------------------------------------------------------
#https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
# disentangling the contribution of geographic and environment

# significance found via randomization permutaions
# get individual regression coefficients and p-values for each dependent variable and a 
  # “coefficient ratio,” which is the ratio between regression coefficients, 
   # which thus provides an idea of the relative contributions of IBD and IBE in explaining variation in the genetic distances in your data.

# assumptions:
    # (1) the coordinates and genetic distance files MUST have the same ordering of individuals;
    # (2) this function assumes that each individual has its own sampling coordinates (even if population-based sampling was performed).

set.seed(01)
results_full <- mmrr_run(gendist2, distmat, nperm = 500, stdz = TRUE, model = "full")
# either need to stdz here or the environmental data directly. I currently have done it above.
# The results from running the “full” MMRR model contains four elements:
  # 1. coeff_df: a dataframe with statistics relating to each variable’s distance related to genetic distance, including coefficient values for each environmental variable and geographic distance
  # 2 mod: a dataframe containing statistics for the results of the model, including an R^2 value, and F statistics
  # 3/4: X and Y: the input data


Y <- gendist2
X<-  distmat
nperm = 4
# figure out the actual MMRR manually
# Compute regression coefficients and test statistics
nrowsY <- nrow(Y)
y <- unfold(Y, FALSE)
if (is.null(names(X))) names(X) <- paste("X", 1:length(X), sep = "")
Xmats <- sapply(X, unfold, scale = TRUE)
fit <- stats::lm(y ~ Xmats)
coeffs <- fit$coefficients
summ <- summary(fit)
r.squared <- summ$r.squared
tstat <- summ$coefficients[, "t value"]
Fstat <- summ$fstatistic[1]
tprob <- rep(1, length(tstat))
Fprob <- 1

# Perform permutations
for (i in 1:nperm) {
  rand <- sample(1:nrowsY)
  Yperm <- Y[rand, rand]
  yperm <- unfold(Yperm, scale)
  fit <- stats::lm(yperm ~ Xmats)
  summ <- summary(fit)
  Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
  tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
}


# Get confidence interval
conf_df <- stats::confint(fit, names(fit$coefficients), level = 0.90)
rownames(conf_df) <- c("Intercept", names(X))


# Perform permutations
for (i in 1:nperm) {
  rand <- sample(1:nrowsY)
  Yperm <- Y[rand, rand]
  yperm <- unfold(Yperm, FALSE)
  fit <- stats::lm(yperm ~ Xmats)
  summ <- summary(fit)
  Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
  tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
}

# Return values
tp <- tprob / (nperm + 1)
Fp <- Fprob / (nperm + 1)
names(r.squared) <- "r.squared"
names(coeffs) <- c("Intercept", names(X))
names(tstat) <- paste(c("Intercept", names(X)), "(t)", sep = "")
names(tp) <- paste(c("Intercept", names(X)), "(p)", sep = "")
names(Fstat) <- "F-statistic"
names(Fp) <- "F p-value"















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
covariance_full <- mmrr_plot(y=gendist2, x= distmat, mod = results_full$mod, plot_type = "cov", stdz = TRUE)
ggsave(filename="figures/mmrr_covariance_full.png", h=7, w=7)



#### best model only #######
# Run MMRR with all variables and select best model
set.seed(01)
results_best <- mmrr_run(gendist2, distmat, nperm = 200, stdz = TRUE, model = "best")
results_best$coeff_df
table_best <- mmrr_table(results_best, digits = 2, summary_stats = TRUE)
gt::gtsave(table_best, filename = "figures/mmrr_best_results_table.png")   # Save as PNG


variables_best <- mmrr_plot(gendist2, distmat, mod = results_best$mod, plot_type = "vars", stdz = TRUE)
variables_best
fitted_best <- mmrr_plot(gendist2, distmat, mod = results_best$mod, plot_type = "fitted", stdz = TRUE)

# to do:
## run model with and without geographic distance. can get idea of the contribution of each
## 








#---------------------------------------------------------------------------------
### can I predict anything?
#---------------------------------------------------------------------------------

unfold <- function(X, scale = TRUE) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:i - 1])
  if (scale == TRUE) x <- scale(x, center = TRUE, scale = TRUE)
  return(x)
}

# Calculate original scaling parameters
dist_reference <-unfold(distmat[['depth']], scale=T)

# Get scaling parameters for prediction
scale_mean <- attr(dist_reference, "scaled:center")
scale_sd <- attr(dist_reference, "scaled:scale")


predict_from_depth <- function(ref_depths,      
                               unknown_depth,     
                               beta = 0.63,  
                               scale_mean,        
                               scale_sd) {        
###
  # Calculate raw distances between unknown and reference populations
  depth_diffs <- abs(ref_depths - unknown_depth)
  
  # Scale using parameters from original
  depth_diffs_scaled <- (depth_diffs - scale_mean) / scale_sd
  
  # Predict genetic distances
  pred_gen_dist <- depth_diffs_scaled * beta
  hist(pred_gen_dist)
  
  results_df <- data.frame(
    Pop = all$Pop,
    Distance = pred_gen_dist
  )
  pop_summary <- results_df %>%
    group_by(Pop) %>%
    summarise(
      Mean_Distance = mean(Distance)
    )                    
  
  return(assigned_pop <- pop_summary$Pop[which.min(abs(pop_summary$Mean_Distance))])
  
}


# Leave one out

for(i in 1:nrow(gendist2)){
  gendist2_subset <- gendist2[-i, -i]
  distmat_subset <- lapply(distmat, function(mat) {
    mat[-i, -i]
  })
  results_full <- mmrr_run(gendist2_subset, distmat_subset, nperm = 1, stdz = TRUE, model = "full")
  
  # pull out coefficient, the unique standardization, then assign pop.
  # what parts of this did I automate in the function above?
}
results_full <- mmrr_run(gendist2, distmat, nperm = 1, stdz = TRUE, model = "full")


# ---------------------------------------------------------------------------
# pop assignment based on depth only

results_df <- data.frame()

for(i in 1:nrow(all)){
  tmp_subset <- all[-i,]
  tmp_unknown <- all[i,]
  
  tmp_subset$difference <- abs(abs(tmp_subset$depth) - abs(tmp_unknown$depth))

  tmp_summary <- tmp_subset %>%
    group_by(Pop) %>%
    summarise(
      Mean_Distance = mean(difference)
    )
  tmp_out <- tmp_summary[which.min(tmp_summary$Mean_Distance),]

  tmp_out$True_pop <- tmp_unknown$Pop
  tmp_out$loo_idx <- i
  colnames(tmp_out) <- c("Predicted_pop", "Mean_Distance", "True_pop", "loo_idx")
  tmp_out$depth <- tmp_unknown$depth
  
  
  results_df <- rbind(results_df, tmp_out)
  
}


sum(results_df$Predicted_pop == results_df$True_pop)/nrow(results_df)
# 65%

results_df$True_region <- case_when(
  results_df$True_pop %in% c("NWInner", "EInner", "EOuter") ~ "Inner",
  results_df$True_pop %in% c("NEOFF", "DeepOff", "EastOff") ~ "Offshore",
  results_df$True_pop %in% c("ShelfOff") ~ "ShelfOff"
)

results_df$Predicted_region <- case_when(
  results_df$Predicted_pop %in% c("NWInner", "EInner", "EOuter") ~ "Inner",
  results_df$Predicted_pop %in% c("NEOFF", "DeepOff", "EastOff") ~ "Offshore",
  results_df$Predicted_pop %in% c("ShelfOff") ~ "ShelfOff"
)

sum(results_df$True_region == results_df$Predicted_region)/nrow(results_df)
# 0.9016698

pop_accuracy <- results_df %>%
  group_by(True_pop) %>%
  summarize(
    accuracy = mean(True_pop == Predicted_pop) * 100
  )

overall_pop_accuracy <- mean(results_df$Predicted_pop == results_df$True_pop) * 100

table(subset(results_df, True_pop=="DeepOff")$Predicted_pop)

region_accuracy <- results_df %>%
  group_by(True_region) %>%
  summarize(
    accuracy = mean(True_region == Predicted_region) * 100
  )

overall_region_accuracy <- mean(results_df$True_region == results_df$Predicted_region) * 100

new_order <- c("NWInner", "EInner", "EOuter", "ShelfOff", "NEOFF", "DeepOff", "EastOff")

library(ggplot2)
pop_plot <- ggplot(pop_accuracy, aes(x = factor(True_pop, levels = new_order), 
                                     y = accuracy, fill = factor(True_pop, levels = new_order))) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = ordered_colors) +
  geom_text(aes(label = sprintf("%.1f%%", accuracy)), 
            vjust = -0.8, size = 3.5) +
  geom_hline(yintercept = overall_pop_accuracy, 
             linetype = "dashed", color = "red") +
  annotate("text", x = 2, y = overall_pop_accuracy, 
           label = sprintf("Overall: %.1f%%", overall_pop_accuracy),
           hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "Population", y = "Accuracy (%)") +
  ggtitle("Assignment Accuracy by Population") +
  ylim(0,100)
pop_plot

ggsave(file="figures/accuracy_pop_depth.png",pop_plot, h=4, w=5)


region_accuracy$True_region <- factor(region_accuracy$True_region, 
                                      levels=c("Inner", "ShelfOff", "Offshore"))

region_plot <- ggplot(region_accuracy, aes(x = True_region, y = accuracy, 
                                           fill = True_region )) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = ordered_colors[c(1,4,7)]) +
  geom_text(aes(label = sprintf("%.1f%%", accuracy)), 
            vjust = -0.8, size = 3.5) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "Region", y = "Accuracy (%)") +
  ggtitle("Assignment Accuracy by Region")+
  ylim(0,100)

region_plot

ggsave(file="figures/accuracy_region_depth.png",region_plot, h=4, w=5)

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




