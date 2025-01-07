# prep genetic data

library(adegenet)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(algatr)

# dolphin seascape genetics

# need to find some measure of genetic distance overall between indivs/pops
# could this be pairwise fst? pc space?
# examples:
    # https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4745
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7463313/
    # https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
    # https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/4.Redundancy_analysis/4.redundancy_analysis.R
    # seems very relevant: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9873513/
# https://www.int-res.com/abstracts/meps/v477/p107-121/
  # - use linearized fst. 
# compairson of individual based genetic distances: https://doi.org/10.1111/1755-0998.12684
  # . We calculated the PCA-based GD metrics in the R statistical environment 
      # by first calculating principal components (PC) from allele
      # usage (0, 1 or 2) for all alleles in the population and then creating
      # distance matrices from the Euclidean distance among varying numbers of PC axes (1, 4, 16 or 64; Shirk et al., 2010). 
# also suggests pca is good: https://onlinelibrary.wiley.com/doi/pdf/10.1111/1755-0998.13831

# how exactly to summarize the variation across the PCs?

# generalized dissimilarity modelling (GDM; Ferrier, 2002; Ferrier et al., 2002, 2007) \
  # multiple matrix regression with randomization (MMRR; Wang, 2013; Figure 1).
  # MMRR performs linear matrix regression on genetic and environmental distances 
      # and allows for multiple independent variables (environmental and geographic distances) 
      # to be examined simultaneously
      # gives a p-value for each explanatory variable and the whole model

# really good pipeline and paper, for genomics but many things apply:
  # https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13884

# I think I need to get everything to population. 
  # Then summarize environment. for each. 
  # then do the GEA. could run GLM? but maybe RDA makes sense. 


df <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_MicroSats-adegenet.csv",
         row.names=1)

genin <- df2genind(df, sep="/")


#pca:

pca1 <- dudi.pca(genin,scannf=FALSE,scale=FALSE)
## basic plot
plot(pca1$li, pch=21, bg="grey45", color="black",  cex=1)


# how much variation is explained?
barplot((100*(pca1$eig/sum(pca1$eig))[1:20]))

distgenEUCL <- adegenet::dist(genin, method = "euclidean", 
                    diag = FALSE, upper = FALSE, p = 2)
hist(distgenEUCL)

# snps
df <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_SNPs-adegenet.csv",
               row.names=1)
head(df)

# change snps to correct format: each locus one entry with #/#
# note that microsats were already in this format.
out <- as.data.frame(matrix(ncol=ncol(df)/2, nrow=nrow(df)))

df_tracker <- 1
outnm <- NA
for(i in 1:ncol(out)){
  out[,i] <- paste(df[,df_tracker], df[,df_tracker+1], sep="/")  
  outnm[i] <- colnames(df)[df_tracker]
  df_tracker <- df_tracker + 2
}

colnames(out) <- gsub("\\.","_",outnm)
genin <- df2genind(out, sep="/")


#pca:

pca1 <- dudi.pca(genin,scannf=FALSE,scale=FALSE)
## basic plot
plot(pca1$li, pch=21, bg="grey45",  cex=1)

# read in the pop labels.
dinfo <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

head(pca1$li)

# how much variation is explained?
barplot((100*(pca1$eig/sum(pca1$eig))[1:20]))



#-----------------------------------------------------------------------------
# combined snps and microsats
# read in micros
df2 <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_MicroSats-adegenet.csv",
               row.names=1)
# merge the micros and snps
dfall <- cbind(out, df2)
genin <- df2genind(dfall, sep="/")

pca1 <- dudi.pca(genin,scannf=FALSE,scale=FALSE)
## basic plot
plot(pca1$li, pch=21, bg="grey45",  cex=1)

# read in the pop labels.
dinfo <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

head(pca1$li)

# how much variation is explained?
barplot((100*(pca1$eig/sum(pca1$eig))[1:20]))

dinfo$Pop <- as.factor(dinfo$Pop)

Colorsdf <-
  with(dinfo,
       data.frame(population = levels(dinfo$Pop),
                  color = I(brewer.pal(nlevels(Pop), name = 'Dark2'))))
cols <- Colorsdf$color[match(dinfo$Pop, Colorsdf$population)]

pltdat <- cbind(dinfo$Pop, pca1$li)
colnames(pltdat) <- c("Pop", "PC1", "PC2")
pout <-ggplot(pltdat, aes(x=PC1, y=PC2, fill=Pop),
       color="black")+
  geom_point(size=3,pch=21) +
  theme_classic(base_size=14) +
  scale_fill_manual(values = Colorsdf$color,guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  theme(legend.title = element_blank())

ggsave("figures/pca_pops.png", pout, h=4, w=5)
ggsave("figures/pca_pops.pdf", pout, h=4, w=5)





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# get groups based on DAPC
# and number of pcs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# use merged df
df <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_SNPs-adegenet.csv",
               row.names=1)
head(df)

out <- as.data.frame(matrix(ncol=ncol(df)/2, nrow=nrow(df)))

df_tracker <- 1
outnm <- NA
for(i in 1:ncol(out)){
  out[,i] <- paste(df[,df_tracker], df[,df_tracker+1], sep="/")  
  outnm[i] <- colnames(df)[df_tracker]
  df_tracker <- df_tracker + 2
}

colnames(out) <- gsub("\\.","_",outnm)
genin <- df2genind(out, sep="/")

df2 <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_MicroSats-adegenet.csv",
                row.names=1)
dfall <- cbind(out, df2)
genin <- df2genind(dfall, sep="/")

#grp <- find.clusters(genin, max.n.clust=40, n.pca=7)
grp <- find.clusters(genin, max.n.clust=40, n.pca=200)
# keep 6
dapc1 <- dapc(genin, grp$grp, n.pca=100, n.da=50)


#table(pop(x), grp$grp)

xval <- xvalDapc(genin, grp=grp$grp, n.pca.max = 100, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc2 <- dapc(genin,grp$grp, n.da=50, n.pca=100)

tmp <- optim.a.score(dapc2)
# 13 pcs

# use k-1
dapc3 <- dapc(genin,grp$grp, n.da=5, n.pca=13)
myCol <- rainbow(15)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc3, lab="", posi=list(x=12,y=-.01), cleg=.7)


round(head(dapc3$posterior),3)
summary(dapc3 )

plot(dapc3$posterior, dapc1$posterior)
assignplot(dapc3, subset=1:50)

compoplot(dapc3, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))

popgroup <- as.data.frame(matrix(ncol=2, nrow=nrow(dapc3$posterior)))
colnames(popgroup) <- c("indiv", "pop")
popgroup$indiv <- row.names(dapc3$posterior)
popgroup$pop[dapc3$posterior[,1]> 0.6] <- "pop_1"
popgroup$pop[dapc3$posterior[,2]> 0.6] <- "pop_2"
popgroup$pop[dapc3$posterior[,3]> 0.6] <- "pop_3"
popgroup$pop[dapc3$posterior[,4]> 0.6] <- "pop_4"
popgroup$pop[dapc3$posterior[,5]> 0.6] <- "pop_5"
popgroup$pop[dapc3$posterior[,6]> 0.6] <- "pop_6"

popgroup$pred.post <- paste0("pop_",apply(dapc3$posterior, 1, which.max))

# figure out if these group assignments correspond to the "pops" in the data:

popgroup$geo_pop <-dinfo$Pop
head(popgroup)
popgroup$comp <- paste(popgroup$pred.post, popgroup$geo_pop, sep=":")
table(popgroup$comp)

#sum(is.na(popgroup$pop))

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black")
scatter(dapc3,
        scree.da=TRUE,  # Shows scree plot of eigenvalues
        bg="white",
        pch=20,         
        col=colors,     
        #legend=TRUE,   
        solid=0.9)      

# color PCA by new pops:

table(popgroup$pop)
sum(is.na(popgroup$pop))

popgroup$newpop <- NA
popgroup$newpop[popgroup$pred.post == "pop_1"] <- "NW_InnerShelf"
popgroup$newpop[popgroup$pred.post == "pop_2"] <- "NW_Outer_shelf"
popgroup$newpop[popgroup$pred.post == "pop_3"] <- "East_Oceanic"
popgroup$newpop[popgroup$pred.post == "pop_4"] <- "NW_Oceanic"
popgroup$newpop[popgroup$pred.post == "pop_5"] <- "NE_Oceanic"
popgroup$newpop[popgroup$pred.post == "pop_6"] <- "East_shelf"

table(popgroup$newpop)

pltdat <- cbind(popgroup$comp, popgroup$newpop, popgroup$geo_pop, pca1$li)
colnames(pltdat) <- c("Pop_comparison","DAPC_pop", "geo_pop", "PC1", "PC2")

p1 <- ggplot(pltdat, aes(x=PC1, y=PC2, color=geo_pop))+
  geom_point(size=3) +
  ggtitle("geo_populations")

p1

#pltdat <- cbind(popgroup$pop, pca1$li)
#colnames(pltdat) <- c("Pop", "PC1", "PC2")

p2 <- ggplot(pltdat, aes(x=PC1, y=PC2, color=DAPC_pop))+
  geom_point(size=3) +
  ggtitle("dapc_populations")

ggsave(ggpubr::ggarrange(p1, p2), filename="figures/pca_comp.png", h=5, w=9)



write.csv(popgroup, file="analysis/population_assignments.csv", row.names=F)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# pca genetic distance matrix, based on tutorial pipeline https://github.com/TheWangLab/algatr/blob/558519684f5d346ca295f68b184c7ed99cf9c97e/R/Gen_dist.R#L17
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# base decisions on DAPC, above.


#gl <- as.genlight(dfall)
#gl <- as.genlight(genin)
#df2 <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore_MicroSats-adegenet.csv",
#                row.names=1)
#dfall <- cbind(out, df2)
genin <- df2genind(dfall, sep="/")

pc <- dudi.pca(genin,scale=FALSE,nf=80,scannf = FALSE )

eig <- pc$eig
# Run Tracy-Widom test
# to determine how many pcs to include in the distance calculation
# NOTE: critical point corresponds to significance level.
# If the significance level is 0.05, 0.01, 0.005, or 0.001,
# the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively.
# The default is 2.0234.
tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 2.0234)
npc <- tw_result$SigntEigenL
npc

# Calculate PC-based distance
dists_20pcs <- as.matrix(dist(pc$li[, 1:npc], diag = TRUE, upper = TRUE))

# repeat with screen plot results to see how consistent
# 6 pcs, bc 7 putative pops (K-1)
dists_3pcs <- as.matrix(dist(pc$li[, 1:13], diag = TRUE, upper = TRUE))

dist_x <- as.data.frame(dists_20pcs)
dist_y <- as.data.frame(dists_3pcs)

# Check to ensure sample IDs match ----------------------------------------
if (all(rownames(dist_x) == rownames(dist_y)) == FALSE) {
  stop("Sample IDs do not match")
}

# Melt data from square to long -------------------------------------------
# Assign NAs to upper triangle of square matrix
dist_x1 <- as.matrix(dist_x)
dist_x[upper.tri(dist_x, diag = FALSE)] <- NA

melt_x <- dist_x %>%
  tibble::rownames_to_column(var = "comparison") %>%
  tidyr::pivot_longer(cols = -(comparison)) %>%
  na.omit() %>%
  dplyr::filter(comparison != name)%>%
  dplyr::rename(pcs20 = value)

# plot results
gen_dist_hm(dist_x)

# now do k-1 pcs
dist_y[upper.tri(dist_y, diag = FALSE)] <- NA

melt_y <- dist_y %>%
  tibble::rownames_to_column(var = "comparison") %>%
  tidyr::pivot_longer(cols = -(comparison)) %>%
  na.omit() %>%
  dplyr::filter(comparison != name) %>%
  dplyr::rename(pcs3 = value)


# compare the two approaches:
joined <- dplyr::full_join(melt_x, melt_y)
joined %>%
  ggplot2::ggplot(ggplot2::aes(x = pcs20, y = pcs3)) +
  ggplot2::geom_abline(ggplot2::aes(intercept = 0.0, slope = 1), color = "gray") +
  ggplot2::geom_point(color = "black", size = .2, alpha = .2)+
  geom_smooth(method=lm)


# k-1 underestimates relative to using 20
dist_y1 <- as.matrix(dist_y)

write.csv(melt_y, file="genetic_distances_pairwise_indivs.csv", quote=F, row.names=F)
write.csv(dist_y1, file="genetic_distances_pairwise_matrix.csv", quote=F, row.names=F)



#-----------------------------------------------------------------
# nj tree. genetic distances
#-----------------------------------------------------------------
library(ape)
library(Matrix)
library(tidyverse)
library(RColorBrewer)

gendist <- as.matrix(read.csv("genetic_distances_pairwise_matrix.csv"))
gendist2 <- as.matrix(Matrix::forceSymmetric(gendist,uplo="L"))


all <- read.csv("population_assignments.csv")

all$newpop <- as.factor(all$newpop)

my_colors <- c("#1B9E77","firebrick3","#7570B3","#E6AB02","#8B4513","black")

Colorsdf <-
  with(all,
       data.frame(population = levels(all$newpop),
                  color = I(my_colors)))
cols <- Colorsdf$color[match(all$newpop, Colorsdf$population)]


nj(gendist2) %>% plot(.,"unrooted", tip.color = cols)
nj(gendist2) %>% plot(., tip.color = cols)

legend("bottomright", 
       legend = Colorsdf$population,
       col = Colorsdf$color,
       pch = 19,
       title = "Population")



library(ggtree)
out_tree <- nj(gendist2)
out <- as_tibble(nj(gendist2))
out$label <- gsub("X", "", out$label)
dm <- left_join(out, all, by=c("label" ="indiv"))

ggtree(out_tree) + 
  theme_tree()

#p <- ggtree(out_tree,layout="daylight") + theme_tree()
#p <- ggtree(out_tree, layout="unrooted") + theme_tree()
p <- ggtree(out_tree, layout="ape") + theme_tree()


pout <- p %<+% dm + 
  #geom_tiplab(aes(color=newpop), size=0.9) +
  theme(legend.position="right")+ 
  #geom_text( show.legend  = F ) +
  geom_tippoint(aes(color=newpop), size=4, alpha=0.7) +
  scale_color_manual(values=my_colors)

pout
ggsave(pout, filename="figures/nj_tree_unrooted.png", h=5, w=6)
ggsave(pout, filename="figures/nj_tree_unrooted.pdf", h=5, w=6)






#------------------------------------------------------------------------
# map with new pops;
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

dat$new_pops <- popgroup$newpop
# 
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = dat, aes(x = long, y = lat, fill=new_pops, shape=new_pops), 
             shape= 21, color="black", size =3,
             alpha=1) +
  coord_sf() +
  theme_bw(base_size = 14) +
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
  annotation_scale() +
  scale_fill_manual(values = my_colors)

p

ggsave("figures/map_newpops.pdf", p, h=5, w=7)
ggsave("figures/map_newpops.png", p, h=5, w=7)
