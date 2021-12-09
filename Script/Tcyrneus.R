#

# Isaia et al., in prep.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Author: Stefano Mammola

# Loading R packages ------------------------------------------------------

library("Amelia") 
library("BAT")
library("ggplot2")
library("hypervolume")
library("psych")
library("reshape")
library("summarytools")
library("tidyverse")
library("vegan")

# Functions ---------------------------------------------------------------

source("Script/functions.r")

#################################################################
# Loading the data ----------------------------------------------
#################################################################

### loading the trait data ###
intraspecific <- read.csv("Data/db_morpho_intraspecific.csv", sep ="\t", dec = ",", as.is = FALSE)
trait         <- read.csv("Data/Database_Traits_Only_Troglohyphantes.csv", sep="\t", dec=",", as.is = FALSE) 

### Analysing Intraspecific data
str(intraspecific)

table(intraspecific$species)

# Generate traits
# see table 2 in:
# Mammola, S., Arnedo, M. A., Fišer, C., Cardoso, P., Dejanaz, A. J., & Isaia, M. (2020). 
# Environmental filtering and convergent evolution determine the ecological specialization of subterranean spiders.
# Functional Ecology, 34(5), 1064-1077

m_intra2 <- intraspecific %>% mutate(species,
                                    Pigmentation = pigment,
                                    Leg1 = rowSums(across(4:7)),
                                    Leg2 = rowSums(across(8:11)),
                                    Leg3 = rowSums(across(12:15)),
                                    Leg4 = rowSums(across(16:19)),
                                    Sternum_ratio = Sterno_L/Sterno_W,
                                    Leg_elong = rowSums(across(4:7))/Sterno_L,
                                    Cepha_red = HEIGHT/WIDTH,
                                    Anterior_eye_reg = (AME + ALE) / Ant_LINE,
                                    Posterior_eye_reg = (PLE + PME) / Post_LINE,
                                    Fang = FANG) %>% 
  dplyr::select(species,
                Pigmentation,Leg1,Leg2,Leg3,Leg4,Sternum_ratio,Leg_elong,Cepha_red,Anterior_eye_reg,Posterior_eye_reg,Fang)

m_intra <- m_intra2[,2:(ncol(m_intra2))] #remove species

# Data exploration and preparation ----------------------------------------

# Checking continous variables
continous <- c(2:ncol(m_intra)) #Id continous variables

par(mfrow= c(3,3))
for (i in continous)
  dotchart(m_intra[,i], main= colnames(m_intra)[i]) ; rm(i)

# Missing data 
Amelia::missmap(m_intra)

# Standardize continuous traits
m_intra <- BAT::standard(trait = m_intra, method = "standard", convert = continous)

# Collinearity
psych::pairs.panels(m_intra[,colnames(m_intra[continous])]) 

#Remove sternum ratio---
m_intra <- m_intra %>% dplyr::select(-c(Sternum_ratio))

# Estimating Gower distance  --------------------------------------
colnames(m_intra)
# Creating groups for traits to obtain optimization
groups_traits <-
  c(rep(1), 
    rep(2,4),
    rep(1,4),
    3)

#group 1 ---> adaptive traits
#group 2 ---> body size
#group 3 ---> diet

length(groups_traits) == ncol(m_intra) #check if groups have same length(). Should be TRUE

# DO NOT RUN (takes several minutes)
fdist <- gawdis::gawdis(data.frame(m_intra),
                        groups = groups_traits, # grouping variable defined above
                        w.type = "optimized",
                        opti.maxiter = 300,
                        ord="podani",
                        groups.weight = TRUE)

saveRDS(fdist,"functional_distance.rds") #storing the data

fdist <- readRDS("Data/functional_distance.rds") ; rm(groups_traits)


# Principal Coordinate Analysis -------------------------------------------

# PCoA
ord <- cmdscale(fdist, k = 3)
HV <- data.frame(ord, species = m_intra2$species) ; colnames(HV)[1:3] <- c("PC1", "PC2","PC3")

# PC fit
(fit12 <- vegan::envfit(ord, m_intra, choices=c(1:2), na.rm = TRUE))
(fit13 <- vegan::envfit(ord, m_intra, choices=c(1:3), na.rm = TRUE))
(fit23 <- vegan::envfit(ord, m_intra, choices=c(2:3), na.rm = TRUE))

dev.off()
par(mfrow=c(1,1), mar=c(rep(2,4)))
fit <- vegan::envfit(ord, m_intra, choices=c(1,2), na.rm = TRUE)
{
  plot(ord)
  trait_pos <- get_position(fit, add = TRUE)
}

trait_pos <- data.frame(trait_pos) %>%
  rownames_to_column("Trait")

# Plot
myCol<-viridis::viridis(n=6, option="B") #Make Color palette

(plot_traits <-   ggplot(HV, aes(PC1, PC2)) +
    stat_density_2d(
      aes(fill = ..level..),
      geom = "polygon",
      colour = NA,
      alpha = .5,
     # h = .55
    ) +
    geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
    geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
    geom_point(
      data = trait_pos,
      aes(x = Dim1, y = Dim2),
      shape = 21,
      fill = "white",
      size = 2.5
    ) +
    geom_point(
      data = trait_pos,
      aes(x = Dim1, y = Dim2),
      shape = 19,
      colour = "black",
      size = 1
    ) +
    scale_fill_gradientn(colours = rev(myCol)) +
    ggrepel::geom_text_repel(data = trait_pos, aes(x = Dim1, y = Dim2, label = Trait)) +
    theme(
      panel.background = element_rect(
        fill = NA,
        colour = "black",
        size = 1,
        linetype = "solid"
      ),
      panel.grid = element_blank(),
      legend.position = "none"
    )  +
    labs(x = "PCoA 1", y = "PCoA 2") +
    #ylim(-.46, .36) + xlim(-.46, .36) +
    coord_fixed())

# Hypervolume analysis ----------------------------------------------------

## Hypervolume construction (all species)

n.var = 3 #pcoa axes

HV$species = factor(HV$species, levels = c("cyrnaeus", levels(HV$species)[levels(HV$species) != "cyrnaeus"]))
levels(HV$species)

hv = hypervolume_gaussian(subset(HV,
                                 species == levels(HV$species)[[max(nlevels(HV$species))]])[,c(1:n.var)],
                                 name = levels(HV$species)[[max(nlevels(HV$species))]],
                                 kde.bandwidth = estimate_bandwidth(subset(HV,species==levels(HV$species)[[max(nlevels(HV$species))]])[,c(1:n.var)],method="cross-validation"))

for (i in (nlevels(factor(HV$species))-1):1) {
  
  hv2 = hypervolume_gaussian(subset(HV,species==levels(HV$species)[[i]])[,c(1:n.var)],name=levels(HV$species)[[i]],kde.bandwidth =estimate_bandwidth(subset(HV,species==levels(HV$species)[[i]])[,c(1:n.var)],method="cross-validation"))
  hv  = hypervolume_join(hv2,hv)
  
}

## Visualizing the hypervolumes

levels(HV$species)
# Pair plot
color_custom = c("red","orange","blue","orange","purple","purple","blue","blue","black","black","orange","blue","black","orange")

plot(hv,
     num.points.max.random = 2000,
     colors=color_custom,
     show.data=F,
     pairplot=T,
     show.3d=F,
     showdensity=T,
     showrandom=F,
     contour.lwd=0.5,
     cex.names=0.8,
     cex.random=0.3,
     cex.legend=1,
     cex.axis=1,
     contour.filled=T,
     contour.filled.alpha=0.2,
     show.centroid=TRUE,
     cex.centroid=2)

volume <- BAT::kernel.alpha(hv)*1000
beta_diversity <- BAT::kernel.beta(hv)
centroid <- as.matrix(BAT::kernel.similarity(hv)$Distance_centroids)

overlap_matrix <- as.matrix(beta_diversity$Btotal)

overlap_matrix[upper.tri(overlap_matrix)] <- centroid[upper.tri(centroid)] #replace above the diagonal with centroid
diag(overlap_matrix) <- volume #replace diagonal with volume

write.table(round(overlap_matrix,2), "Table_X.csv", sep = "\t")

# Mapping traits in relation to all Troglohyphantes  ----------------------

trait_m3 <- trait %>% dplyr::select(
    Genus_species,
    Ecological_classification,
    Leg_elongation = Femur_I_length_avg/Prosoma_length_avg,
    Body_length_avg,
    Eye_reduction_ratio)

trait_m2 <- trait_m3[,3:ncol(trait_m3)]

# Data exploration and preparation ----------------------------------------

# Checking continous variables
continous <- c(1:ncol(trait_m2)) #Id continous variables
dev.off()
par(mfrow= c(2,2), mar = c(rep(2,4)))
for (i in continous)
  dotchart(trait_m2[,i], main= colnames(trait_m2)[i]) ; rm(i)
 
# Missing data 
Amelia::missmap(trait_m2)

# Collinearity
psych::pairs.panels(trait_m2[,colnames(trait_m2[continous])]) 

# Plot
db_ggplot2 <- data.frame(Species = rep(trait_m3$Genus_species,3),
  
                         Value = c( trait_m2$Body_length_avg, 
                                    trait_m2$Leg_elongation,
                                    trait_m2$Eye_reduction_ratio),
                         
                         Trait = c(rep("Body length",   nrow(trait_m2)),
                                   rep("Femur elongation",nrow(trait_m2)),
                                   rep("Eye development",nrow(trait_m2))))

db_ggplot2 %>% ggplot(aes(y=Value, x =1)) + facet_wrap( ~ Trait, nrow = 1, ncol = 3, scale ="free") +
  geom_violin(alpha =0.25,outlier.alpha = 0,width=0.2, col = "grey70")+
  labs(xlab = "",ylab = "")+
  geom_point(size = 3, alpha = 0.15) + 
  geom_point(data= db_ggplot2[db_ggplot2$Species == "Troglohyphantes cyrnaeus",], 
             aes(y=Value, x =1), size = 5, shape = 21,col = "black",fill = "red", alpha = 0.9) + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
