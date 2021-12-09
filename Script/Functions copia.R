## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

#Custom function for trait position
get_position<-function (x, choices = c(1, 2), labels, arrow.mul, at = c(0, 0), axis = FALSE, p.max = NULL, col = "blue", bg, add = TRUE, 
                        ...) {
  formals(arrows) <- c(formals(arrows), alist(... = ))
  labs <- list(v = rownames(x$vectors$arrows), f = rownames(x$factors$centroids))
  if (!missing(labels)) {
    if (is.list(labels)) {
      if (!is.null(labs$v) && !is.null(labels$vectors)) 
        labs$v <- labels$vectors
      if (!is.null(labs$f) && !is.null(labels$factors)) 
        labs$f <- labels$factors
    }
    else {
      if (!is.null(x$vectors) && !is.null(x$factors)) 
        stop("needs a list with both 'vectors' and 'factors' labels")
      if (!is.null(x$factors)) 
        labs$f <- labels
      else labs$v <- labels
    }
  }
  vect <- NULL
  if (!is.null(p.max)) {
    if (!is.null(x$vectors)) {
      take <- x$vectors$pvals <= p.max
      x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
      labs$v <- labs$v[take]
      x$vectors$r <- x$vectors$r[take]
      if (nrow(x$vectors$arrows) == 0) 
        x$vectors <- vect <- NULL
    }
    if (!is.null(x$factors)) {
      tmp <- x$factors$pvals <= p.max
      nam <- names(tmp)[tmp]
      take <- x$factors$var.id %in% nam
      x$factors$centroids <- x$factors$centroids[take, 
                                                 , drop = FALSE]
      labs$f <- labs$f[take]
      if (nrow(x$factors$centroids) == 0) 
        x$factors <- NULL
    }
  }
  if (!is.null(x$vectors)) {
    vect <- sqrt(x$vectors$r) * x$vectors$arrows[, choices, 
                                                 drop = FALSE]
    if (missing(arrow.mul)) {
      if (!add) 
        arrow.mul <- 1
      else arrow.mul <- ordiArrowMul(vect, at = at)
    }
    if (axis) {
      maxarr <- round(sqrt(max(x$vectors$r)), 1)
      ax <- -c(-1, 0, 1) * arrow.mul * maxarr
    }
    vect <- arrow.mul * vect
    vect <- sweep(vect, 2, at, "+")
    if (add) {
      vtext <- ordiArrowTextXY(vect, labs$v, rescale = FALSE, 
                               at = at)
    }
  }
  if (!add) {
    plot.new()
    if (is.null(vect) || is.null(x$factors)) {
      xstack <- rbind(vect, x$factors$centroids)
      plot.window(xlim = range(xstack[, 1], at[1]), ylim = range(xstack[, 
                                                                        2], at[2]), asp = 1, ...)
    }
    else {
      plot.window(xlim = range(x$factors$centroids[, 1], 
                               at[1]), ylim = range(x$factors$centroids[, 2], 
                                                    at[2]), asp = 1)
      vfill <- 0.75
      arrow.mul <- ordiArrowMul(vect, at = at, fill = 1)
      vect <- arrow.mul * vect
    }
    sw <- strwidth(c(labs$v, labs$f), ...)/2
    sh <- strheight(c(labs$v, labs$f), ...)
    xstack <- rbind(x$factors$centroids, vect)
    xlim <- range(xstack[, 1] + sw, xstack[, 2] - sw)
    ylim <- range(xstack[, 2] + sh, xstack[, 2] - sh)
    plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
    if (!is.null(vect)) {
      arrow.mul <- ordiArrowMul(vect, at = at, fill = 1)
      vect <- arrow.mul * vect
      vtext <- ordiArrowTextXY(vect, labs$v, at = at, 
                               rescale = FALSE, ...)
      sw <- strwidth(labs$v, ...)/2
      sh <- strheight(labs$v, ...)
      xlim <- range(xlim, vtext[, 1] + sw, vtext[, 1] - 
                      sw)
      ylim <- range(xlim, vtext[, 2] + sh, vtext[, 2] - 
                      sh)
      plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
    }
    axis(side = 1, ...)
    axis(side = 2, ...)
    box(...)
    alabs <- colnames(vect)
    title(..., ylab = alabs[2], xlab = alabs[1])
  }
  if (!is.null(vect)) {
    arrows(at[1], at[2], vect[, 1], vect[, 2], len = 0.05, 
           col = col)
    if (missing(bg)) 
      text(vtext, labs$v, col = col, ...)
    else ordilabel(vtext, labels = labs$v, col = col, fill = bg, 
                   ...)
  }
  if (!is.null(x$factors)) {
    if (missing(bg)) 
      text(x$factors$centroids[, choices, drop = FALSE], 
           labs$f, col = col, ...)
    else ordilabel(x$factors$centroids[, choices, drop = FALSE], 
                   labels = labs$f, col = col, fill = bg, ...)
  }
  if (axis && !is.null(vect)) {
    axis(3, at = ax + at[1], labels = c(maxarr, 0, maxarr), 
         col = col)
    axis(4, at = ax + at[2], labels = c(maxarr, 0, maxarr), 
         col = col)
  }
  return(rbind(vect,x$factors$centroids))
  
}

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)

# Setting working directory -----------------------------------------------

setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/Isaia_et_al_Troglohyphantes_cyrneus") #change with your directory

# Loading R packages ------------------------------------------------------

library("Amelia") 
library("BAT")
library("ggplot2")
library("hypervolume")
library("psych")
library("summarytools")
library("tidyverse")
library("vegan")

#################################################################
# Loading the data ----------------------------------------------
#################################################################

### loading the trait data ###
intraspecific <- read.csv("db_morpho_intraspecific.csv", sep ="\t", dec = ",", as.is = FALSE)
trait  <- read.csv("Database_Traits_Only_Troglohyphantes.csv", sep="\t", dec=",", as.is = FALSE) 
#ita    <- read.csv("troglo_italy.csv", sep="\t", dec=",", as.is = FALSE) 

### Analysing Intraspecific data
str(intraspecific)

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
fdist <- readRDS("functional_distance.rds") ; rm(groups_traits)

# PCoA
ord <- cmdscale(fdist, k = 3)
HV <- data.frame(ord, species = m_intra2$species)

colnames(HV)[1:3] <- c("PC1", "PC2","PC3")

dev.off()
par(mfrow=c(1,1), mar=c(rep(2,4)))
fit <- vegan::envfit(ord, m_intra, na.rm = TRUE)
{
  plot(ord)
  trait_pos <- get_position(fit, add = TRUE)
}

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

library("reshape")
volume <- BAT::kernel.alpha(hv)*1000
beta_diversity <- BAT::kernel.beta(hv)
centroid <- as.matrix(BAT::kernel.similarity(hv)$Distance_centroids)

overlap_matrix <- as.matrix(beta_diversity$Btotal)

overlap_matrix[upper.tri(overlap_matrix)] <- centroid[upper.tri(centroid)] #replace above the diagonal with centroid
diag(overlap_matrix) <- volume #replace diagonal with volume

write.table(round(overlap_matrix,2), "Table_X.csv", sep = "\t")

# Mapping traits in relation to all Troglohyphantes  ----------------------

trait_m3 <- trait %>% select(
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
                                   rep("Eye reduction (ratio)",nrow(trait_m2))))


db_ggplot2 %>% ggplot(aes(y=Value, x =1)) + facet_wrap( ~ Trait, nrow = 1, ncol = 3, scale ="free") +
  geom_violin(alpha =0.25,outlier.alpha = 0,width=0.2, col = "grey70")+
  labs(xlab = "",ylab = "")+
  geom_point(size = 3, alpha = 0.15) + 
  geom_point(data= db_ggplot2[db_ggplot2$Species == "Troglohyphantes cyrnaeus",], 
             aes(y=Value, x =1), size = 3, col = "red",fill = "red", alpha = 0.9) + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




# 
# # Estimating Gower distance  --------------------------------------
# 
# str(trait_m2)
# # Creating groups for traits to obtain optimization
# groups_traits <-
#   c(rep(1,9), 
#     rep(2,5))
# 
# #group 1 ---> adaptive traits
# #group 2 ---> body size
# 
# length(groups_traits) == ncol(trait_m2) #check if groups have same length(). Should be TRUE
# 
# # DO NOT RUN (takes several minutes)
# fdist <- gawdis::gawdis(data.frame(trait_m2),
#                 groups = groups_traits, # grouping variable defined above
#                  w.type = "optimized",
#                  opti.maxiter = 300,
#                  ord="podani",
#                  groups.weight = TRUE)
# 
# 
# trait_m2 <- trait_m3[3:ncol(trait_m3)]
# continous <- c(4:ncol(trait_m2)) #Id continous variables
# trait_m2 <- BAT::standard(trait = trait_m2, method = "standard", convert = continous)
# 
# fdist <- gawdis::gawdis(data.frame(trait_m2),
#                         #groups = groups_traits, # grouping variable defined above
#                         w.type = "optimized",
#                         opti.maxiter = 300,
#                         ord="podani",
#                         groups.weight = TRUE)
# 
# saveRDS(fdist,"functional_distance.rds") #storing the data
# fdist <- readRDS("functional_distance.rds") ; rm(groups_traits)
# 
# #principal coordinate analyses
# ord <- cmdscale(fdist, k = 2)
# coordinates <- data.frame(ord, Genus_species = trait_m3$Genus_species,Ecological_classification = trait_m3$Ecological_classification)
# colnames(coordinates)[1:2] <- c("PC1", "PC2")
# 
# dev.off()
# par(mfrow=c(1,1), mar=c(rep(2,4)))
# fit <- vegan::envfit(ord, trait_m2, na.rm = TRUE)
# {
#   plot(ord)
#   trait_pos <- get_position(fit, add = TRUE)
# }
# 
# trait_pos <- data.frame(trait_position) %>%
#   rownames_to_column("Trait")
# 
# ## PLOT
# 
myCol<-viridis::viridis(n=6, option="B") #Make Color palette
pie(seq_along(myCol),myCol,col= myCol)

ita$Genus_species <- paste(ita$Genus, ita$species)
coordinates <- dplyr::left_join(coordinates, ita, by = "Genus_species")
coordinates_ita <- na.omit(coordinates)

levels(coordinates_ita$CAT_ECOL) <- c("High","Intermediate","Low")

(plot_space_trait <-  ggplot(data = coordinates[coordinates$Ecological_classification == "Troglobiont",], aes(x=PC1, y=PC2)) +
  # stat_density_2d(
  #   aes(fill = ..level..),
  #   geom = "polygon",
  #   colour = NA,
  #   alpha = .5,
  #   h = .25
  # ) +
  geom_hline(aes(yintercept = 0), linetype = 3, colour = "gray70") +
  geom_vline(aes(xintercept = 0), linetype = 3, colour = "gray70") +
  geom_point(shape = 21,
             size = 2,
             fill = "white",
             colour = "grey20") +

  geom_point(data= coordinates_ita, aes(fill = CAT_ECOL), shape = 21,
               size = 2,
               colour = "grey20") +

  geom_point(data = coordinates[coordinates$Genus_species=="Troglohyphantes cyrnaeus",],
             aes(PC1, PC2),
             shape = 21,
               size = 3,
               fill = "red",
               colour = "black") +

  geom_point(data = coordinates[coordinates$Genus_species=="Troglohyphantes cyrnaeus",],
               aes(PC1, PC2),
               shape = 16,
               size = 2,
               fill = "black",
               colour = "black") +

    geom_point(data = coordinates[coordinates$Genus_species=="Troglohyphantes cyrnaeus",],
               aes(PC1, PC2),
               shape = 16,
               size = 1,
               fill = "red",
               colour = "black") +
  #
  #scale_fill_gradientn(colours = rev(myCol)) +
  scale_fill_manual(values=c("blue","orange","grey30")) +
  theme(
    panel.background = element_rect(
      fill = NA,
      colour = "black",
      size = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    #legend.position = "none"
  ) +
  labs(x = "PCoA 1", y = "PCoA 2", fill = "Subterranean\nadaptation") +
  #ylim(-.46, .36) + xlim(-.46, .36) +
  coord_fixed())



(plot_traits <-   ggplot(coordinates, aes(PC1, PC2)) +
  # stat_density_2d(
  #   aes(fill = ..level..),
  #   geom = "polygon",
  #   colour = NA,
  #   alpha = .5,
  #   h = .25
  # ) +
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
  #scale_fill_gradientn(colours = rev(myCol)) +
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

# #################################################################
# # Calculating functional diversity ------------------------------
# #################################################################
# 
# # Generating full hypervolume --------------------------------------------
# hv_tot <- BAT::kernel.build(comm = rep(1, nrow(trait_m)), trait = trait_m[,1:2],
#                           method = "gaussian", axes = 0, abund = FALSE, cores = 1)
# 
# #rownames(hv_tot@Data) <- trait$Genus_species
# 
# dist(rbind(hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes apenninicus",1:2],
#            hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes lucifuga",1:2],
#            hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes lucifer",1:2],
#            hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes fagei",1:2],
#            hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes cerberus",1:2],
#            hv_tot@Data[rownames(hv_tot@Data) == "Troglohyphantes pedemontanus",1:2]))
#       
# # Visualize
# plot(hv_tot, 
#      show.data=TRUE,
#      show.random=TRUE,
#      show.centroid=FALSE,
#      col="black", 
#      num.points.max.random = 3000,
#       plot.function.additional=function(i,j){
#        points( x= trait_m[trait_m3$Ecological_classification=="Troglobiont", i], y=trait_m[trait_m3$Ecological_classification =="Troglobiont", j],cex=1,pch=16, col ="blue")
#        points( x= trait_m[trait_m3$Ecological_classification=="Troglophile", i], y=trait_m[trait_m3$Ecological_classification =="Troglophile", j],cex=1,pch=16, col ="orange")
#        
#        points( x= trait_m[trait_m3$Genus_species=="Troglohyphantes cyrneus", i], y=trait_m[trait_m3$Genus_species=="Troglohyphantes cyrneus", j],cex=2,pch=16, col ="black")
#        points( x= trait_m[trait_m3$Genus_species=="Troglohyphantes cyrneus", i], y=trait_m[trait_m3$Genus_species=="Troglohyphantes cyrneus", j],cex=1,pch=16, col ="red")
#        
#        # points( x= data[, i], y=data[, j],cex=1.2,pch=16,col="black")
#        # 
#        # 
#        # points( x= data[data$Family=="Agelenidae", i], y=data[data$Family =="Agelenidae", j],cex=1,pch=16,col="brown")
#        # points( x= data[data$Family=="Linyphiidae", i], y=data[data$Family =="Linyphiidae", j],cex=1,pch=16,col="red")
#        # #points( x= data[data$Family=="Dysderidae", i], y=data[data$Family =="Dysderidae", j],cex=1,pch=16,col="orange")
#        # 
#        # points( x= data[data$Family=="Pholcidae", i], y=data[data$Family =="Pholcidae", j],cex=1,pch=16,col="orange")
#        # 
#        # points( x= data[data$Family=="Leptonetidae", i], y=data[data$Family =="Leptonetidae", j],cex=1,pch=16,col="blue")
#        # points( x= data[data$Family=="Nesticidae", i], y=data[data$Family =="Nesticidae", j],cex=1,pch=16,col="turquoise")
#        # 
#        # 
#        # points( x= data[data$Family=="Tetragnathidae", i], y=data[data$Family =="Tetragnathidae", j],cex=1,pch=16,col="grey")
#        # 
#        # points( x= data[data$Family=="Pimoidae", i], y=data[data$Family =="Pimoidae", j],cex=1,pch=16,col="yellow")
#        # 
#        # points( x= data[data$Family=="Symphytognathidae", i], y=data[data$Family =="Symphytognathidae", j],cex=1,pch=16,col="purple")
# 
#       }
# )
# 
# 
# # test
# db_test <- data.frame(hv_tot@Data, Genus_species = trait$Genus_species)
# db_test <- rbind(db_test[6,],db_test[-6,])
# dist <- as.matrix(dist(db_test[,1:2]))[1,]
# db_test <- data.frame(db_test,dist)
# 
# db_test %>% ggplot(aes(x = dist)) %>% geom_boxplot()
# 
# 
# 
# # test
# db_test <- data.frame(hv_tot@Data, Genus_species = trait$Genus_species)
# ita$Genus_species <- paste(ita$Genus, ita$species)
# db_test <- dplyr::left_join(db_test, ita, by = "Genus_species")
# db_test <- na.omit(db_test)
# 
# dist <- as.matrix(dist(db_test[,1:2]))[1,]
# 
# data.frame(db_test,dist)
# 
# levels(db_test$CAT_ECOL) <- c("red","blue","black")
# 
# par(mfrow=c(1,1), mar=c(rep(6,4)))
# 
# plot(db_test$DISTR,dist, col = as.character(db_test$CAT_ECOL), pch = 16, xlab = "range size", ylab = "Distance from T. albopctus")
# abline(lm(dist ~ db_test$DISTR), col = "black", lwd = 2)
# 
# 
# trait_m3$Leg_elongation
# 
# trait_m3 %>% ggplot(aes(y=Leg_elongation, x =1)) + 
#   geom_violin(alpha =0.25,outlier.alpha = 0,width=0.2, col = "grey70")+
#   labs(xlab = "")+
#   geom_point(size = 3, alpha = 0.15) + 
#   geom_point(data= trait_m3[trait_m3$Genus_species == "Troglohyphantes cyrnaeus",], 
#              aes(y=Leg_elongation, x =1), size = 3, col = "red",fill = "red", alpha = 0.9) + theme_classic()
# max(trait_m3$Leg_elongation,na.rm =TRUE)
# 
# 
# trait_m3 %>% ggplot(aes(y=AME, x =1)) + 
#   geom_violin(alpha =0.25,outlier.alpha = 0,width=0.2, col = "grey70")+
#   labs(xlab = "")+
#   geom_point(size = 3, alpha = 0.15) + 
#   geom_point(data= trait_m3[trait_m3$Genus_species == "Troglohyphantes cyrnaeus",], 
#              aes(y=AME, x =1), size = 3, col = "red",fill = "red", alpha = 0.9) + theme_classic()
# 
# 
# max(trait_m3$Leg_elongation,na.rm =TRUE)
# 
# trait_m3[trait_m3$Leg_elongation == max(trait_m3$Leg_elongation,na.rm =TRUE),]
