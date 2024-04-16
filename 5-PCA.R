#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.5 - PCA whole skull, rostrum and braincase

#LOAD LIBRARIES ----
#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(RColorBrewer) 
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)
library(phytools)
library(car)
library(Rvcg)
library(scales)

#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#PCA COMPLETE DATASET ----

##Whole skull ----
#Run PCA on complete dataset
PCA_all <- gm.prcomp(gdf$coords)

#List of PC components and proportion of variation
PCA_all 

#Save PCA results to file
sink("Output/5-PCA/PCA_all_components.txt")
print("PCA complete dataset")
print(PCA_all)
sink() 

#Change row names to codes to make plot readable
row.names(PCA_all$x) <- gdf$Id

##View plot
plot(PCA_all, main = "PCA all data - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all$x[,1], y = PCA_all$x[,2], labels = rownames(PCA_all$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

##View plot
plot(PCA_all, axis1 = 1, axis2 = 3, main = "PCA all data - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "orange",   #border of points
     bg = "orange",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all$x[,1], y = PCA_all$x[,3], labels = rownames(PCA_all$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all <- PCA_all$x 

#Save shapes of extremes for axes used in plot
PC1min_all <- PCA_all[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all <- PCA_all[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all <- PCA_all[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all <- PCA_all[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_all_points <- spheres3d(PC1min_all, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA/PC1min_all.png")
rgl.snapshot(filename = "Output/5-PCA/PC1min_all1.png") 
clear3d()

#PC1max colors
PC1max_all_points <- spheres3d(PC1max_all, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA/PC1max_all1.png") 
rgl.snapshot(filename = "Output/5-PCA/PC1max_all.png") 
clear3d()

#PC2min colors
PC2min_all_points <- spheres3d(PC2min_all, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA/PC2min_all.png") 
rgl.snapshot(filename = "Output/5-PCA/PC2min_all1.png") 
clear3d()

#PC2max colors
PC2max_all_points <- spheres3d(PC2max_all, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA/PC2max_all1.png") 
rgl.snapshot(filename = "Output/5-PCA/PC2max_all.png")  
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all_df <- as_tibble(pcscores_all)
#Add labels and other attributes to tibble as columns
pcscores_all_df <- pcscores_all_df %>% mutate(specimens = gdf$Id, group = gdf$group, category = gdf$category,
                                              genus = gdf$genus, family = gdf$family, size = gdf$size, 
                                              module = rep("skull", length(gdf$Id)))
glimpse(pcscores_all_df)

#Nice PCA plot with stages and groups
PCA_all_ggplot <- ggplot(pcscores_all_df, aes(x = Comp1, y = Comp2, label = specimens, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all$sdev[1]^2/sum(PCA_all$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all$sdev[2]^2/sum(PCA_all$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_ggplot


#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_all_category_myst <- pcscores_all_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_all_category_odont <- pcscores_all_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_all_category_ggplot <- ggplot(pcscores_all_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_all_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_all_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  theme_bw()+
  ggtitle("Whole skull")+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all$sdev[1]^2/sum(PCA_all$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all$sdev[2]^2/sum(PCA_all$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_category_ggplot

#Add phylopics for groups
PCA_all_category_ggplot <- 
  PCA_all_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.08, y = -0.15, ysize = 0.04, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.13, y = 0.2, ysize = 0.03, fill = "gray50")
PCA_all_category_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_df)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all_df)

#View results and p-value
summary(reg_PC1all_size)
summary(reg_PC2all_size)
anova(reg_PC1all_size)
anova(reg_PC2all_size)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2all_size_lm.txt")
print("PC1")
summary(reg_PC1all_size)
anova(reg_PC1all_size)
print("PC2")
summary(reg_PC2all_size)
anova(reg_PC2all_size)
sink() 

#Calculate regression for each component taking family into account
reg_PC1all_family <- lm(Comp1 ~ family, data = pcscores_all_df)
reg_PC2all_family <- lm(Comp2 ~ family, data = pcscores_all_df)

#View results and p-value
summary(reg_PC1all_family)
summary(reg_PC2all_family)
anova(reg_PC1all_family)
anova(reg_PC2all_family)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2all_family_lm.txt")
print("PC1")
summary(reg_PC1all_family)
anova(reg_PC1all_family)
print("PC2")
summary(reg_PC2all_family)
anova(reg_PC2all_family)
sink() 

#Calculate regression for each component taking category into account
reg_PC1all_category <- lm(Comp1 ~ category, data = pcscores_all_df)
reg_PC2all_category <- lm(Comp2 ~ category, data = pcscores_all_df)

#View results and p-value
summary(reg_PC1all_category)
summary(reg_PC2all_category)
anova(reg_PC1all_category)
anova(reg_PC2all_category)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2all_category_lm.txt")
print("PC1")
summary(reg_PC1all_category)
anova(reg_PC1all_category)
print("PC2")
summary(reg_PC2all_category)
anova(reg_PC2all_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1all_group <- lm(Comp1 ~ group, data = pcscores_all_df)
reg_PC2all_group <- lm(Comp2 ~ group, data = pcscores_all_df)

#View results and p-value
summary(reg_PC1all_group)
summary(reg_PC2all_group)
anova(reg_PC1all_group)
anova(reg_PC2all_group)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2all_group_lm.txt")
print("PC1")
summary(reg_PC1all_group)
anova(reg_PC1all_group)
print("PC2")
summary(reg_PC2all_group)
anova(reg_PC2all_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1all_group_cat <- lm(Comp1 ~ group * category, data = pcscores_all_df)
reg_PC2all_group_cat <- lm(Comp2 ~ group * category, data = pcscores_all_df)

#View results and p-value
summary(reg_PC1all_group_cat)
summary(reg_PC2all_group_cat)
anova(reg_PC1all_group_cat)
anova(reg_PC2all_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2all_group_cat_lm.txt")
print("PC1")
summary(reg_PC1all_group_cat)
anova(reg_PC1all_group_cat)
print("PC2")
summary(reg_PC2all_group_cat)
anova(reg_PC2all_group_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/5-PCA/PC1-2_all_lm.txt")
print("PC1")
anova(reg_PC1all_size)
anova(reg_PC1all_family)
anova(reg_PC1all_category)
anova(reg_PC1all_group)
anova(reg_PC1all_group_cat)
print("PC2")
anova(reg_PC2all_size)
anova(reg_PC2all_family)
anova(reg_PC2all_category)
anova(reg_PC2all_group)
anova(reg_PC2all_group_cat)
sink()




#Import trees in Nexus format - branch lengths needed!!
tree_specs <- "Data/tree_all_specs.txt"   #tree with selected families with category data

##Read the trees for analysis
tree <- read.nexus(tree_specs) #tree with selected families with category data
plot(tree)


#Make sure tip labels match taxa names in data frame
wrong_tips_all <- sort(tree$tip.label)
tree1 <- tree

tree1$tip.label <- levels(as.factor(gdf$Id))[match(tree1$tip.label, wrong_tips_all)]
plot(tree1)
tree <- tree1

#Check names of each tree in object
summary(tree)

coords_phylo <- gdf$coords
dimnames(coords_phylo)[[3]] <- gdf$Id

phylo.tPCA <- gm.prcomp(coords_phylo, phy = tree, 
                        GLS = TRUE, transform = TRUE)
summary(phylo.tPCA)
plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA")
plot(phylo.tPCA, phylo = F, main = "phylo PCA")
text(x = phylo.tPCA$x[,1], y = phylo.tPCA$x[,2], labels = rownames(phylo.tPCA$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)




##Rostrum ----
pca_rostrum <- gm.prcomp(gdf_rostrum$coords)

#List of PC components and proportion of variations
pca_rostrum 

#Save PCA results to file
sink("Output/5-PCA/pca_rostrum results.txt")
print("PCA rostrum")
print(pca_rostrum)
sink() 

#Change row names to codes to make plot readable
row.names(pca_rostrum$x) <- gdf_rostrum$Id

##View plot
plot(pca_rostrum, main = "PCA rostrum - PC1-PC2",  pch = 21, #title and type of point to be used
     col = mypalette_paired[5],    bg = mypalette_paired[5],  cex = 1, font.main = 2)      #improve graphics
#Add quick labels to plot
text(x = pca_rostrum$x[,1], y = pca_rostrum$x[,2], labels = rownames(pca_rostrum$x), 
     pos = 1,   offset = 0.5,  cex = 0.75)    #improve graphics 

#Save PC scores as object to use later
pcscores_rostrum <- pca_rostrum$x

#Save shapes of extremes for axes used in plot
PC1min_rostrum <- pca_rostrum[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_rostrum <- pca_rostrum[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_rostrum <- pca_rostrum[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_rostrum <- pca_rostrum[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_rostrum_points <- spheres3d(PC1min_rostrum, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA/PC1min_rostrum.png") 
rgl.snapshot(filename = "Output/5-PCA/PC1min_rostrum1.png") 
clear3d()

#PC1max colors
PC1max_rostrum_points <- spheres3d(PC1max_rostrum, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA/PC1max_rostrum1.png") 
rgl.snapshot(filename = "Output/5-PCA/PC1max_rostrum.png") 
clear3d()

#PC2min colors
PC2min_rostrum_points <- spheres3d(PC2min_rostrum, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA/PC2min_rostrum.png") 
rgl.snapshot(filename = "Output/5-PCA/PC2min_rostrum1.png") 
clear3d()

#PC2max colors
PC2max_rostrum_points <- spheres3d(PC2max_rostrum, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA/PC2max_rostrum1.png")
rgl.snapshot(filename = "Output/5-PCA/PC2max_rostrum.png") 
clear3d()

###Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_rostrum_df <- as_tibble(pcscores_rostrum)
#Add labels and other attributes to tibble as columns
pcscores_rostrum_df <- pcscores_rostrum_df %>% 
  mutate(specimens = gdf$Id,  family = gdf$family,  group = gdf$group, category = gdf$category, size = gdf$size,  
         genus = gdf$genus, module = rep("rostrum", length(gdf$Id)))
glimpse(pcscores_rostrum_df)

#Nice PCA plot with stages and groups
PCA_rostrum_ggplot <- ggplot(pcscores_rostrum_df, aes(x = Comp1, y = Comp2, label = specimens, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 60)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  theme_bw()+ 
  xlab(paste0("PC 1 (",round(as.numeric(pca_rostrum$sdev[1]^2/sum(pca_rostrum$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(pca_rostrum$sdev[2]^2/sum(pca_rostrum$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_rostrum_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_rostrum_category_myst <- pcscores_rostrum_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_rostrum_category_odont <- pcscores_rostrum_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_rostrum_category_ggplot <- ggplot(pcscores_rostrum_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_rostrum_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_rostrum_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  theme_bw()+
  ggtitle("Rostrum")+
  xlab(paste0("PC 1 (",round(as.numeric(pca_rostrum$sdev[1]^2/sum(pca_rostrum$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(pca_rostrum$sdev[2]^2/sum(pca_rostrum$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_rostrum_category_ggplot

#Add phylopics for groups
PCA_rostrum_category_ggplot <- 
  PCA_rostrum_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.12, y = -0.15, ysize = 0.04, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.2, y = 0.15, ysize = 0.03, fill = "gray50")
PCA_rostrum_category_ggplot

###Regression PC1 and PC2 rostrum ----

#Calculate regression for each component for size
reg_PC1rostrum_size <- lm(Comp1 ~ size, data = pcscores_rostrum_df)
reg_PC2rostrum_size <- lm(Comp2 ~ size, data = pcscores_rostrum_df)

#View results and p-value
summary(reg_PC1rostrum_size)
summary(reg_PC2rostrum_size)
anova(reg_PC1rostrum_size)
anova(reg_PC2rostrum_size)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2rostrum_size_lm.txt")
print("PC1")
summary(reg_PC1rostrum_size)
anova(reg_PC1rostrum_size)
print("PC2")
summary(reg_PC2rostrum_size)
anova(reg_PC2rostrum_size)
sink() 

#Calculate regression for each component taking family into account
reg_PC1rostrum_family <- lm(Comp1 ~ family, data = pcscores_rostrum_df)
reg_PC2rostrum_family <- lm(Comp2 ~ family, data = pcscores_rostrum_df)

#View results and p-value
summary(reg_PC1rostrum_family)
summary(reg_PC2rostrum_family)
anova(reg_PC1rostrum_family)
anova(reg_PC2rostrum_family)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2rostrum_family_lm.txt")
print("PC1")
summary(reg_PC1rostrum_family)
anova(reg_PC1rostrum_family)
print("PC2")
summary(reg_PC2rostrum_family)
anova(reg_PC2rostrum_family)
sink() 

#Calculate regression for each component taking category into account
reg_PC1rostrum_category <- lm(Comp1 ~ category, data = pcscores_rostrum_df)
reg_PC2rostrum_category <- lm(Comp2 ~ category, data = pcscores_rostrum_df)

#View results and p-value
summary(reg_PC1rostrum_category)
summary(reg_PC2rostrum_category)
anova(reg_PC1rostrum_category)
anova(reg_PC2rostrum_category)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2rostrum_category_lm.txt")
print("PC1")
summary(reg_PC1rostrum_category)
anova(reg_PC1rostrum_category)
print("PC2")
summary(reg_PC2rostrum_category)
anova(reg_PC2rostrum_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1rostrum_group <- lm(Comp1 ~ group, data = pcscores_rostrum_df)
reg_PC2rostrum_group <- lm(Comp2 ~ group, data = pcscores_rostrum_df)

#View results and p-value
summary(reg_PC1rostrum_group)
summary(reg_PC2rostrum_group)
anova(reg_PC1rostrum_group)
anova(reg_PC2rostrum_group)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2rostrum_group_lm.txt")
print("PC1")
summary(reg_PC1rostrum_group)
anova(reg_PC1rostrum_group)
print("PC2")
summary(reg_PC2rostrum_group)
anova(reg_PC2rostrum_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1rostrum_group_cat <- lm(Comp1 ~ group * category, data = pcscores_rostrum_df)
reg_PC2rostrum_group_cat <- lm(Comp2 ~ group * category, data = pcscores_rostrum_df)

#View results and p-value
summary(reg_PC1rostrum_group_cat)
summary(reg_PC2rostrum_group_cat)
anova(reg_PC1rostrum_group_cat)
anova(reg_PC2rostrum_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2rostrum_group_cat_lm.txt")
print("PC1")
summary(reg_PC1rostrum_group_cat)
anova(reg_PC1rostrum_group_cat)
print("PC2")
summary(reg_PC2rostrum_group_cat)
anova(reg_PC2rostrum_group_cat)
sink() 

#Save results of rostrum regressions to 1 file
sink("Output/5-PCA/PC1-2_rostrum_lm.txt")
print("PC1")
anova(reg_PC1rostrum_size)
anova(reg_PC1rostrum_family)
anova(reg_PC1rostrum_category)
anova(reg_PC1rostrum_group)
anova(reg_PC1rostrum_group_cat)
print("PC2")
anova(reg_PC2rostrum_size)
anova(reg_PC2rostrum_family)
anova(reg_PC2rostrum_category)
anova(reg_PC2rostrum_group)
anova(reg_PC2rostrum_group_cat)
sink()


##Braincase ----
pca_braincase <- gm.prcomp(gdf_braincase$coords)

#List of PC components and proportion of variations
pca_braincase 

#Save PCA results to file
sink("Output/5-PCA/pca_braincase results.txt")
print("PCA braincase")
print(pca_braincase)
sink() 

#Change row names to codes to make plot readable
row.names(pca_braincase$x) <- gdf_braincase$Id

##View plot
plot(pca_braincase, main = "PCA braincase - PC1-PC2",  pch = 21, #title and type of point to be used
     col = mypalette_paired[1],    bg = mypalette_paired[1],  cex = 1, font.main = 2)      #improve graphics
#Add quick labels to plot
text(x = pca_braincase$x[,1], y = pca_braincase$x[,2], labels = rownames(pca_braincase$x), 
     pos = 1,   offset = 0.5,  cex = 0.75)    #improve graphics 

#Save PC scores as object to use later
pcscores_braincase <- pca_braincase$x

#Save shapes of extremes for axes used in plot
PC1min_braincase <- pca_braincase[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_braincase <- pca_braincase[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_braincase <- pca_braincase[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_braincase <- pca_braincase[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_braincase_points <- spheres3d(PC1min_braincase, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA/PC1min_braincase.png") 
rgl.snapshot(filename = "Output/5-PCA/PC1min_braincase1.png") 
clear3d()

#PC1max colors
PC1max_braincase_points <- spheres3d(PC1max_braincase, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA/PC1max_braincase1.png") 
rgl.snapshot(filename = "Output/5-PCA/PC1max_braincase.png") 
clear3d()

#PC2min colors
PC2min_braincase_points <- spheres3d(PC2min_braincase, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA/PC2min_braincase.png") 
rgl.snapshot(filename = "Output/5-PCA/PC2min_braincase1.png") 
clear3d()

#PC2max colors
PC2max_braincase_points <- spheres3d(PC2max_braincase, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA/PC2max_braincase1.png")
rgl.snapshot(filename = "Output/5-PCA/PC2max_braincase.png") 
clear3d()

###Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_braincase_df <- as_tibble(pcscores_braincase)
#Add labels and other attributes to tibble as columns
pcscores_braincase_df <- pcscores_braincase_df %>% 
  mutate(specimens = gdf$Id,  family = gdf$family,  group = gdf$group, category = gdf$category, size = gdf$size, 
         genus = gdf$genus, module = rep("braincase", length(gdf$Id)))
glimpse(pcscores_braincase_df)

#Nice PCA plot with stages and groups
PCA_braincase_ggplot <- ggplot(pcscores_braincase_df, aes(x = Comp1, y = Comp2, label = specimens, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 60)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  scale_y_reverse()+ #to match other modules
  theme_bw()+ 
  xlab(paste0("PC 1 (",round(as.numeric(pca_braincase$sdev[1]^2/sum(pca_braincase$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(pca_braincase$sdev[2]^2/sum(pca_braincase$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_braincase_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_braincase_category_myst <- pcscores_braincase_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_braincase_category_odont <- pcscores_braincase_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_braincase_category_ggplot <- ggplot(pcscores_braincase_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_braincase_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_braincase_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_y_reverse()+
  theme_bw()+
  ggtitle("Braincase")+
  xlab(paste0("PC 1 (",round(as.numeric(pca_braincase$sdev[1]^2/sum(pca_braincase$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(pca_braincase$sdev[2]^2/sum(pca_braincase$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_braincase_category_ggplot

#Add phylopics for groups
PCA_braincase_category_ggplot <- 
  PCA_braincase_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.22, y = 0.22, ysize = 0.025, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.22, y = -0.12, ysize = 0.02, fill = "gray50")
PCA_braincase_category_ggplot

ggarrange(PCA_rostrum_category_ggplot, PCA_braincase_category_ggplot, ncol =2 , nrow =1, common.legend = T, legend = "bottom")

#Arrange all PCA plots
ggarrange(PCA_all_category_ggplot,
          PCA_rostrum_category_ggplot, PCA_braincase_category_ggplot, ncol =3 , nrow =1, common.legend = T, legend = "bottom")

###Regression PC1 and PC2 braincase ----

#Calculate regression for each component for size
reg_PC1braincase_size <- lm(Comp1 ~ size, data = pcscores_braincase_df)
reg_PC2braincase_size <- lm(Comp2 ~ size, data = pcscores_braincase_df)

#View results and p-value
summary(reg_PC1braincase_size)
summary(reg_PC2braincase_size)
anova(reg_PC1braincase_size)
anova(reg_PC2braincase_size)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2braincase_size_lm.txt")
print("PC1")
summary(reg_PC1braincase_size)
anova(reg_PC1braincase_size)
print("PC2")
summary(reg_PC2braincase_size)
anova(reg_PC2braincase_size)
sink() 

#Calculate regression for each component taking family into account
reg_PC1braincase_family <- lm(Comp1 ~ family, data = pcscores_braincase_df)
reg_PC2braincase_family <- lm(Comp2 ~ family, data = pcscores_braincase_df)

#View results and p-value
summary(reg_PC1braincase_family)
summary(reg_PC2braincase_family)
anova(reg_PC1braincase_family)
anova(reg_PC2braincase_family)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2braincase_family_lm.txt")
print("PC1")
summary(reg_PC1braincase_family)
anova(reg_PC1braincase_family)
print("PC2")
summary(reg_PC2braincase_family)
anova(reg_PC2braincase_family)
sink() 

#Calculate regression for each component taking category into account
reg_PC1braincase_category <- lm(Comp1 ~ category, data = pcscores_braincase_df)
reg_PC2braincase_category <- lm(Comp2 ~ category, data = pcscores_braincase_df)

#View results and p-value
summary(reg_PC1braincase_category)
summary(reg_PC2braincase_category)
anova(reg_PC1braincase_category)
anova(reg_PC2braincase_category)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2braincase_category_lm.txt")
print("PC1")
summary(reg_PC1braincase_category)
anova(reg_PC1braincase_category)
print("PC2")
summary(reg_PC2braincase_category)
anova(reg_PC2braincase_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1braincase_group <- lm(Comp1 ~ group, data = pcscores_braincase_df)
reg_PC2braincase_group <- lm(Comp2 ~ group, data = pcscores_braincase_df)

#View results and p-value
summary(reg_PC1braincase_group)
summary(reg_PC2braincase_group)
anova(reg_PC1braincase_group)
anova(reg_PC2braincase_group)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2braincase_group_lm.txt")
print("PC1")
summary(reg_PC1braincase_group)
anova(reg_PC1braincase_group)
print("PC2")
summary(reg_PC2braincase_group)
anova(reg_PC2braincase_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1braincase_group_cat <- lm(Comp1 ~ group * category, data = pcscores_braincase_df)
reg_PC2braincase_group_cat <- lm(Comp2 ~ group * category, data = pcscores_braincase_df)

#View results and p-value
summary(reg_PC1braincase_group_cat)
summary(reg_PC2braincase_group_cat)
anova(reg_PC1braincase_group_cat)
anova(reg_PC2braincase_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA/PC1-2braincase_group_cat_lm.txt")
print("PC1")
summary(reg_PC1braincase_group_cat)
anova(reg_PC1braincase_group_cat)
print("PC2")
summary(reg_PC2braincase_group_cat)
anova(reg_PC2braincase_group_cat)
sink() 

#Save results of braincase regressions to 1 file
sink("Output/5-PCA/PC1-2_braincase_lm.txt")
print("PC1")
anova(reg_PC1braincase_size)
anova(reg_PC1braincase_family)
anova(reg_PC1braincase_category)
anova(reg_PC1braincase_group)
anova(reg_PC1braincase_group_cat)
print("PC2")
anova(reg_PC2braincase_size)
anova(reg_PC2braincase_family)
anova(reg_PC2braincase_category)
anova(reg_PC2braincase_group)
anova(reg_PC2braincase_group_cat)
sink()


###### 

#Next - ch. 6 - Disparity analyses