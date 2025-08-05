#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.5 - PCA whole skull, rostrum and braincase mean shapes

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
library(abind)

#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#Create GDF rostrum and braincase ----
#Extract the rostrum and braincase from common alignment

#Create rostrum and briancase partitions based on developmental hypothesis in Goswami et al. 2022
#Rostrum -> Neural crest
#Neural crest - maxilla, premaxilla, palatine, vomer, nasal, squamosal, orbit
#Briancase  -> Mesoderm
#Mesoderm - supraoccipital, exoccipital, interparietal, condyle, basioccipital
landmarks <- 1:dim(gdf$coords)[[1]]

rostrum <- landmarks[which(modules_all %in% c("maxilla", "premaxilla", "palatine", "nasal", "squamosal", "orbit"))]  

braincase <- landmarks[-rostrum]

#Plot on surface to check assigment
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[rostrum,,41], col =  mypalette_paired[5], type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(shape_array[braincase,,41], col =  mypalette_paired[1], type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)


##Make data frame for analyses in geomorph separate for rostrum and braincase
gdf_mean_rostrum <- geomorph.data.frame(coords = gpa_means$coords[rostrum,,],
                                    genus = gdf_mean_all$genus, category = gdf_mean_all$category,
                                    family =  gdf_mean_all$family, group = gdf_mean_all$group,
                                    size = gdf_mean_all$size) #overall skull size
glimpse(gdf_mean_rostrum)

gdf_mean_braincase <- geomorph.data.frame(coords = gpa_means$coords[braincase,,],
                                        genus = gdf_mean_all$genus, category = gdf_mean_all$category,
                                        family =  gdf_mean_all$family, group = gdf_mean_all$group,
                                        size = gdf_mean_all$size) #overall skull size
glimpse(gdf_mean_braincase)


#PCA COMPLETE DATASET ----

##Whole skull ----
#Run PCA on complete dataset
PCA_all_means <- gm.prcomp(gdf_mean_all$coords)

#List of PC components and proportion of variation
PCA_all_means

#Save PCA results to file
sink("Output/5-PCA means/PCA_all_means_components.txt")
print("PCA complete dataset")
print(PCA_all_means)
sink() 

##View plot
plot(PCA_all_means, main = "PCA all data mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all_means$x[,1], y = PCA_all_means$x[,2], labels = rownames(PCA_all_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

##View plot
plot(PCA_all_means, axis1 = 1, axis2 =3, main = "PCA all data mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "orange",   #border of points
     bg = "orange",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all_means$x[,1], y = PCA_all_means$x[,3], labels = rownames(PCA_all_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all_means <- PCA_all_means$x 

#Save shapes of extremes for axes used in plot
PC1min_all_means <- PCA_all_means[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all_means <- PCA_all_means[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all_means <- PCA_all_means[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all_means <- PCA_all_means[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_all_means_points <- spheres3d(PC1min_all_means, radius=.0008, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_all_means.png")
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_all_means1.png") 
clear3d()

#PC1max colors
PC1max_all_means_points <- spheres3d(PC1max_all_means, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_all_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_all_means.png") 
clear3d()

#PC2min colors
PC2min_all_means_points <- spheres3d(PC2min_all_means, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_all_means.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_all_means1.png") 
clear3d()

#PC2max colors
PC2max_all_means_points <- spheres3d(PC2max_all_means, radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_all_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_all_means.png")  
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all_means_df <- as_tibble(pcscores_all_means)
#Add labels and other attributes to tibble as columns
pcscores_all_means_df <- pcscores_all_means_df %>% mutate(group = gdf_mean_all$group, category = gdf_mean_all$category,
                                              genus = gdf_mean_all$genus, size = gdf_mean_all$size, family = gdf_mean_all$family,
                                              module = rep("skull", length(gdf_mean_all$genus)))
glimpse(pcscores_all_means_df)

#Nice PCA plot with stages and groups
PCA_all_means_ggplot <- ggplot(pcscores_all_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_means$sdev[1]^2/sum(PCA_all_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_means$sdev[2]^2/sum(PCA_all_means$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_all_means_category_myst <- pcscores_all_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_all_means_category_odont <- pcscores_all_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_all_means_category_ggplot <- ggplot(pcscores_all_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_all_means_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_all_means_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  scale_y_reverse()+
  theme_bw()+
  ggtitle("Whole skull")+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_means$sdev[1]^2/sum(PCA_all_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_means$sdev[2]^2/sum(PCA_all_means$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_means_category_ggplot

#Add phylopics for groups
PCA_all_means_category_ggplot <- 
  PCA_all_means_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.25, y = -0.06, height = 0.03, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.11, y = -0.22, height = 0.02, fill = "gray50")
PCA_all_means_category_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_means_size <- lm(Comp1 ~ size, data = pcscores_all_means_df)
reg_PC2all_means_size <- lm(Comp2 ~ size, data = pcscores_all_means_df)

#View results and p-value
summary(reg_PC1all_means_size)
summary(reg_PC2all_means_size)
anova(reg_PC1all_means_size)
anova(reg_PC2all_means_size)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_size_lm.txt")
print("PC1")
summary(reg_PC1all_means_size)
anova(reg_PC1all_means_size)
print("PC2")
summary(reg_PC2all_means_size)
anova(reg_PC2all_means_size)
sink() 

#Calculate regression for each component for size and category
reg_PC1all_means_size_cat <- lm(Comp1 ~ size * category, data = pcscores_all_means_df)
reg_PC2all_means_size_cat <- lm(Comp2 ~ size * category, data = pcscores_all_means_df)

#View results and p-value
summary(reg_PC1all_means_size_cat)
summary(reg_PC2all_means_size_cat)
anova(reg_PC1all_means_size_cat)
anova(reg_PC2all_means_size_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_size_cat_lm.txt")
print("PC1")
summary(reg_PC1all_means_size_cat)
anova(reg_PC1all_means_size_cat)
print("PC2")
summary(reg_PC2all_means_size_cat)
anova(reg_PC2all_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_PC1all_means_category <- lm(Comp1 ~ category, data = pcscores_all_means_df)
reg_PC2all_means_category <- lm(Comp2 ~ category, data = pcscores_all_means_df)

#View results and p-value
summary(reg_PC1all_means_category)
summary(reg_PC2all_means_category)
anova(reg_PC1all_means_category)
anova(reg_PC2all_means_category)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_category_lm.txt")
print("PC1")
summary(reg_PC1all_means_category)
anova(reg_PC1all_means_category)
print("PC2")
summary(reg_PC2all_means_category)
anova(reg_PC2all_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1all_means_group <- lm(Comp1 ~ group, data = pcscores_all_means_df)
reg_PC2all_means_group <- lm(Comp2 ~ group, data = pcscores_all_means_df)

#View results and p-value
summary(reg_PC1all_means_group)
summary(reg_PC2all_means_group)
anova(reg_PC1all_means_group)
anova(reg_PC2all_means_group)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_group_lm.txt")
print("PC1")
summary(reg_PC1all_means_group)
anova(reg_PC1all_means_group)
print("PC2")
summary(reg_PC2all_means_group)
anova(reg_PC2all_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1all_means_group_cat <- lm(Comp1 ~ group * category, data = pcscores_all_means_df)
reg_PC2all_means_group_cat <- lm(Comp2 ~ group * category, data = pcscores_all_means_df)

#View results and p-value
summary(reg_PC1all_means_group_cat)
summary(reg_PC2all_means_group_cat)
anova(reg_PC1all_means_group_cat)
anova(reg_PC2all_means_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_group_cat_lm.txt")
print("PC1")
summary(reg_PC1all_means_group_cat)
anova(reg_PC1all_means_group_cat)
print("PC2")
summary(reg_PC2all_means_group_cat)
anova(reg_PC2all_means_group_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/5-PCA means/PC1-2_all_means_lm.txt")
print("PC1")
anova(reg_PC1all_means_size)
anova(reg_PC1all_means_size_cat)
anova(reg_PC1all_means_category)
anova(reg_PC1all_means_group)
anova(reg_PC1all_means_group_cat)
print("PC2")
anova(reg_PC2all_means_size)
anova(reg_PC2all_means_size_cat)
anova(reg_PC2all_means_category)
anova(reg_PC2all_means_group)
anova(reg_PC2all_means_group_cat)
sink()


##Rostrum ----
#Run PCA on rostrum dataset
PCA_rostrum_means <- gm.prcomp(gdf_mean_rostrum$coords)

#List of PC components and proportion of variation
PCA_rostrum_means

#Save PCA results to file
sink("Output/5-PCA means/PCA_rostrum_means_components.txt")
print("PCA complete dataset")
print(PCA_rostrum_means)
sink() 

##View plot
plot(PCA_rostrum_means, main = "PCA rostrum mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_rostrum_means$x[,1], y = PCA_rostrum_means$x[,2], labels = rownames(PCA_rostrum_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

##View plot
plot(PCA_rostrum_means, axis1 = 1, axis2 =3, main = "PCA rostrum mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "orange",   #border of points
     bg = "orange",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_rostrum_means$x[,1], y = PCA_rostrum_means$x[,3], labels = rownames(PCA_rostrum_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_rostrum_means <- PCA_rostrum_means$x 

#Save shapes of extremes for axes used in plot
PC1min_rostrum_means <- PCA_rostrum_means[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_rostrum_means <- PCA_rostrum_means[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_rostrum_means <- PCA_rostrum_means[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_rostrum_means <- PCA_rostrum_means[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_rostrum_means_points <- spheres3d(PC1min_rostrum_means, radius=.0008, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_rostrum_means.png")
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_rostrum_means1.png") 
clear3d()

#PC1max colors
PC1max_rostrum_means_points <- spheres3d(PC1max_rostrum_means, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_rostrum_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_rostrum_means.png") 
clear3d()

#PC2min colors
PC2min_rostrum_means_points <- spheres3d(PC2min_rostrum_means, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_rostrum_means.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_rostrum_means1.png") 
clear3d()

#PC2max colors
PC2max_rostrum_means_points <- spheres3d(PC2max_rostrum_means, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_rostrum_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_rostrum_means.png")  
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_rostrum_means_df <- as_tibble(pcscores_rostrum_means)
#Add labels and other attributes to tibble as columns
pcscores_rostrum_means_df <- pcscores_rostrum_means_df %>% mutate(group = gdf_mean_rostrum$group, category = gdf_mean_rostrum$category,
                                                          genus = gdf_mean_rostrum$genus, size = gdf_mean_rostrum$size, family = gdf_mean_rostrum$family,
                                                          module = rep("rostrum", length(gdf_mean_rostrum$genus)))
glimpse(pcscores_rostrum_means_df)

#Nice PCA plot with stages and groups
PCA_rostrum_means_ggplot <- ggplot(pcscores_rostrum_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 50)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_rostrum_means$sdev[1]^2/sum(PCA_rostrum_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_rostrum_means$sdev[2]^2/sum(PCA_rostrum_means$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_rostrum_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_rostrum_means_category_myst <- pcscores_rostrum_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_rostrum_means_category_odont <- pcscores_rostrum_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_rostrum_means_category_ggplot <- ggplot(pcscores_rostrum_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_rostrum_means_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_rostrum_means_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  theme_bw()+
  ggtitle("Rostrum")+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_rostrum_means$sdev[1]^2/sum(PCA_rostrum_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_rostrum_means$sdev[2]^2/sum(PCA_rostrum_means$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_rostrum_means_category_ggplot

#Add phylopics for groups
PCA_rostrum_means_category_ggplot <- 
  PCA_rostrum_means_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.25, y = 0.17, height = 0.03, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.08, y = -0.22, height = 0.02, fill = "gray50")
PCA_rostrum_means_category_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_means_size <- lm(Comp1 ~ size, data = pcscores_rostrum_means_df)
reg_PC2all_means_size <- lm(Comp2 ~ size, data = pcscores_rostrum_means_df)

#View results and p-value
summary(reg_PC1all_means_size)
summary(reg_PC2all_means_size)
anova(reg_PC1all_means_size)
anova(reg_PC2all_means_size)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_size_lm.txt")
print("PC1")
summary(reg_PC1all_means_size)
anova(reg_PC1all_means_size)
print("PC2")
summary(reg_PC2all_means_size)
anova(reg_PC2all_means_size)
sink() 

#Calculate regression for each component for size and category
reg_PC1all_means_size_cat <- lm(Comp1 ~ size * category, data = pcscores_rostrum_means_df)
reg_PC2all_means_size_cat <- lm(Comp2 ~ size * category, data = pcscores_rostrum_means_df)

#View results and p-value
summary(reg_PC1all_means_size_cat)
summary(reg_PC2all_means_size_cat)
anova(reg_PC1all_means_size_cat)
anova(reg_PC2all_means_size_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_size_cat_lm.txt")
print("PC1")
summary(reg_PC1all_means_size_cat)
anova(reg_PC1all_means_size_cat)
print("PC2")
summary(reg_PC2all_means_size_cat)
anova(reg_PC2all_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_PC1all_means_category <- lm(Comp1 ~ category, data = pcscores_rostrum_means_df)
reg_PC2all_means_category <- lm(Comp2 ~ category, data = pcscores_rostrum_means_df)

#View results and p-value
summary(reg_PC1all_means_category)
summary(reg_PC2all_means_category)
anova(reg_PC1all_means_category)
anova(reg_PC2all_means_category)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_category_lm.txt")
print("PC1")
summary(reg_PC1all_means_category)
anova(reg_PC1all_means_category)
print("PC2")
summary(reg_PC2all_means_category)
anova(reg_PC2all_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1all_means_group <- lm(Comp1 ~ group, data = pcscores_rostrum_means_df)
reg_PC2all_means_group <- lm(Comp2 ~ group, data = pcscores_rostrum_means_df)

#View results and p-value
summary(reg_PC1all_means_group)
summary(reg_PC2all_means_group)
anova(reg_PC1all_means_group)
anova(reg_PC2all_means_group)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_group_lm.txt")
print("PC1")
summary(reg_PC1all_means_group)
anova(reg_PC1all_means_group)
print("PC2")
summary(reg_PC2all_means_group)
anova(reg_PC2all_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1all_means_group_cat <- lm(Comp1 ~ group * category, data = pcscores_rostrum_means_df)
reg_PC2all_means_group_cat <- lm(Comp2 ~ group * category, data = pcscores_rostrum_means_df)

#View results and p-value
summary(reg_PC1all_means_group_cat)
summary(reg_PC2all_means_group_cat)
anova(reg_PC1all_means_group_cat)
anova(reg_PC2all_means_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2all_means_group_cat_lm.txt")
print("PC1")
summary(reg_PC1all_means_group_cat)
anova(reg_PC1all_means_group_cat)
print("PC2")
summary(reg_PC2all_means_group_cat)
anova(reg_PC2all_means_group_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/5-PCA means/PC1-2_rostrum_means_lm.txt")
print("PC1")
anova(reg_PC1all_means_size)
anova(reg_PC1all_means_size_cat)
anova(reg_PC1all_means_category)
anova(reg_PC1all_means_group)
anova(reg_PC1all_means_group_cat)
print("PC2")
anova(reg_PC2all_means_size)
anova(reg_PC2all_means_size_cat)
anova(reg_PC2all_means_category)
anova(reg_PC2all_means_group)
anova(reg_PC2all_means_group_cat)
sink()

##Braincase ----
#Run PCA on braincase dataset
PCA_braincase_means <- gm.prcomp(gdf_mean_braincase$coords)

#List of PC components and proportion of variation
PCA_braincase_means

#Save PCA results to file
sink("Output/5-PCA means/PCA_braincase_means_components.txt")
print("PCA complete dataset")
print(PCA_braincase_means)
sink() 

##View plot
plot(PCA_braincase_means, main = "PCA braincase mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_braincase_means$x[,1], y = PCA_braincase_means$x[,2], labels = rownames(PCA_braincase_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

##View plot
plot(PCA_braincase_means, axis1 = 1, axis2 =3, main = "PCA braincase mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "orange",   #border of points
     bg = "orange",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_braincase_means$x[,1], y = PCA_braincase_means$x[,3], labels = rownames(PCA_braincase_means$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_braincase_means <- PCA_braincase_means$x 

#Save shapes of extremes for axes used in plot
PC1min_braincase_means <- PCA_braincase_means[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_braincase_means <- PCA_braincase_means[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_braincase_means <- PCA_braincase_means[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_braincase_means <- PCA_braincase_means[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_braincase_means_points <- spheres3d(PC1min_braincase_means, radius=.0008, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_braincase_means.png")
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1min_braincase_means1.png") 
clear3d()

#PC1max colors
PC1max_braincase_means_points <- spheres3d(PC1max_braincase_means, radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_braincase_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC1max_braincase_means.png") 
clear3d()

#PC2min colors
PC2min_braincase_means_points <- spheres3d(PC2min_braincase_means, radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_braincase_means.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2min_braincase_means1.png") 
clear3d()

#PC2max colors
PC2max_braincase_means_points <- spheres3d(PC2max_braincase_means, radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_braincase_means1.png") 
rgl.snapshot(filename = "Output/5-PCA means/min-max shapes/PC2max_braincase_means.png")  
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_braincase_means_df <- as_tibble(pcscores_braincase_means)
#Add labels and other attributes to tibble as columns
pcscores_braincase_means_df <- pcscores_braincase_means_df %>% mutate(group = gdf_mean_all$group, category = gdf_mean_all$category,
                                                          genus = gdf_mean_all$genus, size = gdf_mean_all$size, family = gdf_mean_all$family,
                                                          module = rep("braincase", length(gdf_mean_all$genus)))
glimpse(pcscores_braincase_means_df)

#Nice PCA plot with stages and groups
PCA_braincase_means_ggplot <- ggplot(pcscores_braincase_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = family, fill = family))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 50)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Families", labels = levels(families), #copy from as.factor(genera)
                      values = mypalette_families, aesthetics = c("colour","fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_braincase_means$sdev[1]^2/sum(PCA_braincase_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_braincase_means$sdev[2]^2/sum(PCA_braincase_means$sdev^2)*100), digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_braincase_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_braincase_means_category_myst <- pcscores_braincase_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_braincase_means_category_odont <- pcscores_braincase_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_braincase_means_category_ggplot <- ggplot(pcscores_braincase_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_braincase_means_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_braincase_means_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, linewidth = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  theme_bw()+
  ggtitle("Braincase")+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_braincase_means$sdev[1]^2/sum(PCA_braincase_means$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_braincase_means$sdev[2]^2/sum(PCA_braincase_means$sdev^2)*100), digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_braincase_means_category_ggplot

#Add phylopics for groups
PCA_braincase_means_category_ggplot <- 
  PCA_braincase_means_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.05, y = -0.13, height = 0.015, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.1, y = 0.05, height = 0.01, fill = "gray50")
PCA_braincase_means_category_ggplot

#Arrange all PCA plots
ggarrange(PCA_all_means_category_ggplot,
          PCA_rostrum_means_category_ggplot, PCA_braincase_means_category_ggplot, ncol =3 , nrow =1, common.legend = T, legend = "bottom")


###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1braincase_means_size <- lm(Comp1 ~ size, data = pcscores_braincase_means_df)
reg_PC2braincase_means_size <- lm(Comp2 ~ size, data = pcscores_braincase_means_df)

#View results and p-value
summary(reg_PC1braincase_means_size)
summary(reg_PC2braincase_means_size)
anova(reg_PC1braincase_means_size)
anova(reg_PC2braincase_means_size)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2braincase_means_size_lm.txt")
print("PC1")
summary(reg_PC1braincase_means_size)
anova(reg_PC1braincase_means_size)
print("PC2")
summary(reg_PC2braincase_means_size)
anova(reg_PC2braincase_means_size)
sink() 

#Calculate regression for each component for size and category
reg_PC1braincase_means_size_cat <- lm(Comp1 ~ size * category, data = pcscores_braincase_means_df)
reg_PC2braincase_means_size_cat <- lm(Comp2 ~ size * category, data = pcscores_braincase_means_df)

#View results and p-value
summary(reg_PC1braincase_means_size_cat)
summary(reg_PC2braincase_means_size_cat)
anova(reg_PC1braincase_means_size_cat)
anova(reg_PC2braincase_means_size_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2braincase_means_size_cat_lm.txt")
print("PC1")
summary(reg_PC1braincase_means_size_cat)
anova(reg_PC1braincase_means_size_cat)
print("PC2")
summary(reg_PC2braincase_means_size_cat)
anova(reg_PC2braincase_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_PC1braincase_means_category <- lm(Comp1 ~ category, data = pcscores_braincase_means_df)
reg_PC2braincase_means_category <- lm(Comp2 ~ category, data = pcscores_braincase_means_df)

#View results and p-value
summary(reg_PC1braincase_means_category)
summary(reg_PC2braincase_means_category)
anova(reg_PC1braincase_means_category)
anova(reg_PC2braincase_means_category)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2braincase_means_category_lm.txt")
print("PC1")
summary(reg_PC1braincase_means_category)
anova(reg_PC1braincase_means_category)
print("PC2")
summary(reg_PC2braincase_means_category)
anova(reg_PC2braincase_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1braincase_means_group <- lm(Comp1 ~ group, data = pcscores_braincase_means_df)
reg_PC2braincase_means_group <- lm(Comp2 ~ group, data = pcscores_braincase_means_df)

#View results and p-value
summary(reg_PC1braincase_means_group)
summary(reg_PC2braincase_means_group)
anova(reg_PC1braincase_means_group)
anova(reg_PC2braincase_means_group)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2braincase_means_group_lm.txt")
print("PC1")
summary(reg_PC1braincase_means_group)
anova(reg_PC1braincase_means_group)
print("PC2")
summary(reg_PC2braincase_means_group)
anova(reg_PC2braincase_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1braincase_means_group_cat <- lm(Comp1 ~ group * category, data = pcscores_braincase_means_df)
reg_PC2braincase_means_group_cat <- lm(Comp2 ~ group * category, data = pcscores_braincase_means_df)

#View results and p-value
summary(reg_PC1braincase_means_group_cat)
summary(reg_PC2braincase_means_group_cat)
anova(reg_PC1braincase_means_group_cat)
anova(reg_PC2braincase_means_group_cat)

#Save results of significant regression to file
sink("Output/5-PCA means/PC1-2braincase_means_group_cat_lm.txt")
print("PC1")
summary(reg_PC1braincase_means_group_cat)
anova(reg_PC1braincase_means_group_cat)
print("PC2")
summary(reg_PC2braincase_means_group_cat)
anova(reg_PC2braincase_means_group_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/5-PCA means/PC1-2_braincase_means_lm.txt")
print("PC1")
anova(reg_PC1braincase_means_size)
anova(reg_PC1braincase_means_size_cat)
anova(reg_PC1braincase_means_category)
anova(reg_PC1braincase_means_group)
anova(reg_PC1braincase_means_group_cat)
print("PC2")
anova(reg_PC2braincase_means_size)
anova(reg_PC2braincase_means_size_cat)
anova(reg_PC2braincase_means_category)
anova(reg_PC2braincase_means_group)
anova(reg_PC2braincase_means_group_cat)
sink()

###### 

#Next - ch. 6 - Disparity analyses mean shapes