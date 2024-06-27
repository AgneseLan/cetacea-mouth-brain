#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.5b - PCA phylogenetically transformed whole skull, rostrum and braincase mean shapes

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

phyloPCA_all_means <- list()

for (c in 1:length(categories_list)){    
phyloPCA_all_means[[c]] <- gm.prcomp(gdf_mean_shapes[[c]]$coords, phy = trees_list[[c]], 
                          GLS = TRUE, transform = TRUE)
}

#List of PC components and proportion of variation
phyloPCA_all_means

#Save PCA results to file
sink("Output/5b-phyloPCA means/phyloPCA_all_means_components.txt")
print("phyloPCA complete dataset")
print(phyloPCA_all_means)
sink() 

##View plot
for (c in 1:length(categories_list)){    
plot(phyloPCA_all_means[[c]], phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA all data mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = phyloPCA_all_means[[c]]$x[,1], y = phyloPCA_all_means$x[,2], labels = rownames(phyloPCA_all_means[[c]]$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)
}

##View plot
for (c in 1:length(categories_list)){   
  plot(phyloPCA_all_means[[c]], axis1 = 1, axis2 = 3, phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA all data mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "orange",   #border of points
     bg = "orange",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = phyloPCA_all_means[[c]]$x[,1], y = phyloPCA_all_means[[c]]$x[,3], labels = rownames(phyloPCA_all_means[[c]]$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)
}

#Save PC scores as object to use later
phylo_pcscores_all_means <- list()

for (c in 1:length(categories_list)){   
  phylo_pcscores_all_means[[c]] <- phyloPCA_all_means[[c]]$x 
}

#Save shapes of extremes for axes used in plot
phyloPC1min_all_means <- list()
phyloPC1max_all_means <- list()
phyloPC2min_all_means <- list()
phyloPC2max_all_means <- list()

for (c in 1:length(categories_list)){   
phyloPC1min_all_means[[c]] <- phyloPCA_all_means[[c]][["shapes"]][["shapes.comp1"]][["min"]]
phyloPC1max_all_means[[c]] <- phyloPCA_all_means[[c]][["shapes"]][["shapes.comp1"]][["max"]] 
phyloPC2min_all_means[[c]] <- phyloPCA_all_means[[c]][["shapes"]][["shapes.comp2"]][["min"]] 
phyloPC2max_all_means[[c]] <- phyloPCA_all_means[[c]][["shapes"]][["shapes.comp2"]][["max"]] 
}

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#Early
#phyloPC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_all_means_points <- spheres3d(phyloPC1min_all_means[[1]], radius=.0008, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means_early.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means1_early.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_all_means_points <- spheres3d(phyloPC1max_all_means[[1]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means_early.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_all_means_points <- spheres3d(phyloPC2min_all_means[[1]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means1_early.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_all_means_points <- spheres3d(phyloPC2max_all_means[[1]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means_early.png")  
clear3d()

#Late/New
#phyloPC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_all_means_points <- spheres3d(phyloPC1min_all_means[[2]], radius=.0008, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means_late_new.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means1_late_new.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_all_means_points <- spheres3d(phyloPC1max_all_means[[2]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means_late_new.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_all_means_points <- spheres3d(phyloPC2min_all_means[[2]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means1_late_new.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_all_means_points <- spheres3d(phyloPC2max_all_means[[2]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means_late_new.png")  
clear3d()

#Immature
#phyloPC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_all_means_points <- spheres3d(phyloPC1min_all_means[[3]], radius=.0008, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means_immature.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means1_immature.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_all_means_points <- spheres3d(phyloPC1max_all_means[[3]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means_immature.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_all_means_points <- spheres3d(phyloPC2min_all_means[[3]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means1_immature.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_all_means_points <- spheres3d(phyloPC2max_all_means[[3]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means_immature.png")  
clear3d()

#Adult
#phyloPC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_all_means_points <- spheres3d(phyloPC1min_all_means[[4]], radius=.0008, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means_adult.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_all_means1_adult.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_all_means_points <- spheres3d(phyloPC1max_all_means[[4]], radius=.001, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_all_means_adult.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_all_means_points <- spheres3d(phyloPC2min_all_means[[4]], radius=.0007, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_all_means1_adult.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_all_means_points <- spheres3d(phyloPC2max_all_means[[4]], radius=.002, color = col_modules)
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_all_means_adult.png")  
clear3d()

##Make better PCA plot using ggplot
#Check how many components explain >95% variation and compare with total components each category
phyloPCA_all_means_prop <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_all_means_prop[[c]] <- phyloPCA_all_means[[c]]$sdev^2/ sum(phyloPCA_all_means[[c]]$sdev^2)
}

phyloPCA_all_means_cum <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_all_means_cum[[c]] <- cumsum(phyloPCA_all_means_prop[[c]])
}

#Check that at the same component as the shortest category cum variance explained is >0.95
phyloPCA_all_means_cum

#Keep same number of components for each category and make data frame
phylo_pcscores_all_means_df <- data.frame(rbind(phylo_pcscores_all_means[[1]], phylo_pcscores_all_means[[2]][,1:14], phylo_pcscores_all_means[[3]][,1:14],phylo_pcscores_all_means[[4]][,1:14]))

#Add labels and other attributes to tibble as columns
group_all_list <- list()
category_all_list <- list()
genus_all_list <- list()
size_all_list <- list()

for (c in 1:length(categories_list)){   
  group_all_list[[c]] <- gdf_mean_shapes[[c]]$group
  category_all_list[[c]] <- gdf_mean_shapes[[c]]$category
  genus_all_list[[c]] <- gdf_mean_shapes[[c]]$genus
  size_all_list[[c]] <- gdf_mean_shapes[[c]]$size
}

phylo_pcscores_all_means_df <- phylo_pcscores_all_means_df %>% mutate(group = unlist(group_all_list), category = unlist(category_all_list),
                                                                      genus = unlist(genus_all_list), size = unlist(size_all_list),
                                                                      module = rep("skull", times = length(unlist(group_all_list))))
glimpse(phylo_pcscores_all_means_df)

#Make dataframe with sdev values
phyloPCA_all_means_prop_df <- as.data.frame(cbind(phyloPCA_all_means_prop[[1]],phyloPCA_all_means_prop[[2]][1:14],
                                               phyloPCA_all_means_prop[[3]][1:14],phyloPCA_all_means_prop[[4]][1:14]))
colnames(phyloPCA_all_means_prop_df) <- categories_list

multiplier <- 100

phyloPCA_all_means_prop_df <- phyloPCA_all_means_prop_df %>% mutate(across(everything(), ~ . * multiplier, .names = "{.col}_100"))
phyloPCA_all_means_prop_df <- phyloPCA_all_means_prop_df %>% select(5:8) %>%
                          mutate(mean_100 = rowMeans(.))
phyloPCA_all_means_prop_df <- as.data.frame(t(phyloPCA_all_means_prop_df))
phyloPCA_all_means_prop_df$category <- c(categories_list, "mean")
phyloPCA_all_means_prop_df

#Nice PCA plot with stages and groups
phyloPCA_all_means_ggplot <- ggplot(phylo_pcscores_all_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = group, fill = group))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 80)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Groups", labels = levels(groups), #copy from as.factor(genera)
                      values = mypalette_groups, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_all_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_phylo_means_all_category_myst <- phylo_pcscores_all_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_phylo_means_all_category_odont <- phylo_pcscores_all_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_all_means_category_ggplot <- ggplot(phylo_pcscores_all_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_all_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle("Whole skull")+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_all_means_category_ggplot

#Divide by category
#Make categories labels for plots
categories_labels <- as_labeller(c("1-early" = "Early Fetus", "2-late/new" = "Late Fetus/Neonate", "3-immature" = "Juvenile" ,"4-adult"   = "Adult"))
  
#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_all_means_div_category_ggplot <- ggplot(phylo_pcscores_all_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_all_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  facet_wrap(vars(category), scales = "free", labeller = categories_labels)+
  theme_bw()+
  ggtitle("Whole skull")+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_all_means_div_category_ggplot

#Separate plots for each category

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_all_means_category_early_ggplot <- ggplot(phylo_pcscores_all_means_df[phylo_pcscores_all_means_df$category == categories_list[1],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[1], fill = mypalette_category[1])+
  geom_polygon(data = hulls_phylo_means_all_category_myst[hulls_phylo_means_all_category_myst$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[1], fill = mypalette_category[1], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont[hulls_phylo_means_all_category_odont$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[1], fill = mypalette_category[1],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[1])+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[1], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[1], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_all_means_category_late_new_ggplot <- ggplot(phylo_pcscores_all_means_df[phylo_pcscores_all_means_df$category == categories_list[2],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[2], fill = mypalette_category[2])+
  geom_polygon(data = hulls_phylo_means_all_category_myst[hulls_phylo_means_all_category_myst$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[2], fill = mypalette_category[2], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont[hulls_phylo_means_all_category_odont$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[2], fill = mypalette_category[2],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[2])+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[2], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[2], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_all_means_category_immature_ggplot <- ggplot(phylo_pcscores_all_means_df[phylo_pcscores_all_means_df$category == categories_list[3],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[3], fill = mypalette_category[3])+
  geom_polygon(data = hulls_phylo_means_all_category_myst[hulls_phylo_means_all_category_myst$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[3], fill = mypalette_category[3], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont[hulls_phylo_means_all_category_odont$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[3], fill = mypalette_category[3],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[3])+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[3], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[3], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_all_means_category_adult_ggplot <- ggplot(phylo_pcscores_all_means_df[phylo_pcscores_all_means_df$category == categories_list[4],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[4], fill = mypalette_category[4])+
  geom_polygon(data = hulls_phylo_means_all_category_myst[hulls_phylo_means_all_category_myst$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[4], fill = mypalette_category[4], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_all_category_odont[hulls_phylo_means_all_category_odont$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[4], fill = mypalette_category[4],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[4])+
  xlab(paste0("PC 1 (",round(phyloPCA_all_means_prop_df$Comp1[4], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_all_means_prop_df$Comp2[4], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")


plotPC1 <- ggarrange(phyloPCA_all_means_category_early_ggplot, phyloPCA_all_means_category_late_new_ggplot, phyloPCA_all_means_category_immature_ggplot,phyloPCA_all_means_category_adult_ggplot,
          nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
plotPC1 <- annotate_figure(plotPC1, top = text_grob("Whole skull", face = "bold", size = 17))
plotPC1

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_phyloPC1all_means_size <- lm(Comp1 ~ size, data = phylo_pcscores_all_means_df)
reg_phyloPC2all_means_size <- lm(Comp2 ~ size, data = phylo_pcscores_all_means_df)

#View results and p-value
summary(reg_phyloPC1all_means_size)
summary(reg_phyloPC2all_means_size)
anova(reg_phyloPC1all_means_size)
anova(reg_phyloPC2all_means_size)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2all_means_size_lm.txt")
print("PC1")
summary(reg_phyloPC1all_means_size)
anova(reg_phyloPC1all_means_size)
print("PC2")
summary(reg_phyloPC2all_means_size)
anova(reg_phyloPC2all_means_size)
sink() 

#Calculate regression for each component for size and category
reg_phyloPC1all_means_size_cat <- lm(Comp1 ~ size * category, data = phylo_pcscores_all_means_df)
reg_phyloPC2all_means_size_cat <- lm(Comp2 ~ size * category, data = phylo_pcscores_all_means_df)

#View results and p-value
summary(reg_phyloPC1all_means_size_cat)
summary(reg_phyloPC2all_means_size_cat)
anova(reg_phyloPC1all_means_size_cat)
anova(reg_phyloPC2all_means_size_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2all_means_size_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1all_means_size_cat)
anova(reg_phyloPC1all_means_size_cat)
print("PC2")
summary(reg_phyloPC2all_means_size_cat)
anova(reg_phyloPC2all_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_phyloPC1all_means_category <- lm(Comp1 ~ category, data = phylo_pcscores_all_means_df)
reg_phyloPC2all_means_category <- lm(Comp2 ~ category, data = phylo_pcscores_all_means_df)

#View results and p-value
summary(reg_phyloPC1all_means_category)
summary(reg_phyloPC2all_means_category)
anova(reg_phyloPC1all_means_category)
anova(reg_phyloPC2all_means_category)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2all_means_category_lm.txt")
print("PC1")
summary(reg_phyloPC1all_means_category)
anova(reg_phyloPC1all_means_category)
print("PC2")
summary(reg_phyloPC2all_means_category)
anova(reg_phyloPC2all_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_phyloPC1all_means_group <- lm(Comp1 ~ group, data = phylo_pcscores_all_means_df)
reg_phyloPC2all_means_group <- lm(Comp2 ~ group, data = phylo_pcscores_all_means_df)

#View results and p-value
summary(reg_phyloPC1all_means_group)
summary(reg_phyloPC2all_means_group)
anova(reg_phyloPC1all_means_group)
anova(reg_phyloPC2all_means_group)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2all_means_group_lm.txt")
print("PC1")
summary(reg_phyloPC1all_means_group)
anova(reg_phyloPC1all_means_group)
print("PC2")
summary(reg_phyloPC2all_means_group)
anova(reg_phyloPC2all_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_phyloPC1all_means_group_cat <- lm(Comp1 ~ group * category, data = phylo_pcscores_all_means_df)
reg_phyloPC2all_means_group_cat <- lm(Comp2 ~ group * category, data = phylo_pcscores_all_means_df)

#View results and p-value
summary(reg_phyloPC1all_means_group_cat)
summary(reg_phyloPC2all_means_group_cat)
anova(reg_phyloPC1all_means_group_cat)
anova(reg_phyloPC2all_means_group_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2all_means_group_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1all_means_group_cat)
anova(reg_phyloPC1all_means_group_cat)
print("PC2")
summary(reg_phyloPC2all_means_group_cat)
anova(reg_phyloPC2all_means_group_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/5b-phyloPCA means/PC1-2_all_means_lm.txt")
print("PC1")
anova(reg_phyloPC1all_means_size)
anova(reg_phyloPC1all_means_size_cat)
anova(reg_phyloPC1all_means_category)
anova(reg_phyloPC1all_means_group)
anova(reg_phyloPC1all_means_group_cat)
print("PC2")
anova(reg_phyloPC2all_means_size)
anova(reg_phyloPC2all_means_size_cat)
anova(reg_phyloPC2all_means_category)
anova(reg_phyloPC2all_means_group)
anova(reg_phyloPC2all_means_group_cat)
sink()


##Rostrum ----
#Run PCA on rostrum dataset

phyloPCA_rostrum_means <- list()

for (c in 1:length(categories_list)){    
  phyloPCA_rostrum_means[[c]] <- gm.prcomp(gdf_mean_rostrum_shapes[[c]]$coords, phy = trees_list[[c]], 
                                       GLS = TRUE, transform = TRUE)
}

#List of PC components and proportion of variation
phyloPCA_rostrum_means

#Save PCA results to file
sink("Output/5b-phyloPCA means/phyloPCA_rostrum_means_components.txt")
print("phyloPCA complete dataset")
print(phyloPCA_rostrum_means)
sink() 

##View plot
for (c in 1:length(categories_list)){    
  plot(phyloPCA_rostrum_means[[c]], phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA rostrum data mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
       col = "deeppink",   #border of points
       bg = "deeppink",    #fill of points
       cex = 1,            #size of points (1=regular)
       font.main = 2)       #bold font for title
  #Add quick labels to plot
  text(x = phyloPCA_rostrum_means[[c]]$x[,1], y = phyloPCA_rostrum_means$x[,2], labels = rownames(phyloPCA_rostrum_means[[c]]$x), 
       pos = 1,       #position relative to data point
       offset = 0.5,  #distance from data point
       cex = 0.75)    #font size (1=regular)
}

##View plot
for (c in 1:length(categories_list)){   
  plot(phyloPCA_rostrum_means[[c]], axis1 = 1, axis2 = 3, phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA rostrum data mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
       col = "orange",   #border of points
       bg = "orange",    #fill of points
       cex = 1,            #size of points (1=regular)
       font.main = 2)       #bold font for title
  #Add quick labels to plot
  text(x = phyloPCA_rostrum_means[[c]]$x[,1], y = phyloPCA_rostrum_means[[c]]$x[,3], labels = rownames(phyloPCA_rostrum_means[[c]]$x), 
       pos = 1,       #position relative to data point
       offset = 0.5,  #distance from data point
       cex = 0.75)    #font size (1=regular)
}

#Save PC scores as object to use later
phylo_pcscores_rostrum_means <- list()

for (c in 1:length(categories_list)){   
  phylo_pcscores_rostrum_means[[c]] <- phyloPCA_rostrum_means[[c]]$x 
}

#Save shapes of extremes for axes used in plot
phyloPC1min_rostrum_means <- list()
phyloPC1max_rostrum_means <- list()
phyloPC2min_rostrum_means <- list()
phyloPC2max_rostrum_means <- list()

for (c in 1:length(categories_list)){   
  phyloPC1min_rostrum_means[[c]] <- phyloPCA_rostrum_means[[c]][["shapes"]][["shapes.comp1"]][["min"]]
  phyloPC1max_rostrum_means[[c]] <- phyloPCA_rostrum_means[[c]][["shapes"]][["shapes.comp1"]][["max"]] 
  phyloPC2min_rostrum_means[[c]] <- phyloPCA_rostrum_means[[c]][["shapes"]][["shapes.comp2"]][["min"]] 
  phyloPC2max_rostrum_means[[c]] <- phyloPCA_rostrum_means[[c]][["shapes"]][["shapes.comp2"]][["max"]] 
}

#Show 3D deformation from mean with points overlay, do this for rostrum 4 extremes - using spheres3D for points
#Early
#phyloPC1min colors
#spheres3d(mean_rostrum_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_rostrum_means_points <- spheres3d(phyloPC1min_rostrum_means[[1]], radius=.0008, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means_early.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means1_early.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_rostrum_means_points <- spheres3d(phyloPC1max_rostrum_means[[1]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means_early.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_rostrum_means_points <- spheres3d(phyloPC2min_rostrum_means[[1]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means1_early.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_rostrum_means_points <- spheres3d(phyloPC2max_rostrum_means[[1]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means_early.png")  
clear3d()

#Late/New
#phyloPC1min colors
#spheres3d(mean_rostrum_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_rostrum_means_points <- spheres3d(phyloPC1min_rostrum_means[[2]], radius=.0008, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means_late_new.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means1_late_new.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_rostrum_means_points <- spheres3d(phyloPC1max_rostrum_means[[2]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means_late_new.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_rostrum_means_points <- spheres3d(phyloPC2min_rostrum_means[[2]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means1_late_new.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_rostrum_means_points <- spheres3d(phyloPC2max_rostrum_means[[2]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means_late_new.png")  
clear3d()

#Immature
#phyloPC1min colors
#spheres3d(mean_rostrum_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_rostrum_means_points <- spheres3d(phyloPC1min_rostrum_means[[3]], radius=.0008, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means_immature.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means1_immature.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_rostrum_means_points <- spheres3d(phyloPC1max_rostrum_means[[3]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means_immature.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_rostrum_means_points <- spheres3d(phyloPC2min_rostrum_means[[3]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means1_immature.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_rostrum_means_points <- spheres3d(phyloPC2max_rostrum_means[[3]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means_immature.png")  
clear3d()

#Adult
#phyloPC1min colors
#spheres3d(mean_rostrum_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_rostrum_means_points <- spheres3d(phyloPC1min_rostrum_means[[4]], radius=.0008, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means_adult.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_rostrum_means1_adult.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_rostrum_means_points <- spheres3d(phyloPC1max_rostrum_means[[4]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_rostrum_means_adult.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_rostrum_means_points <- spheres3d(phyloPC2min_rostrum_means[[4]], radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_rostrum_means1_adult.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_rostrum_means_points <- spheres3d(phyloPC2max_rostrum_means[[4]], radius=.0006, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_rostrum_means_adult.png")  
clear3d()

##Make better PCA plot using ggplot
#Check how many components explain >95% variation and compare with total components each category
phyloPCA_rostrum_means_prop <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_rostrum_means_prop[[c]] <- phyloPCA_rostrum_means[[c]]$sdev^2/ sum(phyloPCA_rostrum_means[[c]]$sdev^2)
}

phyloPCA_rostrum_means_cum <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_rostrum_means_cum[[c]] <- cumsum(phyloPCA_rostrum_means_prop[[c]])
}

#Check that at the same component as the shortest category cum variance explained is >0.95
phyloPCA_rostrum_means_cum

#Keep same number of components for each category and make data frame
phylo_pcscores_rostrum_means_df <- data.frame(rbind(phylo_pcscores_rostrum_means[[1]], phylo_pcscores_rostrum_means[[2]][,1:14], phylo_pcscores_rostrum_means[[3]][,1:14],phylo_pcscores_rostrum_means[[4]][,1:14]))

#Add labels and other attributes to tibble as columns
group_rostrum_list <- list()
category_rostrum_list <- list()
genus_rostrum_list <- list()
size_rostrum_list <- list()

for (c in 1:length(categories_list)){   
  group_rostrum_list[[c]] <- gdf_mean_rostrum_shapes[[c]]$group
  category_rostrum_list[[c]] <- gdf_mean_rostrum_shapes[[c]]$category
  genus_rostrum_list[[c]] <- gdf_mean_rostrum_shapes[[c]]$genus
  size_rostrum_list[[c]] <- gdf_mean_rostrum_shapes[[c]]$size
}

phylo_pcscores_rostrum_means_df <- phylo_pcscores_rostrum_means_df %>% mutate(group = unlist(group_rostrum_list), category = unlist(category_rostrum_list),
                                                                      genus = unlist(genus_rostrum_list), size = unlist(size_rostrum_list),
                                                                      module = rep("rostrum", times = length(unlist(group_rostrum_list))))
glimpse(phylo_pcscores_rostrum_means_df)

#Make dataframe with sdev values
phyloPCA_rostrum_means_prop_df <- as.data.frame(cbind(phyloPCA_rostrum_means_prop[[1]],phyloPCA_rostrum_means_prop[[2]][1:14],
                                                  phyloPCA_rostrum_means_prop[[3]][1:14],phyloPCA_rostrum_means_prop[[4]][1:14]))
colnames(phyloPCA_rostrum_means_prop_df) <- categories_list

phyloPCA_rostrum_means_prop_df <- phyloPCA_rostrum_means_prop_df %>% mutate(across(everything(), ~ . * multiplier, .names = "{.col}_100"))
phyloPCA_rostrum_means_prop_df <- phyloPCA_rostrum_means_prop_df %>% select(5:8) %>%
  mutate(mean_rostrum_100 = rowMeans(.))
phyloPCA_rostrum_means_prop_df <- as.data.frame(t(phyloPCA_rostrum_means_prop_df))
phyloPCA_rostrum_means_prop_df$category <- c(categories_list, "mean")
phyloPCA_rostrum_means_prop_df

#Nice PCA plot with stages and groups
phyloPCA_rostrum_means_ggplot <- ggplot(phylo_pcscores_rostrum_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = group, fill = group))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 80)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Groups", labels = levels(groups), #copy from as.factor(genera)
                      values = mypalette_groups, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_rostrum_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_phylo_means_rostrum_category_myst <- phylo_pcscores_rostrum_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_phylo_means_rostrum_category_odont <- phylo_pcscores_rostrum_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_rostrum_means_category_ggplot <- ggplot(phylo_pcscores_rostrum_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle("Rostrum")+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_rostrum_means_category_ggplot

#Divide by category
#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_rostrum_means_div_category_ggplot <- ggplot(phylo_pcscores_rostrum_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  facet_wrap(vars(category), scales = "free", labeller = categories_labels)+
  theme_bw()+
  ggtitle("Rostrum")+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_rostrum_means_div_category_ggplot

#Separate plots for each category

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_rostrum_means_category_early_ggplot <- ggplot(phylo_pcscores_rostrum_means_df[phylo_pcscores_rostrum_means_df$category == categories_list[1],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[1], fill = mypalette_category[1])+
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst[hulls_phylo_means_rostrum_category_myst$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[1], fill = mypalette_category[1], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont[hulls_phylo_means_rostrum_category_odont$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[1], fill = mypalette_category[1],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[1])+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[1], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[1], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_rostrum_means_category_late_new_ggplot <- ggplot(phylo_pcscores_rostrum_means_df[phylo_pcscores_rostrum_means_df$category == categories_list[2],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[2], fill = mypalette_category[2])+
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst[hulls_phylo_means_rostrum_category_myst$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[2], fill = mypalette_category[2], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont[hulls_phylo_means_rostrum_category_odont$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[2], fill = mypalette_category[2],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[2])+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[2], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[2], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_rostrum_means_category_immature_ggplot <- ggplot(phylo_pcscores_rostrum_means_df[phylo_pcscores_rostrum_means_df$category == categories_list[3],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[3], fill = mypalette_category[3])+
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst[hulls_phylo_means_rostrum_category_myst$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[3], fill = mypalette_category[3], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont[hulls_phylo_means_rostrum_category_odont$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[3], fill = mypalette_category[3],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[3])+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[3], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[3], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_rostrum_means_category_adult_ggplot <- ggplot(phylo_pcscores_rostrum_means_df[phylo_pcscores_rostrum_means_df$category == categories_list[4],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[4], fill = mypalette_category[4])+
  geom_polygon(data = hulls_phylo_means_rostrum_category_myst[hulls_phylo_means_rostrum_category_myst$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[4], fill = mypalette_category[4], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_rostrum_category_odont[hulls_phylo_means_rostrum_category_odont$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[4], fill = mypalette_category[4],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[4])+
  xlab(paste0("PC 1 (",round(phyloPCA_rostrum_means_prop_df$Comp1[4], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_rostrum_means_prop_df$Comp2[4], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")


plotPC2 <- ggarrange(phyloPCA_rostrum_means_category_early_ggplot, phyloPCA_rostrum_means_category_late_new_ggplot, phyloPCA_rostrum_means_category_immature_ggplot,phyloPCA_rostrum_means_category_adult_ggplot,
                     nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
plotPC2 <- annotate_figure(plotPC2, top = text_grob("Rostrum", face = "bold", size = 17))
plotPC2

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_phyloPC1rostrum_means_size <- lm(Comp1 ~ size, data = phylo_pcscores_rostrum_means_df)
reg_phyloPC2rostrum_means_size <- lm(Comp2 ~ size, data = phylo_pcscores_rostrum_means_df)

#View results and p-value
summary(reg_phyloPC1rostrum_means_size)
summary(reg_phyloPC2rostrum_means_size)
anova(reg_phyloPC1rostrum_means_size)
anova(reg_phyloPC2rostrum_means_size)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2rostrum_means_size_lm.txt")
print("PC1")
summary(reg_phyloPC1rostrum_means_size)
anova(reg_phyloPC1rostrum_means_size)
print("PC2")
summary(reg_phyloPC2rostrum_means_size)
anova(reg_phyloPC2rostrum_means_size)
sink() 

#Calculate regression for each component for size and category
reg_phyloPC1rostrum_means_size_cat <- lm(Comp1 ~ size * category, data = phylo_pcscores_rostrum_means_df)
reg_phyloPC2rostrum_means_size_cat <- lm(Comp2 ~ size * category, data = phylo_pcscores_rostrum_means_df)

#View results and p-value
summary(reg_phyloPC1rostrum_means_size_cat)
summary(reg_phyloPC2rostrum_means_size_cat)
anova(reg_phyloPC1rostrum_means_size_cat)
anova(reg_phyloPC2rostrum_means_size_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2rostrum_means_size_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1rostrum_means_size_cat)
anova(reg_phyloPC1rostrum_means_size_cat)
print("PC2")
summary(reg_phyloPC2rostrum_means_size_cat)
anova(reg_phyloPC2rostrum_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_phyloPC1rostrum_means_category <- lm(Comp1 ~ category, data = phylo_pcscores_rostrum_means_df)
reg_phyloPC2rostrum_means_category <- lm(Comp2 ~ category, data = phylo_pcscores_rostrum_means_df)

#View results and p-value
summary(reg_phyloPC1rostrum_means_category)
summary(reg_phyloPC2rostrum_means_category)
anova(reg_phyloPC1rostrum_means_category)
anova(reg_phyloPC2rostrum_means_category)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2rostrum_means_category_lm.txt")
print("PC1")
summary(reg_phyloPC1rostrum_means_category)
anova(reg_phyloPC1rostrum_means_category)
print("PC2")
summary(reg_phyloPC2rostrum_means_category)
anova(reg_phyloPC2rostrum_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_phyloPC1rostrum_means_group <- lm(Comp1 ~ group, data = phylo_pcscores_rostrum_means_df)
reg_phyloPC2rostrum_means_group <- lm(Comp2 ~ group, data = phylo_pcscores_rostrum_means_df)

#View results and p-value
summary(reg_phyloPC1rostrum_means_group)
summary(reg_phyloPC2rostrum_means_group)
anova(reg_phyloPC1rostrum_means_group)
anova(reg_phyloPC2rostrum_means_group)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2rostrum_means_group_lm.txt")
print("PC1")
summary(reg_phyloPC1rostrum_means_group)
anova(reg_phyloPC1rostrum_means_group)
print("PC2")
summary(reg_phyloPC2rostrum_means_group)
anova(reg_phyloPC2rostrum_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_phyloPC1rostrum_means_group_cat <- lm(Comp1 ~ group * category, data = phylo_pcscores_rostrum_means_df)
reg_phyloPC2rostrum_means_group_cat <- lm(Comp2 ~ group * category, data = phylo_pcscores_rostrum_means_df)

#View results and p-value
summary(reg_phyloPC1rostrum_means_group_cat)
summary(reg_phyloPC2rostrum_means_group_cat)
anova(reg_phyloPC1rostrum_means_group_cat)
anova(reg_phyloPC2rostrum_means_group_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2rostrum_means_group_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1rostrum_means_group_cat)
anova(reg_phyloPC1rostrum_means_group_cat)
print("PC2")
summary(reg_phyloPC2rostrum_means_group_cat)
anova(reg_phyloPC2rostrum_means_group_cat)
sink() 

#Save results of rostrum regressions to 1 file
sink("Output/5b-phyloPCA means/PC1-2_rostrum_means_lm.txt")
print("PC1")
anova(reg_phyloPC1rostrum_means_size)
anova(reg_phyloPC1rostrum_means_size_cat)
anova(reg_phyloPC1rostrum_means_category)
anova(reg_phyloPC1rostrum_means_group)
anova(reg_phyloPC1rostrum_means_group_cat)
print("PC2")
anova(reg_phyloPC2rostrum_means_size)
anova(reg_phyloPC2rostrum_means_size_cat)
anova(reg_phyloPC2rostrum_means_category)
anova(reg_phyloPC2rostrum_means_group)
anova(reg_phyloPC2rostrum_means_group_cat)
sink()

##Braincase ----
#Run PCA on braincase dataset

phyloPCA_braincase_means <- list()

for (c in 1:length(categories_list)){    
  phyloPCA_braincase_means[[c]] <- gm.prcomp(gdf_mean_braincase_shapes[[c]]$coords, phy = trees_list[[c]], 
                                           GLS = TRUE, transform = TRUE)
}

#List of PC components and proportion of variation
phyloPCA_braincase_means

#Save PCA results to file
sink("Output/5b-phyloPCA means/phyloPCA_braincase_means_components.txt")
print("phyloPCA complete dataset")
print(phyloPCA_braincase_means)
sink() 

##View plot
for (c in 1:length(categories_list)){    
  plot(phyloPCA_braincase_means[[c]], phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA braincase data mean shapes - PC1-PC2",  pch = 21, #title and type of point to be used
       col = "deeppink",   #border of points
       bg = "deeppink",    #fill of points
       cex = 1,            #size of points (1=regular)
       font.main = 2)       #bold font for title
  #Add quick labels to plot
  text(x = phyloPCA_braincase_means[[c]]$x[,1], y = phyloPCA_braincase_means$x[,2], labels = rownames(phyloPCA_braincase_means[[c]]$x), 
       pos = 1,       #position relative to data point
       offset = 0.5,  #distance from data point
       cex = 0.75)    #font size (1=regular)
}

##View plot
for (c in 1:length(categories_list)){   
  plot(phyloPCA_braincase_means[[c]], axis1 = 1, axis2 = 3, phylo = T, phylo.par = list(node.labels = F), main = "phyloPCA braincase data mean shapes - PC1-PC3",  pch = 21, #title and type of point to be used
       col = "orange",   #border of points
       bg = "orange",    #fill of points
       cex = 1,            #size of points (1=regular)
       font.main = 2)       #bold font for title
  #Add quick labels to plot
  text(x = phyloPCA_braincase_means[[c]]$x[,1], y = phyloPCA_braincase_means[[c]]$x[,3], labels = rownames(phyloPCA_braincase_means[[c]]$x), 
       pos = 1,       #position relative to data point
       offset = 0.5,  #distance from data point
       cex = 0.75)    #font size (1=regular)
}

#Save PC scores as object to use later
phylo_pcscores_braincase_means <- list()

for (c in 1:length(categories_list)){   
  phylo_pcscores_braincase_means[[c]] <- phyloPCA_braincase_means[[c]]$x 
}

#Save shapes of extremes for axes used in plot
phyloPC1min_braincase_means <- list()
phyloPC1max_braincase_means <- list()
phyloPC2min_braincase_means <- list()
phyloPC2max_braincase_means <- list()

for (c in 1:length(categories_list)){   
  phyloPC1min_braincase_means[[c]] <- phyloPCA_braincase_means[[c]][["shapes"]][["shapes.comp1"]][["min"]]
  phyloPC1max_braincase_means[[c]] <- phyloPCA_braincase_means[[c]][["shapes"]][["shapes.comp1"]][["max"]] 
  phyloPC2min_braincase_means[[c]] <- phyloPCA_braincase_means[[c]][["shapes"]][["shapes.comp2"]][["min"]] 
  phyloPC2max_braincase_means[[c]] <- phyloPCA_braincase_means[[c]][["shapes"]][["shapes.comp2"]][["max"]] 
}

#Show 3D deformation from mean with points overlay, do this for braincase 4 extremes - using spheres3D for points
#Early
#phyloPC1min colors
#spheres3d(mean_braincase_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_braincase_means_points <- spheres3d(phyloPC1min_braincase_means[[1]], radius=.0008, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means_early.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means1_early.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_braincase_means_points <- spheres3d(phyloPC1max_braincase_means[[1]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means_early.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_braincase_means_points <- spheres3d(phyloPC2min_braincase_means[[1]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means1_early.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_braincase_means_points <- spheres3d(phyloPC2max_braincase_means[[1]], radius=.002, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means1_early.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means_early.png")  
clear3d()

#Late/New
#phyloPC1min colors
#spheres3d(mean_braincase_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_braincase_means_points <- spheres3d(phyloPC1min_braincase_means[[2]], radius=.0008, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means_late_new.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means1_late_new.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_braincase_means_points <- spheres3d(phyloPC1max_braincase_means[[2]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means_late_new.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_braincase_means_points <- spheres3d(phyloPC2min_braincase_means[[2]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means1_late_new.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_braincase_means_points <- spheres3d(phyloPC2max_braincase_means[[2]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means1_late_new.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means_late_new.png")  
clear3d()

#Immature
#phyloPC1min colors
#spheres3d(mean_braincase_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_braincase_means_points <- spheres3d(phyloPC1min_braincase_means[[3]], radius=.0008, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means_immature.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means1_immature.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_braincase_means_points <- spheres3d(phyloPC1max_braincase_means[[3]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means_immature.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_braincase_means_points <- spheres3d(phyloPC2min_braincase_means[[3]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means1_immature.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_braincase_means_points <- spheres3d(phyloPC2max_braincase_means[[3]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means1_immature.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means_immature.png")  
clear3d()

#Adult
#phyloPC1min colors
#spheres3d(mean_braincase_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
phyloPC1min_braincase_means_points <- spheres3d(phyloPC1min_braincase_means[[4]], radius=.0008, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means_adult.png")
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1min_braincase_means1_adult.png") 
clear3d()

#phyloPC1max colors
phyloPC1max_braincase_means_points <- spheres3d(phyloPC1max_braincase_means[[4]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC1max_braincase_means_adult.png") 
clear3d()

#phyloPC2min colors
phyloPC2min_braincase_means_points <- spheres3d(phyloPC2min_braincase_means[[4]], radius=.001, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2min_braincase_means1_adult.png") 
clear3d()

#phyloPC2max colors
phyloPC2max_braincase_means_points <- spheres3d(phyloPC2max_braincase_means[[4]], radius=.0006, color = col_modules[braincase])
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means1_adult.png") 
rgl.snapshot(filename = "Output/5b-phyloPCA means/min-max shapes/phyloPC2max_braincase_means_adult.png")  
clear3d()

##Make better PCA plot using ggplot
#Check how many components explain >95% variation and compare with total components each category
phyloPCA_braincase_means_prop <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_braincase_means_prop[[c]] <- phyloPCA_braincase_means[[c]]$sdev^2/ sum(phyloPCA_braincase_means[[c]]$sdev^2)
}

phyloPCA_braincase_means_cum <- list()

for (c in 1:length(categories_list)){   
  phyloPCA_braincase_means_cum[[c]] <- cumsum(phyloPCA_braincase_means_prop[[c]])
}

#Check that at the same component as the shortest category cum variance explained is >0.95
phyloPCA_braincase_means_cum

#Keep same number of components for each category and make data frame
phylo_pcscores_braincase_means_df <- data.frame(rbind(phylo_pcscores_braincase_means[[1]], phylo_pcscores_braincase_means[[2]][,1:14], phylo_pcscores_braincase_means[[3]][,1:14],phylo_pcscores_braincase_means[[4]][,1:14]))

#Add labels and other attributes to tibble as columns
group_braincase_list <- list()
category_braincase_list <- list()
genus_braincase_list <- list()
size_braincase_list <- list()

for (c in 1:length(categories_list)){   
  group_braincase_list[[c]] <- gdf_mean_braincase_shapes[[c]]$group
  category_braincase_list[[c]] <- gdf_mean_braincase_shapes[[c]]$category
  genus_braincase_list[[c]] <- gdf_mean_braincase_shapes[[c]]$genus
  size_braincase_list[[c]] <- gdf_mean_braincase_shapes[[c]]$size
}

phylo_pcscores_braincase_means_df <- phylo_pcscores_braincase_means_df %>% mutate(group = unlist(group_braincase_list), category = unlist(category_braincase_list),
                                                                              genus = unlist(genus_braincase_list), size = unlist(size_braincase_list),
                                                                              module = rep("braincase", times = length(unlist(group_braincase_list))))
glimpse(phylo_pcscores_braincase_means_df)

#Make dataframe with sdev values
phyloPCA_braincase_means_prop_df <- as.data.frame(cbind(phyloPCA_braincase_means_prop[[1]],phyloPCA_braincase_means_prop[[2]][1:14],
                                                      phyloPCA_braincase_means_prop[[3]][1:14],phyloPCA_braincase_means_prop[[4]][1:14]))
colnames(phyloPCA_braincase_means_prop_df) <- categories_list

phyloPCA_braincase_means_prop_df <- phyloPCA_braincase_means_prop_df %>% mutate(across(everything(), ~ . * multiplier, .names = "{.col}_100"))
phyloPCA_braincase_means_prop_df <- phyloPCA_braincase_means_prop_df %>% select(5:8) %>%
  mutate(mean_braincase_100 = rowMeans(.))
phyloPCA_braincase_means_prop_df <- as.data.frame(t(phyloPCA_braincase_means_prop_df))
phyloPCA_braincase_means_prop_df$category <- c(categories_list, "mean")
phyloPCA_braincase_means_prop_df

#Nice PCA plot with stages and groups
phyloPCA_braincase_means_ggplot <- ggplot(phylo_pcscores_braincase_means_df, aes(x = Comp1, y = Comp2, label = genus, colour = group, fill = group))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 80)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Groups", labels = levels(groups), #copy from as.factor(genera)
                      values = mypalette_groups, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0)))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_braincase_means_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_phylo_means_braincase_category_myst <- phylo_pcscores_braincase_means_df %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_phylo_means_braincase_category_odont <- phylo_pcscores_braincase_means_df %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_braincase_means_category_ggplot <- ggplot(phylo_pcscores_braincase_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_braincase_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle("Braincase")+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
phyloPCA_braincase_means_category_ggplot

#Divide by category
#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_braincase_means_div_category_ggplot <- ggplot(phylo_pcscores_braincase_means_df, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_phylo_means_braincase_category_myst, aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont, aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  facet_wrap(vars(category), scales = "free", labeller = categories_labels)+
  theme_bw()+
  ggtitle("Braincase")+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[5], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[5], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_braincase_means_div_category_ggplot

#Separate plots for each category

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
phyloPCA_braincase_means_category_early_ggplot <- ggplot(phylo_pcscores_braincase_means_df[phylo_pcscores_braincase_means_df$category == categories_list[1],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[1], fill = mypalette_category[1])+
  geom_polygon(data = hulls_phylo_means_braincase_category_myst[hulls_phylo_means_braincase_category_myst$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[1], fill = mypalette_category[1], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont[hulls_phylo_means_braincase_category_odont$category == categories_list[1],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[1], fill = mypalette_category[1],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[1])+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[1], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[1], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_braincase_means_category_late_new_ggplot <- ggplot(phylo_pcscores_braincase_means_df[phylo_pcscores_braincase_means_df$category == categories_list[2],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[2], fill = mypalette_category[2])+
  geom_polygon(data = hulls_phylo_means_braincase_category_myst[hulls_phylo_means_braincase_category_myst$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[2], fill = mypalette_category[2], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont[hulls_phylo_means_braincase_category_odont$category == categories_list[2],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[2], fill = mypalette_category[2],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[2])+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[2], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[2], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_braincase_means_category_immature_ggplot <- ggplot(phylo_pcscores_braincase_means_df[phylo_pcscores_braincase_means_df$category == categories_list[3],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[3], fill = mypalette_category[3])+
  geom_polygon(data = hulls_phylo_means_braincase_category_myst[hulls_phylo_means_braincase_category_myst$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[3], fill = mypalette_category[3], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont[hulls_phylo_means_braincase_category_odont$category == categories_list[3],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[3], fill = mypalette_category[3],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[3])+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[3], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[3], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")

phyloPCA_braincase_means_category_adult_ggplot <- ggplot(phylo_pcscores_braincase_means_df[phylo_pcscores_braincase_means_df$category == categories_list[4],], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = group), colour = mypalette_category[4], fill = mypalette_category[4])+
  geom_polygon(data = hulls_phylo_means_braincase_category_myst[hulls_phylo_means_braincase_category_myst$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 1,
               alpha = .1, linewidth = 0.7,colour = mypalette_category[4], fill = mypalette_category[4], show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_phylo_means_braincase_category_odont[hulls_phylo_means_braincase_category_odont$category == categories_list[4],], aes(x = x, y = y, fill = category), linetype = 2,
               alpha = .1, linewidth = 0.7, colour = mypalette_category[4], fill = mypalette_category[4],show.legend = FALSE)+ #colored hulls with transparency
  scale_shape_manual(name = "Groups", labels = levels(groups),
                     values = shapes)+
  theme_bw()+
  ggtitle(levels(categories)[4])+
  xlab(paste0("PC 1 (",round(phyloPCA_braincase_means_prop_df$Comp1[4], digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(phyloPCA_braincase_means_prop_df$Comp2[4], digits = 2),"%)"))+
  theme(plot.title = element_text(size = 15, face = 3, hjust = 0.5), legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  #Remove legend for a scale_ using guide
  guides(fill = "none")


plotPC3 <- ggarrange(phyloPCA_braincase_means_category_early_ggplot, phyloPCA_braincase_means_category_late_new_ggplot, phyloPCA_braincase_means_category_immature_ggplot,phyloPCA_braincase_means_category_adult_ggplot,
                     nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
plotPC3 <- annotate_figure(plotPC3, top = text_grob("Braincase", face = "bold", size = 17))
plotPC3

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_phyloPC1braincase_means_size <- lm(Comp1 ~ size, data = phylo_pcscores_braincase_means_df)
reg_phyloPC2braincase_means_size <- lm(Comp2 ~ size, data = phylo_pcscores_braincase_means_df)

#View results and p-value
summary(reg_phyloPC1braincase_means_size)
summary(reg_phyloPC2braincase_means_size)
anova(reg_phyloPC1braincase_means_size)
anova(reg_phyloPC2braincase_means_size)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2braincase_means_size_lm.txt")
print("PC1")
summary(reg_phyloPC1braincase_means_size)
anova(reg_phyloPC1braincase_means_size)
print("PC2")
summary(reg_phyloPC2braincase_means_size)
anova(reg_phyloPC2braincase_means_size)
sink() 

#Calculate regression for each component for size and category
reg_phyloPC1braincase_means_size_cat <- lm(Comp1 ~ size * category, data = phylo_pcscores_braincase_means_df)
reg_phyloPC2braincase_means_size_cat <- lm(Comp2 ~ size * category, data = phylo_pcscores_braincase_means_df)

#View results and p-value
summary(reg_phyloPC1braincase_means_size_cat)
summary(reg_phyloPC2braincase_means_size_cat)
anova(reg_phyloPC1braincase_means_size_cat)
anova(reg_phyloPC2braincase_means_size_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2braincase_means_size_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1braincase_means_size_cat)
anova(reg_phyloPC1braincase_means_size_cat)
print("PC2")
summary(reg_phyloPC2braincase_means_size_cat)
anova(reg_phyloPC2braincase_means_size_cat)
sink() 


#Calculate regression for each component taking category into account
reg_phyloPC1braincase_means_category <- lm(Comp1 ~ category, data = phylo_pcscores_braincase_means_df)
reg_phyloPC2braincase_means_category <- lm(Comp2 ~ category, data = phylo_pcscores_braincase_means_df)

#View results and p-value
summary(reg_phyloPC1braincase_means_category)
summary(reg_phyloPC2braincase_means_category)
anova(reg_phyloPC1braincase_means_category)
anova(reg_phyloPC2braincase_means_category)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2braincase_means_category_lm.txt")
print("PC1")
summary(reg_phyloPC1braincase_means_category)
anova(reg_phyloPC1braincase_means_category)
print("PC2")
summary(reg_phyloPC2braincase_means_category)
anova(reg_phyloPC2braincase_means_category)
sink() 

#Calculate regression for each component taking group into account
reg_phyloPC1braincase_means_group <- lm(Comp1 ~ group, data = phylo_pcscores_braincase_means_df)
reg_phyloPC2braincase_means_group <- lm(Comp2 ~ group, data = phylo_pcscores_braincase_means_df)

#View results and p-value
summary(reg_phyloPC1braincase_means_group)
summary(reg_phyloPC2braincase_means_group)
anova(reg_phyloPC1braincase_means_group)
anova(reg_phyloPC2braincase_means_group)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2braincase_means_group_lm.txt")
print("PC1")
summary(reg_phyloPC1braincase_means_group)
anova(reg_phyloPC1braincase_means_group)
print("PC2")
summary(reg_phyloPC2braincase_means_group)
anova(reg_phyloPC2braincase_means_group)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_phyloPC1braincase_means_group_cat <- lm(Comp1 ~ group * category, data = phylo_pcscores_braincase_means_df)
reg_phyloPC2braincase_means_group_cat <- lm(Comp2 ~ group * category, data = phylo_pcscores_braincase_means_df)

#View results and p-value
summary(reg_phyloPC1braincase_means_group_cat)
summary(reg_phyloPC2braincase_means_group_cat)
anova(reg_phyloPC1braincase_means_group_cat)
anova(reg_phyloPC2braincase_means_group_cat)

#Save results of significant regression to file
sink("Output/5b-phyloPCA means/PC1-2braincase_means_group_cat_lm.txt")
print("PC1")
summary(reg_phyloPC1braincase_means_group_cat)
anova(reg_phyloPC1braincase_means_group_cat)
print("PC2")
summary(reg_phyloPC2braincase_means_group_cat)
anova(reg_phyloPC2braincase_means_group_cat)
sink() 

#Save results of braincase regressions to 1 file
sink("Output/5b-phyloPCA means/PC1-2_braincase_means_lm.txt")
print("PC1")
anova(reg_phyloPC1braincase_means_size)
anova(reg_phyloPC1braincase_means_size_cat)
anova(reg_phyloPC1braincase_means_category)
anova(reg_phyloPC1braincase_means_group)
anova(reg_phyloPC1braincase_means_group_cat)
print("PC2")
anova(reg_phyloPC2braincase_means_size)
anova(reg_phyloPC2braincase_means_size_cat)
anova(reg_phyloPC2braincase_means_category)
anova(reg_phyloPC2braincase_means_group)
anova(reg_phyloPC2braincase_means_group_cat)
sink()


###### 

#Next - ch. 6b - Disparity analyses - phylogenetically transformed components mean shapes