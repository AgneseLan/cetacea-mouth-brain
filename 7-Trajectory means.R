#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH. 7 - Trajectory analyses rostrum vs braincase - mean shapes

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
library(abind)
library(reshape2)
library(scales)
library(mcp)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#remotes::install_github("r-lib/rray")


#TRAJECTORY ANALYSIS ROSTRUM AND BRAINCASE SEPARATE ----

##Rostrum ----

##By group and category ----
#Need multiple observations per group/category and all have to be represented (all categories in each group)
#Check in ch. 3

#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups_rostrum (Mysticeti, Odontoceti)
fit_shape_rostrum_means_group_category <- procD.lm(coords ~ group * category, iter = 999, data = gdf_mean_rostrum, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_rostrum_means_group_category)

#Save results to file
sink("Output/7-Trajectory means/fit_shape_rostrum_group_category.txt")
summary(fit_shape_rostrum_means_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_means_groups_rostrum <- trajectory.analysis(fit_shape_rostrum_means_group_category, groups = gdf_mean_rostrum$group, traj.pts = gdf_mean_rostrum$category, 
                                                       pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_means_groups_rostrum_MD <- summary(trajectory_means_groups_rostrum, show.trajectories = T, attribute = "MD") 
trajectory_means_groups_rostrum_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_means_groups_rostrum_TC <- summary(trajectory_means_groups_rostrum, show.trajectories = T, attribute = "TC", angle.type = "deg")
trajectory_means_groups_rostrum_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_means_groups_rostrum_SD <- summary(trajectory_means_groups_rostrum, show.trajectories = T, attribute = "SD") 
trajectory_means_groups_rostrum_SD 

#Save results to file
sink("Output/7-Trajectory means/trajectory_means_groups_rostrum.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_means_groups_rostrum_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_means_groups_rostrum_TC
print("Shape differences between trajectory vectors - shape")
trajectory_means_groups_rostrum_SD 

print("Summary tables")
trajectory_means_groups_rostrum_MD$summary.table
trajectory_means_groups_rostrum_TC$summary.table
trajectory_means_groups_rostrum_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_means_groups_rostrum_plot <- plot(trajectory_means_groups_rostrum, main = "Trajectories of growth mean shapes by group rostrum",  pch = shapes, #title and type of point to be used
                                             col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_rostrum
add.trajectories(trajectory_means_groups_rostrum_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= -0.25, y = -0.2, legend = levels(groups), 
       pch =  shapes, pt.bg = 1, cex = 1)


##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_means_groups_rostrum_pcscores <- as_tibble(trajectory_means_groups_rostrum_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_means_groups_rostrum_pcscores <- trajectory_means_groups_rostrum_pcscores %>% mutate(category = gdf_mean_rostrum$category,  
                                                                                                group = gdf_mean_rostrum$group, size = gdf_mean_rostrum$size,
                                                                                                module = "rostrum")
glimpse(trajectory_means_groups_rostrum_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_means_groups_rostrum_pcscores_means <- trajectory_means_groups_rostrum_pcscores %>% group_by(group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_means_groups_rostrum_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_means_groups_rostrum_pcscores_means <- trajectory_means_groups_rostrum_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_means_groups_rostrum_pcscores_means)

#Nice plot
trajectory_means_groups_rostrum_ggplot <- ggplot(trajectory_means_groups_rostrum_pcscores, aes(x = PC1, y = PC2, shape = group, group = group))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_means_groups_rostrum_pcscores_means, aes(x = x, y = y, fill = group, shape = group, alpha = category, group = group), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_means_groups_rostrum_pcscores_means, aes(x = x, y = y, colour = group, group = group, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Group", labels = levels(groups), #copy from as.factor(groups_rostrum)
                      values = mypalette_groups, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_means_groups_rostrum[["pca"]][["sdev"]][1]^2/sum(trajectory_means_groups_rostrum[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_means_groups_rostrum[["pca"]][["sdev"]][2]^2/sum(trajectory_means_groups_rostrum[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  ggtitle("Rostrum")+
  guides(shape = "none", colour = "none", fill = "none",
         alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_means_groups_rostrum_ggplot

#Add silhouettes groups_rostrum
trajectory_means_groups_rostrum_ggplot <-   
  trajectory_means_groups_rostrum_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.28, y = 0.14, ysize = 0.025, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 0.07, y = -0.18, ysize = 0.02, fill = mypalette_groups[2])
#Annotate adult and early fetus categories
trajectory_means_groups_rostrum_ggplot <-
  trajectory_means_groups_rostrum_ggplot   + 
  annotate("text", x = -0.19, y = 0.13, label = "adult", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = 0.1, y = 0.06, label = "adult", color = mypalette_groups[2], fontface = 3, size = 6)+
  annotate("text", x = -0.21, y = -0.04, label = "early\nfetus", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = -0.055, y = -0.15, label = "early\nfetus", color = mypalette_groups[2], fontface = 3, size = 6)
trajectory_means_groups_rostrum_ggplot

###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_means_groups_rostrum_length <- trajectory_means_groups_rostrum_MD[["pairwise.tables"]][["D"]]
trajectory_means_groups_rostrum_length_p <- trajectory_means_groups_rostrum_MD[["pairwise.tables"]][["P"]]
trajectory_means_groups_rostrum_direction <- trajectory_means_groups_rostrum_TC[["pairwise.tables"]][["angle"]]
trajectory_means_groups_rostrum_direction_p <- trajectory_means_groups_rostrum_TC[["pairwise.tables"]][["P"]]
trajectory_means_groups_rostrum_shape <- trajectory_means_groups_rostrum_SD[["pairwise.tables"]][["D"]]
trajectory_means_groups_rostrum_shape_rostrum_p <- trajectory_means_groups_rostrum_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_means_groups_rostrum_list <- list(trajectory_means_groups_rostrum_length, trajectory_means_groups_rostrum_length_p, trajectory_means_groups_rostrum_direction, trajectory_means_groups_rostrum_direction_p, 
                                             trajectory_means_groups_rostrum_shape, trajectory_means_groups_rostrum_shape_rostrum_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_means_groups_rostrum_list[[l]]) <- levels(groups)
  colnames(trajectory_means_groups_rostrum_list[[l]]) <- levels(groups)
}

#Save only lower triangle for each
trajectory_means_groups_rostrum_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_means_groups_rostrum_lower_tri_list[[l]] <- get_upper_tri(trajectory_means_groups_rostrum_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_means_groups_rostrum_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_means_groups_rostrum_melt[[l]] <- melt(trajectory_means_groups_rostrum_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_means_groups_rostrum_length_melt <- data.frame(trajectory_means_groups_rostrum_melt[[1]], p = trajectory_means_groups_rostrum_melt[[2]][[3]])
trajectory_means_groups_rostrum_direction_melt <- data.frame(trajectory_means_groups_rostrum_melt[[3]], p = trajectory_means_groups_rostrum_melt[[4]][[3]])
trajectory_means_groups_rostrum_shape_rostrum_melt <- data.frame(trajectory_means_groups_rostrum_melt[[5]], p = trajectory_means_groups_rostrum_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_means_groups_rostrum_length_melt <- trajectory_means_groups_rostrum_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                      p_if_sig = ifelse(sig_p, p, NA),
                                                                                                      value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_means_groups_rostrum_direction_melt <- trajectory_means_groups_rostrum_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                            p_if_sig = ifelse(sig_p, p, NA),
                                                                                                            value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_means_groups_rostrum_shape_rostrum_melt <- trajectory_means_groups_rostrum_shape_rostrum_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                    p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                    value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_means_groups_rostrum_length_melt$p_if_sig))

all(is.na(trajectory_means_groups_rostrum_direction_melt$p_if_sig))

all(is.na(trajectory_means_groups_rostrum_shape_rostrum_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_means_groups_rostrum_direction_heatmap_ggplot <- ggplot(data = trajectory_means_groups_rostrum_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low =  mypalette_seq_modules[9], high =  mypalette_seq_modules[2], mid =  mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_purple[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Direction difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 14, vjust = 0.8),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_means_groups_rostrum_direction_heatmap_ggplot

plotR4<-ggarrange(trajectory_means_groups_rostrum_direction_heatmap_ggplot, 
                  ncol = 1, nrow = 1, common.legend = F)
plotR4<-annotate_figure(plotR4, top = text_grob("Rostrum", face = "bold", size = 17))
plotR4

##Braincase ----

##By group and category ----
#Need multiple observations per group/category and all have to be represented (all categories in each group)
#Check in ch. 3

#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups_braincase (Mysticeti, Odontoceti)
fit_shape_braincase_means_group_category <- procD.lm(coords ~ group * category, iter = 999, data = gdf_mean_braincase, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_braincase_means_group_category)

#Save results to file
sink("Output/7-Trajectory means/fit_shape_braincase_group_category.txt")
summary(fit_shape_braincase_means_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_means_groups_braincase <- trajectory.analysis(fit_shape_braincase_means_group_category, groups = gdf_mean_braincase$group, traj.pts = gdf_mean_braincase$category, 
                                                       pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_means_groups_braincase_MD <- summary(trajectory_means_groups_braincase, show.trajectories = T, attribute = "MD") 
trajectory_means_groups_braincase_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_means_groups_braincase_TC <- summary(trajectory_means_groups_braincase, show.trajectories = T, attribute = "TC", angle.type = "deg")
trajectory_means_groups_braincase_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_means_groups_braincase_SD <- summary(trajectory_means_groups_braincase, show.trajectories = T, attribute = "SD") 
trajectory_means_groups_braincase_SD 

#Save results to file
sink("Output/7-Trajectory means/trajectory_means_groups_braincase.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_means_groups_braincase_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_means_groups_braincase_TC
print("Shape differences between trajectory vectors - shape")
trajectory_means_groups_braincase_SD 

print("Summary tables")
trajectory_means_groups_braincase_MD$summary.table
trajectory_means_groups_braincase_TC$summary.table
trajectory_means_groups_braincase_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_means_groups_braincase_plot <- plot(trajectory_means_groups_braincase, main = "Trajectories of growth mean shapes by group braincase",  pch = shapes, #title and type of point to be used
                                             col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_braincase
add.trajectories(trajectory_means_groups_braincase_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= -0.1, y = -0.15, legend = levels(groups), 
       pch =  shapes, pt.bg = 1, cex = 1)


##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_means_groups_braincase_pcscores <- as_tibble(trajectory_means_groups_braincase_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_means_groups_braincase_pcscores <- trajectory_means_groups_braincase_pcscores %>% mutate(category = gdf_mean_braincase$category,  
                                                                                                group = gdf_mean_braincase$group, size = gdf_mean_braincase$size,
                                                                                                module = "braincase")
glimpse(trajectory_means_groups_braincase_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_means_groups_braincase_pcscores_means <- trajectory_means_groups_braincase_pcscores %>% group_by(group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_means_groups_braincase_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_means_groups_braincase_pcscores_means <- trajectory_means_groups_braincase_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_means_groups_braincase_pcscores_means)

#Nice plot
trajectory_means_groups_braincase_ggplot <- ggplot(trajectory_means_groups_braincase_pcscores, aes(x = PC1, y = PC2, shape = group, group = group))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_means_groups_braincase_pcscores_means, aes(x = x, y = y, fill = group, shape = group, alpha = category, group = group), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_means_groups_braincase_pcscores_means, aes(x = x, y = y, colour = group, group = group, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Group", labels = levels(groups), #copy from as.factor(groups_braincase)
                      values = mypalette_groups, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_means_groups_braincase[["pca"]][["sdev"]][1]^2/sum(trajectory_means_groups_braincase[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_means_groups_braincase[["pca"]][["sdev"]][2]^2/sum(trajectory_means_groups_braincase[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  ggtitle("Braincase")+
  guides(shape = "none", colour = "none", fill = "none",
         alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_means_groups_braincase_ggplot

#Add silhouettes groups_braincase
trajectory_means_groups_braincase_ggplot <-   
  trajectory_means_groups_braincase_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.11, y = -0.1, ysize = 0.01, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 0.11, y = -0.03, ysize = 0.009, fill = mypalette_groups[2])
#Annotate adult and early fetus categories
trajectory_means_groups_braincase_ggplot <-
  trajectory_means_groups_braincase_ggplot   + 
  annotate("text", x = -0.11, y = -0.03, label = "adult", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = -0.06, y = 0.05, label = "adult", color = mypalette_groups[2], fontface = 3, size = 6)+
  annotate("text", x = -0.035, y = -0.06, label = "early\nfetus", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = 0.1, y = -0.077, label = "early\nfetus", color = mypalette_groups[2], fontface = 3, size = 6)
trajectory_means_groups_braincase_ggplot

ggarrange(trajectory_means_groups_rostrum_ggplot, trajectory_means_groups_braincase_ggplot, ncol =2 , nrow =1, common.legend = T, legend = "bottom")

###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_means_groups_braincase_length <- trajectory_means_groups_braincase_MD[["pairwise.tables"]][["D"]]
trajectory_means_groups_braincase_length_p <- trajectory_means_groups_braincase_MD[["pairwise.tables"]][["P"]]
trajectory_means_groups_braincase_direction <- trajectory_means_groups_braincase_TC[["pairwise.tables"]][["angle"]]
trajectory_means_groups_braincase_direction_p <- trajectory_means_groups_braincase_TC[["pairwise.tables"]][["P"]]
trajectory_means_groups_braincase_shape <- trajectory_means_groups_braincase_SD[["pairwise.tables"]][["D"]]
trajectory_means_groups_braincase_shape_braincase_p <- trajectory_means_groups_braincase_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_means_groups_braincase_list <- list(trajectory_means_groups_braincase_length, trajectory_means_groups_braincase_length_p, trajectory_means_groups_braincase_direction, trajectory_means_groups_braincase_direction_p, 
                                                   trajectory_means_groups_braincase_shape, trajectory_means_groups_braincase_shape_braincase_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_means_groups_braincase_list[[l]]) <- levels(groups)
  colnames(trajectory_means_groups_braincase_list[[l]]) <- levels(groups)
}

#Save only lower triangle for each
trajectory_means_groups_braincase_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_means_groups_braincase_lower_tri_list[[l]] <- get_upper_tri(trajectory_means_groups_braincase_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_means_groups_braincase_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_means_groups_braincase_melt[[l]] <- melt(trajectory_means_groups_braincase_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_means_groups_braincase_length_melt <- data.frame(trajectory_means_groups_braincase_melt[[1]], p = trajectory_means_groups_braincase_melt[[2]][[3]])
trajectory_means_groups_braincase_direction_melt <- data.frame(trajectory_means_groups_braincase_melt[[3]], p = trajectory_means_groups_braincase_melt[[4]][[3]])
trajectory_means_groups_braincase_shape_braincase_melt <- data.frame(trajectory_means_groups_braincase_melt[[5]], p = trajectory_means_groups_braincase_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_means_groups_braincase_length_melt <- trajectory_means_groups_braincase_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                  value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_means_groups_braincase_direction_melt <- trajectory_means_groups_braincase_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_means_groups_braincase_shape_braincase_melt <- trajectory_means_groups_braincase_shape_braincase_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_means_groups_braincase_length_melt$p_if_sig))

all(is.na(trajectory_means_groups_braincase_direction_melt$p_if_sig))

all(is.na(trajectory_means_groups_braincase_shape_braincase_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_means_groups_braincase_direction_heatmap_ggplot <- ggplot(data = trajectory_means_groups_braincase_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low =  mypalette_seq_modules[9], high =  mypalette_seq_modules[2], mid =  mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_purple[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Direction difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 14, vjust = 0.8),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_means_groups_braincase_direction_heatmap_ggplot

plotB4<-ggarrange(trajectory_means_groups_braincase_direction_heatmap_ggplot, 
                  ncol = 1, nrow = 1, common.legend = F)
plotB4<-annotate_figure(plotB4, top = text_grob("Braincase", face = "bold", size = 17))
plotB4

ggarrange(plotR4, plotB4, ncol = 2, nrow = 1, common.legend = F)

##Plot both modules together for each group ----
#Mysticeti

#Select rows
rows_mysticeti_traj <- which(trajectory_means_groups_rostrum_pcscores$group == "mysticeti")

#PC scores for group for each module
trajectory_means_groups_modules_pcscores_myst <- rbind(trajectory_means_groups_rostrum_pcscores[rows_mysticeti_traj,], trajectory_means_groups_braincase_pcscores[rows_mysticeti_traj,])

glimpse(trajectory_means_groups_modules_pcscores_myst)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_means_groups_modules_pcscores_means_myst <- trajectory_means_groups_modules_pcscores_myst %>% group_by(module, group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_means_groups_modules_pcscores_means_myst)

#Rename columns so they are easier to use for plot
trajectory_means_groups_modules_pcscores_means_myst <- trajectory_means_groups_modules_pcscores_means_myst %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_means_groups_modules_pcscores_means_myst)

#Factor modules for plot legend 
modules_list <- c("braincase", "rostrum")
modules_list_factor <- str_to_sentence(levels(as.factor(modules_list)))

#Palette colors and shape modules
mypalette_modules <- c(mypalette_paired[2],mypalette_paired[5])
shapes_modules <- c(23,24)

#Nice plot
trajectory_means_groups_modules_myst_ggplot <- ggplot(trajectory_means_groups_modules_pcscores_myst, aes(x = PC1, y = PC2, shape = module, group = module))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_means_groups_modules_pcscores_means_myst, aes(x = x, y = y, fill = module, shape = module, alpha = category, group = module), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_means_groups_modules_pcscores_means_myst, aes(x = x, y = y, colour = module, group = module, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Modules", labels = modules_list_factor, values = shapes_modules)+
  scale_colour_manual(name = "Modules", labels = modules_list_factor, #copy from as.factor(groups_rostrum)
                      values = mypalette_modules, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1"))+ #no percentage as different fit plotted together
  ylab(paste0("PC 2"))+
  guides(alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_means_groups_modules_myst_ggplot

#Add silhouettes groups_rostrum
trajectory_means_groups_modules_myst_ggplot <-   
  trajectory_means_groups_modules_myst_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.05, y = 0.15, ysize = 0.02, fill = "gray30")

#Annotate adult and early fetus categories
trajectory_means_groups_modules_myst_ggplot <-
  trajectory_means_groups_modules_myst_ggplot   + 
  annotate("text", x =  -0.2, y = 0.13, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = -0.12, y = -0.01, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = -0.27, y = -0.05, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)+
  annotate("text", x = -0.03, y = -0.05, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)
trajectory_means_groups_modules_myst_ggplot

#Odontoceti

#Select rows
rows_odontoceti_traj <- which(trajectory_means_groups_rostrum_pcscores$group == "odontoceti")

#PC scores for group for each module
trajectory_means_groups_modules_pcscores_odont <- rbind(trajectory_means_groups_rostrum_pcscores[rows_odontoceti_traj,], trajectory_means_groups_braincase_pcscores[rows_odontoceti_traj,])

glimpse(trajectory_means_groups_modules_pcscores_odont)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_means_groups_modules_pcscores_means_odont <- trajectory_means_groups_modules_pcscores_odont %>% group_by(module, group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_means_groups_modules_pcscores_means_odont)

#Rename columns so they are easier to use for plot
trajectory_means_groups_modules_pcscores_means_odont <- trajectory_means_groups_modules_pcscores_means_odont %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_means_groups_modules_pcscores_means_odont)

#Nice plot
trajectory_means_groups_modules_odont_ggplot <- ggplot(trajectory_means_groups_modules_pcscores_odont, aes(x = PC1, y = PC2, shape = module, group = module))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_means_groups_modules_pcscores_means_odont, aes(x = x, y = y, fill = module, shape = module, alpha = category, group = module), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_means_groups_modules_pcscores_means_odont, aes(x = x, y = y, colour = module, group = module, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Modules", labels = modules_list_factor, values = shapes_modules)+
  scale_colour_manual(name = "Modules", labels = modules_list_factor, #copy from as.factor(groups_rostrum)
                      values = mypalette_modules, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1"))+ #no percentage as different fit plotted together
  ylab(paste0("PC 2"))+
  guides(alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_means_groups_modules_odont_ggplot

#Add silhouettes groups_rostrum
trajectory_means_groups_modules_odont_ggplot <-   
  trajectory_means_groups_modules_odont_ggplot   + 
  add_phylopic(odont, alpha = 1, x = -0.07, y = 0.15, ysize = 0.025, fill = "gray40")

#Annotate adult and early fetus categories
trajectory_means_groups_modules_odont_ggplot <-
  trajectory_means_groups_modules_odont_ggplot   + 
  annotate("text", x =  -0.05, y = 0.07, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = 0.09, y = 0.05, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = 0.11, y = -0.09, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)+
  annotate("text", x = -0.045, y = -0.16, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)
trajectory_means_groups_modules_odont_ggplot

ggarrange(trajectory_means_groups_modules_myst_ggplot, trajectory_means_groups_modules_odont_ggplot, ncol = 2, nrow = 1, common.legend = T, legend = "bottom")

###Correlation test on observed distance matrices (from pairwise) by group ----
#Mysticeti
MD_myst <- rbind(trajectory_means_groups_rostrum_MD[["x"]][["LS.means"]][["obs"]][1:4,],trajectory_means_groups_braincase_MD[["x"]][["LS.means"]][["obs"]][1:4,])

rownames(MD_myst)[1:4] <- paste0(rownames(MD_myst)[1:4],".rostrum")
rownames(MD_myst)[5:8] <- paste0(rownames(MD_myst)[5:8],".braincase")

MD_myst <- as.matrix(MD_myst)

correlation_MD_modules_myst <- cor.test(MD_myst[1:4,], MD_myst[5:8,], method = "kendall")
correlation_MD_modules_myst

#Odontoceti
MD_odont <- rbind(trajectory_means_groups_rostrum_MD[["x"]][["LS.means"]][["obs"]][5:8,],trajectory_means_groups_braincase_MD[["x"]][["LS.means"]][["obs"]][5:8,])

rownames(MD_odont)[1:4] <- paste0(rownames(MD_odont)[1:4],".rostrum")
rownames(MD_odont)[5:8] <- paste0(rownames(MD_odont)[5:8],".braincase")

MD_odont <- as.matrix(MD_odont)

correlation_MD_modules_odont <- cor.test(MD_odont[1:4,], MD_odont[5:8,], method = "kendall")
correlation_MD_modules_odont

#Save results to file
sink("Output/7-Trajectory means/correlation_distances_trajectory_groups.txt")
print("Mysticeti")
correlation_MD_modules_myst
print("Odontoceti")
correlation_MD_modules_odont
sink() 

#Next - ch. 8 - Allometry analyses means