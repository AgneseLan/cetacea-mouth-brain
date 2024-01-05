
#===========================================================#
#                                                           #
#     SKULL MODULARITY - MYSTICETI & ODONTOCETI             #
#                                                           #
#===========================================================#

#CH. 7 - Trajectory analyses rostrum vs braincase

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
library(evomap)
library(rray)
library(abind)
library(reshape2)
library(scales)
library(mcp)

#require(devtools)
#instrostrum_github("JeroenSmaers/evomap")
#devtools::instrostrum_github("wabarr/ggphylomorpho")
#devtools::instrostrum_github("aphanotus/borealis")
#remotes::install_github("r-lib/rray")


#TRAJECTORY ANALYSIS ROSTRUM AND BRAINCASE SEPARATE ----

##Rostrum ----

##By group and category ----
#Need multiple observations per group/category and all have to be represented (all categories in each group)
#Check in ch. 3

#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups_rostrum (Mysticeti, Odontoceti)
fit_shape_rostrum_group_category <- procD.lm(coords ~ group * category, iter = 999, data = gdf_rostrum, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_rostrum_group_category)

#Save results to file
sink("Output/fit_shape_rostrum_group_category.txt")
summary(fit_shape_rostrum_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_groups_rostrum <- trajectory.analysis(fit_shape_rostrum_group_category, groups = gdf_rostrum$group, traj.pts = gdf_rostrum$category, 
                                                 pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_groups_rostrum_MD <- summary(trajectory_groups_rostrum, show.trajectories = TRUE, attribute = "MD") 
trajectory_groups_rostrum_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_groups_rostrum_TC <- summary(trajectory_groups_rostrum, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_groups_rostrum_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_groups_rostrum_SD <- summary(trajectory_groups_rostrum, show.trajectories = TRUE, attribute = "SD") 
trajectory_groups_rostrum_SD 

#Save results to file
sink("Output/trajectory_groups_rostrum.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_groups_rostrum_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_groups_rostrum_TC
print("Shape differences between trajectory vectors - shape")
trajectory_groups_rostrum_SD 

print("Summary tables")
trajectory_groups_rostrum_MD$summary.table
trajectory_groups_rostrum_TC$summary.table
trajectory_groups_rostrum_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_groups_rostrum_plot <- plot(trajectory_groups_rostrum, main = "Trajectories of growth by group rostrum",  pch = shapes, #title and type of point to be used
                                       col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_rostrum
add.trajectories(trajectory_groups_rostrum_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= -0.5, y = 0.3, legend = levels(groups), 
       pch =  shapes, pt.bg = 1, cex = 1)


##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_groups_rostrum_pcscores <- as_tibble(trajectory_groups_rostrum_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_groups_rostrum_pcscores <- trajectory_groups_rostrum_pcscores %>% mutate(specimens = gdf_rostrum$Id, category = gdf_rostrum$category,  
                                                                                    group = gdf_rostrum$group, size = gdf_rostrum$size)
glimpse(trajectory_groups_rostrum_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_groups_rostrum_pcscores_means <- trajectory_groups_rostrum_pcscores %>% group_by(group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_groups_rostrum_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_groups_rostrum_pcscores_means <- trajectory_groups_rostrum_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_groups_rostrum_pcscores_means)

#Nice plot
trajectory_groups_rostrum_ggplot <- ggplot(trajectory_groups_rostrum_pcscores, aes(x = PC1, y = PC2, shape = group, group = group))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_groups_rostrum_pcscores_means, aes(x = x, y = y, fill = group, shape = group, alpha = category, group = group), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_groups_rostrum_pcscores_means, aes(x = x, y = y, colour = group, group = group, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Group", labels = levels(groups), #copy from as.factor(groups_rostrum)
                      values = mypalette_groups, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_groups_rostrum[["pca"]][["sdev"]][1]^2/sum(trajectory_groups_rostrum[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_groups_rostrum[["pca"]][["sdev"]][2]^2/sum(trajectory_groups_rostrum[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  ggtitle("Rostrum")+
  guides(shape = "none", colour = "none", fill = "none",
         alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_groups_rostrum_ggplot

#Add silhouettes groups_rostrum
trajectory_groups_rostrum_ggplot <-   
  trajectory_groups_rostrum_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.3, y = -0.25, ysize = 0.04, color = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 0.1, y = 0.28, ysize = 0.035, color = mypalette_groups[2])
#Annotate adult and early fetus categories
trajectory_groups_rostrum_ggplot <-
  trajectory_groups_rostrum_ggplot   + 
  annotate("text", x = -0.175, y = -0.18, label = "adult", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = 0.18, y = -0.06, label = "adult", color = mypalette_groups[2], fontface = 3, size = 6)+
  annotate("text", x = -0.25, y = 0.05, label = "early\nfetus", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = -0.08, y = 0.15, label = "early\nfetus", color = mypalette_groups[2], fontface = 3, size = 6)
trajectory_groups_rostrum_ggplot

#Save mean shapes for each group and stage
#PC points matrix from analysis
traj_pcs_rostrum <- trajectory_groups_rostrum_plot$pc.points[,1:2]
traj_pcs_rostrum_means <- trajectory_groups_rostrum_pcscores_means[,3:4]

#User-picked spots can be anything, but it in this case evenly-spaced PCA coordinates
preds_rostrum <- shape.predictor(gdf_rostrum$coords, x= traj_pcs_rostrum, Intercept = FALSE,
                                 myst1 = as.numeric(as.vector(traj_pcs_rostrum_means[1,])),
                                 myst2 = as.numeric(as.vector(traj_pcs_rostrum_means[2,])),
                                 myst3 = as.numeric(as.vector(traj_pcs_rostrum_means[3,])),
                                 myst4 = as.numeric(as.vector(traj_pcs_rostrum_means[4,])),
                                 odont1 = as.numeric(as.vector(traj_pcs_rostrum_means[5,])),
                                 odont2 = as.numeric(as.vector(traj_pcs_rostrum_means[6,])),
                                 odont3 = as.numeric(as.vector(traj_pcs_rostrum_means[7,])),
                                 odont4 = as.numeric(as.vector(traj_pcs_rostrum_means[8,])))

#Save shapes points
myst_rostrum_early <- spheres3d(preds_rostrum$myst1, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/myst_rostrum_early.png") 
rgl.snapshot(filename = "Output/myst_rostrum_early1.png") 
clear3d()
myst_rostrum_latenew <- spheres3d(preds_rostrum$myst2, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/myst_rostrum_latenew1.png") 
rgl.snapshot(filename = "Output/myst_rostrum_latenew.png") 
clear3d()
myst_rostrum_imm <- spheres3d(preds_rostrum$myst3, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/myst_rostrum_imm.png") 
rgl.snapshot(filename = "Output/myst_rostrum_imm1.png") 
clear3d()
myst_rostrum_adult <- spheres3d(preds_rostrum$myst4, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/myst_rostrum_adult1.png") 
rgl.snapshot(filename = "Output/myst_rostrum_adult.png") 
clear3d()

odont_rostrum_early <- spheres3d(preds_rostrum$odont1, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/odont_rostrum_early.png") 
rgl.snapshot(filename = "Output/odont_rostrum_early1.png") 
clear3d()
odont_rostrum_latenew <- spheres3d(preds_rostrum$odont2, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/odont_rostrum_latenew1.png") 
rgl.snapshot(filename = "Output/odont_rostrum_latenew.png") 
clear3d()
odont_rostrum_imm <- spheres3d(preds_rostrum$odont3, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/odont_rostrum_imm.png") 
rgl.snapshot(filename = "Output/odont_rostrum_imm1.png") 
clear3d()
odont_rostrum_adult <- spheres3d(preds_rostrum$odont4, radius=.001, color = col_modules[rostrum])
rgl.snapshot(filename = "Output/odont_rostrum_adult1.png") 
rgl.snapshot(filename = "Output/odont_rostrum_adult.png") 
clear3d()

###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_groups_rostrum_length <- trajectory_groups_rostrum_MD[["pairwise.tables"]][["D"]]
trajectory_groups_rostrum_length_p <- trajectory_groups_rostrum_MD[["pairwise.tables"]][["P"]]
trajectory_groups_rostrum_direction <- trajectory_groups_rostrum_TC[["pairwise.tables"]][["angle"]]
trajectory_groups_rostrum_direction_p <- trajectory_groups_rostrum_TC[["pairwise.tables"]][["P"]]
trajectory_groups_rostrum_shape <- trajectory_groups_rostrum_SD[["pairwise.tables"]][["D"]]
trajectory_groups_rostrum_shape_rostrum_p <- trajectory_groups_rostrum_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_groups_rostrum_list <- list(trajectory_groups_rostrum_length, trajectory_groups_rostrum_length_p, trajectory_groups_rostrum_direction, trajectory_groups_rostrum_direction_p, 
                                       trajectory_groups_rostrum_shape, trajectory_groups_rostrum_shape_rostrum_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_groups_rostrum_list[[l]]) <- levels(groups)
  colnames(trajectory_groups_rostrum_list[[l]]) <- levels(groups)
}

#Save only lower triangle for each
trajectory_groups_rostrum_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_rostrum_lower_tri_list[[l]] <- get_upper_tri(trajectory_groups_rostrum_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_groups_rostrum_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_rostrum_melt[[l]] <- melt(trajectory_groups_rostrum_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_groups_rostrum_length_melt <- data.frame(trajectory_groups_rostrum_melt[[1]], p = trajectory_groups_rostrum_melt[[2]][[3]])
trajectory_groups_rostrum_direction_melt <- data.frame(trajectory_groups_rostrum_melt[[3]], p = trajectory_groups_rostrum_melt[[4]][[3]])
trajectory_groups_rostrum_shape_rostrum_melt <- data.frame(trajectory_groups_rostrum_melt[[5]], p = trajectory_groups_rostrum_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_groups_rostrum_length_melt <- trajectory_groups_rostrum_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                                          value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_rostrum_direction_melt <- trajectory_groups_rostrum_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_rostrum_shape_rostrum_melt <- trajectory_groups_rostrum_shape_rostrum_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_groups_rostrum_length_melt$p_if_sig))

all(is.na(trajectory_groups_rostrum_direction_melt$p_if_sig))

all(is.na(trajectory_groups_rostrum_shape_rostrum_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_groups_rostrum_direction_heatmap_ggplot <- ggplot(data = trajectory_groups_rostrum_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
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
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_groups_rostrum_direction_heatmap_ggplot

#Nice heatmap plot for each variable
trajectory_groups_rostrum_shape_heatmap_ggplot <- ggplot(data = trajectory_groups_rostrum_shape_rostrum_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low =  mypalette_seq_modules[9], high =  mypalette_seq_modules[2], mid =  mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_purple[1]) + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Shape difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 14, vjust = 0.8),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position = "none")
trajectory_groups_rostrum_shape_heatmap_ggplot

plotR4<-ggarrange(trajectory_groups_rostrum_direction_heatmap_ggplot, trajectory_groups_rostrum_shape_heatmap_ggplot,
                  ncol = 2, nrow = 1, common.legend = F)
plotR4<-annotate_figure(plotR4, top = text_grob("Rostrum", face = "bold", size = 17))
plotR4

##Braincase ----

##By group and category ----
#Need multiple observations per group/category and all have to be represented (all categories in each group)
#Check in ch. 3

#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups_braincase (Mysticeti, Odontoceti)
fit_shape_braincase_group_category <- procD.lm(coords ~ group * category, iter = 999, data = gdf_braincase, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_braincase_group_category)

#Save results to file
sink("Output/fit_shape_braincase_group_category.txt")
summary(fit_shape_braincase_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_groups_braincase <- trajectory.analysis(fit_shape_braincase_group_category, groups = gdf_braincase$group, traj.pts = gdf_braincase$category, 
                                                   pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_groups_braincase_MD <- summary(trajectory_groups_braincase, show.trajectories = TRUE, attribute = "MD") 
trajectory_groups_braincase_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_groups_braincase_TC <- summary(trajectory_groups_braincase, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_groups_braincase_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_groups_braincase_SD <- summary(trajectory_groups_braincase, show.trajectories = TRUE, attribute = "SD") 
trajectory_groups_braincase_SD 

#Save results to file
sink("Output/trajectory_groups_braincase.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_groups_braincase_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_groups_braincase_TC
print("Shape differences between trajectory vectors - shape")
trajectory_groups_braincase_SD 

print("Summary tables")
trajectory_groups_braincase_MD$summary.table
trajectory_groups_braincase_TC$summary.table
trajectory_groups_braincase_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_groups_braincase_plot <- plot(trajectory_groups_braincase, main = "Trajectories of growth by group braincase",  pch = shapes, #title and type of point to be used
                                         col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_braincase
add.trajectories(trajectory_groups_braincase_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= 0.15, y = -0.2, legend = levels(groups), 
       pch =  shapes, pt.bg = 1, cex = 1)


##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_groups_braincase_pcscores <- as_tibble(trajectory_groups_braincase_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_groups_braincase_pcscores <- trajectory_groups_braincase_pcscores %>% mutate(specimens = gdf_braincase$Id, category = gdf_braincase$category,  
                                                                                        group = gdf_braincase$group, size = gdf_braincase$size)
glimpse(trajectory_groups_braincase_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_groups_braincase_pcscores_means <- trajectory_groups_braincase_pcscores %>% group_by(group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_groups_braincase_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_groups_braincase_pcscores_means <- trajectory_groups_braincase_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_groups_braincase_pcscores_means)

#Nice plot
trajectory_groups_braincase_ggplot <- ggplot(trajectory_groups_braincase_pcscores, aes(x = PC1, y = PC2, shape = group, group = group))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_groups_braincase_pcscores_means, aes(x = x, y = y, fill = group, shape = group, alpha = category, group = group), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_groups_braincase_pcscores_means, aes(x = x, y = y, colour = group, group = group, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Group", labels = levels(groups), #copy from as.factor(groups_braincase)
                      values = mypalette_groups, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_groups_braincase[["pca"]][["sdev"]][1]^2/sum(trajectory_groups_braincase[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_groups_braincase[["pca"]][["sdev"]][2]^2/sum(trajectory_groups_braincase[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  ggtitle("Braincase")+
  guides(shape = "none", colour = "none", fill = "none",
         alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_groups_braincase_ggplot

#Add silhouettes groups_braincase
trajectory_groups_braincase_ggplot <-   
  trajectory_groups_braincase_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.25, y = -0.15, ysize = 0.03, color = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 0.15, y = 0.15, ysize = 0.025, color = mypalette_groups[2])
#Annotate adult and early fetus categories
trajectory_groups_braincase_ggplot <-
  trajectory_groups_braincase_ggplot   + 
  annotate("text", x = -0.04, y = -0.22, label = "adult", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = 0.21, y = -0.04, label = "adult", color = mypalette_groups[2], fontface = 3, size = 6)+
  annotate("text", x = -0.22, y = 0.01, label = "early\nfetus", color = mypalette_groups[1], fontface = 3, size = 6)+
  annotate("text", x = -0.11, y = 0.06, label = "early\nfetus", color = mypalette_groups[2], fontface = 3, size = 6)
trajectory_groups_braincase_ggplot


ggarrange(trajectory_groups_rostrum_ggplot, trajectory_groups_braincase_ggplot, ncol =2 , nrow =1, common.legend = T, legend = "bottom")

#Save mean shapes for each group and stage
#PC points matrix from analysis
traj_pcs_braincase <- trajectory_groups_braincase_plot$pc.points[,1:2]
traj_pcs_braincase_means <- trajectory_groups_braincase_pcscores_means[,3:4]

#User-picked spots can be anything, but it in this case evenly-spaced PCA coordinates
preds_braincase <- shape.predictor(gdf_braincase$coords, x= traj_pcs_braincase, Intercept = FALSE,
                                   myst1 = as.numeric(as.vector(traj_pcs_braincase_means[1,])),
                                   myst2 = as.numeric(as.vector(traj_pcs_braincase_means[2,])),
                                   myst3 = as.numeric(as.vector(traj_pcs_braincase_means[3,])),
                                   myst4 = as.numeric(as.vector(traj_pcs_braincase_means[4,])),
                                   odont1 = as.numeric(as.vector(traj_pcs_braincase_means[5,])),
                                   odont2 = as.numeric(as.vector(traj_pcs_braincase_means[6,])),
                                   odont3 = as.numeric(as.vector(traj_pcs_braincase_means[7,])),
                                   odont4 = as.numeric(as.vector(traj_pcs_braincase_means[8,])))

#Save shapes points
myst_braincase_early <- spheres3d(preds_braincase$myst1, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/myst_braincase_early.png") 
rgl.snapshot(filename = "Output/myst_braincase_early1.png") 
clear3d()
myst_braincase_latenew <- spheres3d(preds_braincase$myst2, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/myst_braincase_latenew1.png") 
rgl.snapshot(filename = "Output/myst_braincase_latenew.png") 
clear3d()
myst_braincase_imm <- spheres3d(preds_braincase$myst3, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/myst_braincase_imm.png") 
rgl.snapshot(filename = "Output/myst_braincase_imm1.png") 
clear3d()
myst_braincase_adult <- spheres3d(preds_braincase$myst4, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/myst_braincase_adult1.png") 
rgl.snapshot(filename = "Output/myst_braincase_adult.png") 
clear3d()

odont_braincase_early <- spheres3d(preds_braincase$odont1, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/odont_braincase_early.png") 
rgl.snapshot(filename = "Output/odont_braincase_early1.png") 
clear3d()
odont_braincase_latenew <- spheres3d(preds_braincase$odont2, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/odont_braincase_latenew1.png") 
rgl.snapshot(filename = "Output/odont_braincase_latenew.png") 
clear3d()
odont_braincase_imm <- spheres3d(preds_braincase$odont3, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/odont_braincase_imm.png") 
rgl.snapshot(filename = "Output/odont_braincase_imm1.png") 
clear3d()
odont_braincase_adult <- spheres3d(preds_braincase$odont4, radius=.004, color = col_modules[braincase])
rgl.snapshot(filename = "Output/odont_braincase_adult1.png") 
rgl.snapshot(filename = "Output/odont_braincase_adult.png") 
clear3d()

###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_groups_braincase_length <- trajectory_groups_braincase_MD[["pairwise.tables"]][["D"]]
trajectory_groups_braincase_length_p <- trajectory_groups_braincase_MD[["pairwise.tables"]][["P"]]
trajectory_groups_braincase_direction <- trajectory_groups_braincase_TC[["pairwise.tables"]][["angle"]]
trajectory_groups_braincase_direction_p <- trajectory_groups_braincase_TC[["pairwise.tables"]][["P"]]
trajectory_groups_braincase_shape <- trajectory_groups_braincase_SD[["pairwise.tables"]][["D"]]
trajectory_groups_braincase_shape_braincase_p <- trajectory_groups_braincase_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_groups_braincase_list <- list(trajectory_groups_braincase_length, trajectory_groups_braincase_length_p, trajectory_groups_braincase_direction, trajectory_groups_braincase_direction_p, 
                                         trajectory_groups_braincase_shape, trajectory_groups_braincase_shape_braincase_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_groups_braincase_list[[l]]) <- levels(groups)
  colnames(trajectory_groups_braincase_list[[l]]) <- levels(groups)
}

#Save only lower triangle for each
trajectory_groups_braincase_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_braincase_lower_tri_list[[l]] <- get_upper_tri(trajectory_groups_braincase_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_groups_braincase_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_braincase_melt[[l]] <- melt(trajectory_groups_braincase_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_groups_braincase_length_melt <- data.frame(trajectory_groups_braincase_melt[[1]], p = trajectory_groups_braincase_melt[[2]][[3]])
trajectory_groups_braincase_direction_melt <- data.frame(trajectory_groups_braincase_melt[[3]], p = trajectory_groups_braincase_melt[[4]][[3]])
trajectory_groups_braincase_shape_braincase_melt <- data.frame(trajectory_groups_braincase_melt[[5]], p = trajectory_groups_braincase_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_groups_braincase_length_melt <- trajectory_groups_braincase_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                                              value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_braincase_direction_melt <- trajectory_groups_braincase_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                    p_if_sig = ifelse(sig_p, p, NA),
                                                                                                    value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_braincase_shape_braincase_melt <- trajectory_groups_braincase_shape_braincase_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_groups_braincase_length_melt$p_if_sig))

all(is.na(trajectory_groups_braincase_direction_melt$p_if_sig))

all(is.na(trajectory_groups_braincase_shape_braincase_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_groups_braincase_direction_heatmap_ggplot <- ggplot(data = trajectory_groups_braincase_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low =  mypalette_seq_modules[9], high =  mypalette_seq_modules[2], mid =  mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_purple[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle("Direction difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 14, vjust = 0.8),
        axis.text.y =  element_text(size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_groups_braincase_direction_heatmap_ggplot

plotB4<-annotate_figure(trajectory_groups_braincase_direction_heatmap_ggplot, top = text_grob("Braincase", face = "bold", size = 17))
plotB4

ggarrange(plotR4, plotB4, ncol = 2, nrow = 1, widths = c(2,1),common.legend = F)

#TRAJECTORY ANALYSIS ROSTRUM AND BRAINCASE COMPARE ----
#Mysticeti and Odontoceti compared for each module between each other already
#Now compare modules in each group (all significant differences in modules between the 2 groups)
#Use pcscores since they are the same length, fit same result as coords

#First divide the two groups
pcscores_R_B_mysticeti <- pcscores_R_B %>% filter(group == "mysticeti")
pcscores_R_B_odontoceti <-  pcscores_R_B %>% filter(group == "odontoceti")

##Mysticeti ----
#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 3 modules (R, B, W)
fit_shape_modules_2_myst_group_category <- procD.lm(as.matrix(pcscores_R_B_mysticeti[,c(1:194)]) ~ pcscores_R_B_mysticeti$module * pcscores_R_B_mysticeti$category, iter = 999, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_modules_2_myst_group_category)

#Save results to file
sink("Output/fit_shape_modules_2_myst_group_category.txt")
summary(fit_shape_modules_2_myst_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_groups_modules_2_myst <- trajectory.analysis(fit_shape_modules_2_myst_group_category, groups = pcscores_R_B_mysticeti$module, traj.pts = pcscores_R_B_mysticeti$category, 
                                                        pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_groups_modules_2_myst_MD <- summary(trajectory_groups_modules_2_myst, show.trajectories = TRUE, attribute = "MD") 
trajectory_groups_modules_2_myst_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_groups_modules_2_myst_TC <- summary(trajectory_groups_modules_2_myst, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_groups_modules_2_myst_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_groups_modules_2_myst_SD <- summary(trajectory_groups_modules_2_myst, show.trajectories = TRUE, attribute = "SD") 
trajectory_groups_modules_2_myst_SD 

#Save results to file
sink("Output/trajectory_groups_modules_2_mysticeti.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_groups_modules_2_myst_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_groups_modules_2_myst_TC
print("Shape differences between trajectory vectors - shape")
trajectory_groups_modules_2_myst_SD 

print("Summary tables")
trajectory_groups_modules_2_myst_MD$summary.table
trajectory_groups_modules_2_myst_TC$summary.table
trajectory_groups_modules_2_myst_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_groups_modules_2_myst_plot <- plot(trajectory_groups_modules_2_myst, main = "Trajectories of growth by modules",  pch = c(22,24), #title and type of point to be used
                                              col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_rostrum
add.trajectories(trajectory_groups_modules_2_myst_plot, 
                 traj.pch = c(22,24), traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= 0.15, y = 0.25, legend = str_to_sentence(levels(as.factor(modules_2_list))), 
       pch =  c(22,24), pt.bg = 1, cex = 1)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_groups_modules_2_myst_pcscores <- as_tibble(trajectory_groups_modules_2_myst_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_groups_modules_2_myst_pcscores <- trajectory_groups_modules_2_myst_pcscores %>% mutate(specimens = pcscores_R_B_mysticeti$specimens, category = pcscores_R_B_mysticeti$category,  
                                                                                                  module = pcscores_R_B_mysticeti$module, size = pcscores_R_B_mysticeti$size)
glimpse(trajectory_groups_modules_2_myst_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_groups_modules_2_myst_pcscores_means <- trajectory_groups_modules_2_myst_pcscores %>% group_by(module, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_groups_modules_2_myst_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_groups_modules_2_myst_pcscores_means <- trajectory_groups_modules_2_myst_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_groups_modules_2_myst_pcscores_means)

#Factor modules for plot legend 
modules_2_list_factor <- str_to_sentence(levels(as.factor(modules_2_list)))

#Palette colors and shape modules
mypalette_modules_2 <- c(mypalette_paired[2],mypalette_paired[5])
shapes_modules_2 <- c(23,24)

#Nice plot
trajectory_groups_modules_2_myst_ggplot <- ggplot(trajectory_groups_modules_2_myst_pcscores, aes(x = PC1, y = PC2, shape = module, group = module))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_groups_modules_2_myst_pcscores_means, aes(x = x, y = y, fill = module, shape = module, alpha = category, group = module), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_groups_modules_2_myst_pcscores_means, aes(x = x, y = y, colour = module, group = module, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Modules", labels = modules_2_list_factor, values = shapes_modules_2)+
  scale_colour_manual(name = "Modules", labels = modules_2_list_factor, #copy from as.factor(groups_rostrum)
                      values = mypalette_modules_2, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  scale_y_reverse()+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_groups_modules_2_myst[["pca"]][["sdev"]][1]^2/sum(trajectory_groups_modules_2_myst[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_groups_modules_2_myst[["pca"]][["sdev"]][2]^2/sum(trajectory_groups_modules_2_myst[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  guides(alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_groups_modules_2_myst_ggplot

#Add silhouettes groups_rostrum
trajectory_groups_modules_2_myst_ggplot <-   
  trajectory_groups_modules_2_myst_ggplot   + 
  add_phylopic(myst, alpha = 1, x = -0.2, y = -0.12, ysize = 0.02, color = "gray30")
#Annotate adult and early fetus categories
trajectory_groups_modules_2_myst_ggplot <-
  trajectory_groups_modules_2_myst_ggplot   + 
  annotate("text", x =  0.3, y = 0.1, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = -0.16, y = 0.16, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = -0.01, y = -0.03, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)
trajectory_groups_modules_2_myst_ggplot

PC_myst <- trajectory_groups_modules_2_myst_plot[["pc.points"]][,1:2]

preds_myst_r <- shape.predictor(gdf_rostrum$coords[,,rows_mysticeti], x= NULL, Intercept = FALSE, 
                         pred1 = c(trajectory_groups_modules_2_myst_pcscores_means$x[5],trajectory_groups_modules_2_myst_pcscores_means$y[5]), 
                         pred2 = c(trajectory_groups_modules_2_myst_pcscores_means$x[8], trajectory_groups_modules_2_myst_pcscores_means$y[8]))
preds_myst_b <- shape.predictor(gdf_braincase$coords[,,rows_mysticeti], x= NULL, Intercept = FALSE, 
                                pred1 = c(trajectory_groups_modules_2_myst_pcscores_means$x[1], trajectory_groups_modules_2_myst_pcscores_means$y[1]),
                                pred2 = c(trajectory_groups_modules_2_myst_pcscores_means$x[4], trajectory_groups_modules_2_myst_pcscores_means$y[4]))

                         
###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_groups_modules_2_myst_length <- trajectory_groups_modules_2_myst_MD[["pairwise.tables"]][["D"]]
trajectory_groups_modules_2_myst_length_p <- trajectory_groups_modules_2_myst_MD[["pairwise.tables"]][["P"]]
trajectory_groups_modules_2_myst_direction <- trajectory_groups_modules_2_myst_TC[["pairwise.tables"]][["angle"]]
trajectory_groups_modules_2_myst_direction_p <- trajectory_groups_modules_2_myst_TC[["pairwise.tables"]][["P"]]
trajectory_groups_modules_2_myst_shape <- trajectory_groups_modules_2_myst_SD[["pairwise.tables"]][["D"]]
trajectory_groups_modules_2_myst_shape_p <- trajectory_groups_modules_2_myst_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_groups_modules_2_myst_list <- list(trajectory_groups_modules_2_myst_length, trajectory_groups_modules_2_myst_length_p, trajectory_groups_modules_2_myst_direction, trajectory_groups_modules_2_myst_direction_p, 
                                              trajectory_groups_modules_2_myst_shape, trajectory_groups_modules_2_myst_shape_p)

#Make vector for labels
modules_2_list_factor <-  str_to_title(modules_2_list)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_groups_modules_2_myst_list[[l]]) <- modules_2_list_factor
  colnames(trajectory_groups_modules_2_myst_list[[l]]) <- modules_2_list_factor
}

#Save only lower triangle for each
trajectory_groups_modules_2_myst_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_modules_2_myst_lower_tri_list[[l]] <- get_upper_tri(trajectory_groups_modules_2_myst_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_groups_modules_2_myst_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_modules_2_myst_melt[[l]] <- melt(trajectory_groups_modules_2_myst_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_groups_modules_2_myst_length_melt <- data.frame(trajectory_groups_modules_2_myst_melt[[1]], p = trajectory_groups_modules_2_myst_melt[[2]][[3]])
trajectory_groups_modules_2_myst_direction_melt <- data.frame(trajectory_groups_modules_2_myst_melt[[3]], p = trajectory_groups_modules_2_myst_melt[[4]][[3]])
trajectory_groups_modules_2_myst_shape_modules_2_myst_melt <- data.frame(trajectory_groups_modules_2_myst_melt[[5]], p = trajectory_groups_modules_2_myst_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_groups_modules_2_myst_length_melt <- trajectory_groups_modules_2_myst_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_modules_2_myst_direction_melt <- trajectory_groups_modules_2_myst_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                                                              value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))
trajectory_groups_modules_2_myst_shape_modules_2_myst_melt <- trajectory_groups_modules_2_myst_shape_modules_2_myst_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                                    p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                                    value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_groups_modules_2_myst_length_melt$p_if_sig))

all(is.na(trajectory_groups_modules_2_myst_direction_melt$p_if_sig))

all(is.na(trajectory_groups_modules_2_myst_shape_modules_2_myst_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_groups_modules_2_myst_direction_heatmap_ggplot <- ggplot(data = trajectory_groups_modules_2_myst_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high =  mypalette_seq_groups[2], mid =  mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Direction difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_groups_modules_2_myst_direction_heatmap_ggplot

plotTmys<-annotate_figure(trajectory_groups_modules_2_myst_direction_heatmap_ggplot, 
                          top = text_grob("Mysticeti", face = "bold", size = 15, just = c(0.5,10)))
plotTmys

##Odontoceti ----
#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 3 modules (R, B, W)
fit_shape_modules_2_odont_group_category <- procD.lm(as.matrix(pcscores_R_B_odontoceti[,c(1:194)]) ~ pcscores_R_B_odontoceti$module * pcscores_R_B_odontoceti$category, iter = 999, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_modules_2_odont_group_category)

#Save results to file
sink("Output/fit_shape_modules_2_odont_group_category.txt")
summary(fit_shape_modules_2_odont_group_category)
sink() 

#Use fit to calculate trajectories
trajectory_groups_modules_2_odont <- trajectory.analysis(fit_shape_modules_2_odont_group_category, groups = pcscores_R_B_odontoceti$module, traj.pts = pcscores_R_B_odontoceti$category, 
                                                         pca = TRUE, print.progress = F) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_groups_modules_2_odont_MD <- summary(trajectory_groups_modules_2_odont, show.trajectories = TRUE, attribute = "MD") 
trajectory_groups_modules_2_odont_MD

#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_groups_modules_2_odont_TC <- summary(trajectory_groups_modules_2_odont, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_groups_modules_2_odont_TC 

#Trajectory shape differences - are trajectories different in shape?
trajectory_groups_modules_2_odont_SD <- summary(trajectory_groups_modules_2_odont, show.trajectories = TRUE, attribute = "SD") 
trajectory_groups_modules_2_odont_SD 

#Save results to file
sink("Output/trajectory_groups_modules_2_odontoceti.txt")
print("Magnitude difference (absolute difference between path distances) - length")
trajectory_groups_modules_2_odont_MD 
print("Correlations (angles) between trajectories - direction")
trajectory_groups_modules_2_odont_TC
print("Shape differences between trajectory vectors - shape")
trajectory_groups_modules_2_odont_SD 

print("Summary tables")
trajectory_groups_modules_2_odont_MD$summary.table
trajectory_groups_modules_2_odont_TC$summary.table
trajectory_groups_modules_2_odont_SD$summary.table
sink() 

#Plot results - PCA of fitted values
trajectory_groups_modules_2_odont_plot <- plot(trajectory_groups_modules_2_odont, main = "Trajectories of growth by modules",  pch = c(22,24), #title and type of point to be used
                                               col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups_rostrum
add.trajectories(trajectory_groups_modules_2_odont_plot, 
                 traj.pch = c(22,24), traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= -0.5, y = 0.3, legend = str_to_sentence(levels(as.factor(modules_2_list))), 
       pch =  c(22,24), pt.bg = 1, cex = 1)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_groups_modules_2_odont_pcscores <- as_tibble(trajectory_groups_modules_2_odont_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_groups_modules_2_odont_pcscores <- trajectory_groups_modules_2_odont_pcscores %>% mutate(specimens = pcscores_R_B_odontoceti$specimens, category = pcscores_R_B_odontoceti$category,  
                                                                                                    module = pcscores_R_B_odontoceti$module, size = pcscores_R_B_odontoceti$size)
glimpse(trajectory_groups_modules_2_odont_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
trajectory_groups_modules_2_odont_pcscores_means <- trajectory_groups_modules_2_odont_pcscores %>% group_by(module, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_groups_modules_2_odont_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_groups_modules_2_odont_pcscores_means <- trajectory_groups_modules_2_odont_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_groups_modules_2_odont_pcscores_means)

#Nice plot
trajectory_groups_modules_2_odont_ggplot <- ggplot(trajectory_groups_modules_2_odont_pcscores, aes(x = PC1, y = PC2, shape = module, group = module))+
  geom_point(aes(alpha = category), size = 2, colour = "grey10", fill = "grey10", show.legend = F)+
  geom_point(data = trajectory_groups_modules_2_odont_pcscores_means, aes(x = x, y = y, fill = module, shape = module, alpha = category, group = module), colour = "grey10",
             size = 6, inherit.aes = F)+
  geom_path(data = trajectory_groups_modules_2_odont_pcscores_means, aes(x = x, y = y, colour = module, group = module, linewidth = category), inherit.aes = F,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.02, "npc"), ends = "last", type = "open"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.4, 0.6, 0.8))+            #legend and color adjustments
  scale_shape_manual(name = "Modules", labels = modules_2_list_factor, values = shapes_modules_2)+
  scale_colour_manual(name = "Modules", labels = modules_2_list_factor, #copy from as.factor(groups_rostrum)
                      values = mypalette_modules_2, aesthetics = c("colour", "fill"))+
  scale_linewidth_manual(values = c(0.8, 1.2, 1.5, 1.8))+
  scale_y_reverse()+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_groups_modules_2_odont[["pca"]][["sdev"]][1]^2/sum(trajectory_groups_modules_2_odont[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_groups_modules_2_odont[["pca"]][["sdev"]][2]^2/sum(trajectory_groups_modules_2_odont[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  guides(alpha = guide_legend(override.aes = list(fill = NA, colour = "grey10")))+
  theme(plot.title = element_text(size = 16, face = 3, hjust = 0.5),legend.text = element_text(size = 11), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", 
        legend.direction = "horizontal", legend.justification = c(0,0))
trajectory_groups_modules_2_odont_ggplot

#Add silhouettes groups_rostrum
trajectory_groups_modules_2_odont_ggplot <-   
  trajectory_groups_modules_2_odont_ggplot   + 
  add_phylopic(odont, alpha = 1, x = 0.15, y = -0.28, ysize = 0.03, color = "gray50")
#Annotate adult and early fetus categories
trajectory_groups_modules_2_odont_ggplot <-
  trajectory_groups_modules_2_odont_ggplot   + 
  annotate("text", x = 0.1, y = 0.01, label = "adult", color = "grey10", alpha = 0.8, fontface = 3, size = 6)+
  annotate("text", x = -0.22, y = 0.05, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)+
  annotate("text", x = -0.14, y = -0.21, label = "early\nfetus", color = "grey10", alpha = 0.4, fontface = 3, size = 6)
trajectory_groups_modules_2_odont_ggplot

ggarrange(trajectory_groups_modules_2_myst_ggplot, trajectory_groups_modules_2_odont_ggplot,
          ncol = 2, nrow = 1, common.legend = T, legend = "bottom")

###Heatmaps for pairwise comparison trajectories ----

#Save p-values as object
trajectory_groups_modules_2_odont_length <- trajectory_groups_modules_2_odont_MD[["pairwise.tables"]][["D"]]
trajectory_groups_modules_2_odont_length_p <- trajectory_groups_modules_2_odont_MD[["pairwise.tables"]][["P"]]
trajectory_groups_modules_2_odont_direction <- trajectory_groups_modules_2_odont_TC[["pairwise.tables"]][["angle"]]
trajectory_groups_modules_2_odont_direction_p <- trajectory_groups_modules_2_odont_TC[["pairwise.tables"]][["P"]]
trajectory_groups_modules_2_odont_shape <- trajectory_groups_modules_2_odont_SD[["pairwise.tables"]][["D"]]
trajectory_groups_modules_2_odont_shape_p <- trajectory_groups_modules_2_odont_SD[["pairwise.tables"]][["P"]]

#Make list to change tables faster
trajectory_groups_modules_2_odont_list <- list(trajectory_groups_modules_2_odont_length, trajectory_groups_modules_2_odont_length_p, trajectory_groups_modules_2_odont_direction, trajectory_groups_modules_2_odont_direction_p, 
                                               trajectory_groups_modules_2_odont_shape, trajectory_groups_modules_2_odont_shape_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(trajectory_groups_modules_2_odont_list[[l]]) <- modules_2_list_factor
  colnames(trajectory_groups_modules_2_odont_list[[l]]) <- modules_2_list_factor
}

#Save only lower triangle for each
trajectory_groups_modules_2_odont_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_modules_2_odont_lower_tri_list[[l]] <- get_upper_tri(trajectory_groups_modules_2_odont_list[[l]])
}

#Melt to make table in the format needed for heatmap
trajectory_groups_modules_2_odont_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  trajectory_groups_modules_2_odont_melt[[l]] <- melt(trajectory_groups_modules_2_odont_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
trajectory_groups_modules_2_odont_length_melt <- data.frame(trajectory_groups_modules_2_odont_melt[[1]], p = trajectory_groups_modules_2_odont_melt[[2]][[3]])
trajectory_groups_modules_2_odont_direction_melt <- data.frame(trajectory_groups_modules_2_odont_melt[[3]], p = trajectory_groups_modules_2_odont_melt[[4]][[3]])
trajectory_groups_modules_2_odont_shape_melt <- data.frame(trajectory_groups_modules_2_odont_melt[[5]], p = trajectory_groups_modules_2_odont_melt[[6]][[3]])

#Create columns where only significant values are shown
trajectory_groups_modules_2_odont_length_melt <- trajectory_groups_modules_2_odont_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                                                          value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))
trajectory_groups_modules_2_odont_direction_melt <- trajectory_groups_modules_2_odont_direction_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))
trajectory_groups_modules_2_odont_shape_melt <- trajectory_groups_modules_2_odont_shape_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(trajectory_groups_modules_2_odont_length_melt$p_if_sig))

all(is.na(trajectory_groups_modules_2_odont_direction_melt$p_if_sig))

all(is.na(trajectory_groups_modules_2_odont_shape_melt$p_if_sig))

#Nice heatmap plot for each variable
trajectory_groups_modules_2_odont_length_heatmap_ggplot <- ggplot(data = trajectory_groups_modules_2_odont_length_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high =  mypalette_seq_groups[2], mid =  mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Length difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
trajectory_groups_modules_2_odont_length_heatmap_ggplot

trajectory_groups_modules_2_odont_direction_heatmap_ggplot <- ggplot(data = trajectory_groups_modules_2_odont_direction_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high =  mypalette_seq_groups[2], mid =  mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Direction difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = "none")
trajectory_groups_modules_2_odont_direction_heatmap_ggplot

trajectory_groups_modules_2_odont_shape_heatmap_ggplot <- ggplot(data = trajectory_groups_modules_2_odont_shape_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high =  mypalette_seq_groups[2], mid =  mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1],  name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Shape difference trajectories")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = "none")
trajectory_groups_modules_2_odont_shape_heatmap_ggplot

plotTod<-ggarrange(trajectory_groups_modules_2_odont_length_heatmap_ggplot,
                   trajectory_groups_modules_2_odont_direction_heatmap_ggplot, trajectory_groups_modules_2_odont_shape_heatmap_ggplot,
                   ncol = 3, nrow = 1, common.legend = F)
plotTod<-annotate_figure(plotTod, top = text_grob("Odontoceti", face = "bold", size = 15, just = c(0.5, 10)))
plotTod

ggarrange(plotTmys, plotTod, ncol = 2, nrow = 1, widths = c(1,3), common.legend = F)

###### 

#Next - ch. 8 - Allometry analyses