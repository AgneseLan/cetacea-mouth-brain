#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH. 8a - Allometry analyses rostrum vs braincase, mcp testing  - phylogenetically transformed components

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
#devtools::install_github("kassambara/easyGgplot2")
#remotes::install_github("r-lib/rray")

# sink("Output/session_info.txt")
# sessioninfo::session_info()
# sink()

#ALLOMETRY ANALYSIS - ROSTRUM AND BRAINCASE SEPARATE ----
##Test different allometry between category in each group in rostrum and braincase ----
##For both use overall logCS - allows even comparison and shows progressive growth of skull relative to size

###Rostrum ----

##Regression shape on logCS size rostrum
allometry_phylo_rostrum <- procD.lm(as.matrix(phylo_pcscores_rostrum_df[,c(1:194)]) ~ phylo_pcscores_rostrum_df$size, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_phylo_rostrum)

#Make interaction factor between group and category
phylo_pcscores_rostrum_df$grp_cat <- interaction(phylo_pcscores_rostrum_df$group, phylo_pcscores_rostrum_df$category)

#Allometry by category by group models
allometry_phylo_rostrum_grp_cat_comb <-  procD.lm(as.matrix(phylo_pcscores_rostrum_df[,c(1:194)]) ~ phylo_pcscores_rostrum_df$size + phylo_pcscores_rostrum_df$grp_cat, iter=999, print.progress = TRUE) 
allometry_phylo_rostrum_grp_cat_int <-  procD.lm(as.matrix(phylo_pcscores_rostrum_df[,c(1:194)]) ~ phylo_pcscores_rostrum_df$size * phylo_pcscores_rostrum_df$grp_cat, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_phylo_rostrum_grp_cat_comb)
summary(allometry_phylo_rostrum_grp_cat_int) 

#Save results of significant regression to file
sink("Output/8a-Allometry phylo/allometry_phylo_shape_size_grp_cat_rostrum.txt")
print("Null")
summary(allometry_phylo_rostrum)

print("Combination +")
summary(allometry_phylo_rostrum_grp_cat_comb) 

print("Interaction *")
summary(allometry_phylo_rostrum_grp_cat_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_phylo_models_grp_cat_rostrum <- anova(allometry_phylo_rostrum, allometry_phylo_rostrum_grp_cat_comb, allometry_phylo_rostrum_grp_cat_int)
anova_allometry_phylo_models_grp_cat_rostrum
#If models equivalent, compare between each other
anova_allometry_phylo_models_grp_cat_rostrum1 <- anova(allometry_phylo_rostrum_grp_cat_comb, allometry_phylo_rostrum_grp_cat_int)
anova_allometry_phylo_models_grp_cat_rostrum1

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the allometry trajectory on top of difference in intercept (comb model)
pairwise_allometry_phylo_rostrum_grp_cat <- pairwise(allometry_phylo_rostrum_grp_cat_int, fit.null = allometry_phylo_rostrum_grp_cat_comb,
                                               groups = phylo_pcscores_rostrum_df$grp_cat, 
                                               covariate =  phylo_pcscores_rostrum_df$size, print.progress = FALSE) 
pairwise_allometry_phylo_rostrum_grp_cat

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_phylo_rostrum_grp_cat_dis <- summary(pairwise_allometry_phylo_rostrum_grp_cat, confidence = 0.95, test.type = "dist") 
pairwise_allometry_phylo_rostrum_grp_cat_dis

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_phylo_rostrum_grp_cat_VC <- summary(pairwise_allometry_phylo_rostrum_grp_cat, confidence = 0.95, test.type = "VC",
                                                 angle.type = "deg") 
pairwise_allometry_phylo_rostrum_grp_cat_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_phylo_rostrum_grp_cat_DL <-summary(pairwise_allometry_phylo_rostrum_grp_cat, confidence = 0.95, test.type = "DL") 
pairwise_allometry_phylo_rostrum_grp_cat_DL 

#Compare the dispersion around group slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
pairwise_allometry_phylo_rostrum_grp_cat_var <-summary(pairwise_allometry_phylo_rostrum_grp_cat, confidence = 0.95, test.type = "var")
pairwise_allometry_phylo_rostrum_grp_cat_var

#Save results to file
sink("Output/8a-Allometry phylo/pairwise_allometry_phylo_rostrum_grp_cat.txt")
print("ANOVA models")
print(anova_allometry_phylo_models_grp_cat_rostrum)

print("1-Pairwise absolute distances slopes")
pairwise_allometry_phylo_rostrum_grp_cat_dis 

print("2-Distance between angles (slope directions)")
pairwise_allometry_phylo_rostrum_grp_cat_VC

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
pairwise_allometry_phylo_rostrum_grp_cat_DL

print("4-Difference in dispersion around mean slope")
pairwise_allometry_phylo_rostrum_grp_cat_var
sink()

#Regression score of shape vs logCS and comb or int (best model)- regression method with "RegScore" plotting
allometry_phylo_rostrum_grp_cat_plot_regscore <- plot(allometry_phylo_rostrum_grp_cat_int, type = "regression",predictor = gdf$size, reg.type = "RegScore",
                                                main = "Shape vs logCS * grp_cat phylo corrected",xlab = "logCS", pch = 21, col = mypalette_paired[5], 
                                                bg = mypalette_paired[5], cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf$size, y = allometry_phylo_rostrum_grp_cat_plot_regscore$RegScore, labels = Ids,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

####Plot allometry rostrum by category by group ----
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_phylo_rostrum_grp_cat_plot <- data.frame(logCS = allometry_phylo_rostrum_grp_cat_plot_regscore[["plot_args"]][["x"]], 
                                             RegScores = allometry_phylo_rostrum_grp_cat_plot_regscore[["plot_args"]][["y"]])
#Convert data frame to tibble
allometry_phylo_rostrum_grp_cat_plot <- as_tibble(allometry_phylo_rostrum_grp_cat_plot)
#Add labels and other attributes to tibble as columns
allometry_phylo_rostrum_grp_cat_plot <- allometry_phylo_rostrum_grp_cat_plot %>% 
  mutate(specimens = gdf$Id, family = gdf$family, category = gdf_rostrum$category, group = gdf_rostrum$group,
         grp_cat = phylo_pcscores_rostrum_df$grp_cat)
glimpse(allometry_phylo_rostrum_grp_cat_plot)

#Plot allometry regression by category by group 
allometry_phylo_rostrum_grp_cat_ggplot <- ggplot(allometry_phylo_rostrum_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(aes(colour = category,fill = category), size = 0, alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group,colour = category,fill = category, group = grp_cat), inherit.aes = F,        
              se = F, linewidth = 1.5, alpha = 1)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  ggtitle("Rostrum")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_rostrum_grp_cat_ggplot

#Add phylopic
allometry_phylo_rostrum_grp_cat_ggplot <- 
  allometry_phylo_rostrum_grp_cat_ggplot +
  add_phylopic(myst, alpha = 1, x = 3.8, y = -0.05, ysize = 0.12, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2.7, y = -0.05, ysize = 0.1, fill = "gray50")
allometry_phylo_rostrum_grp_cat_ggplot

####Heatmaps plots for significant differences in pairwise ----

#Save p-values as object
pairwise_allometry_phylo_rostrum_grp_cat_dist <- pairwise_allometry_phylo_rostrum_grp_cat_dis[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_rostrum_grp_cat_dist_p <- pairwise_allometry_phylo_rostrum_grp_cat_dis[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_rostrum_grp_cat_angle <- pairwise_allometry_phylo_rostrum_grp_cat_VC[["pairwise.tables"]][["angle"]]
pairwise_allometry_phylo_rostrum_grp_cat_angle_p <- pairwise_allometry_phylo_rostrum_grp_cat_VC[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_rostrum_grp_cat_length <- pairwise_allometry_phylo_rostrum_grp_cat_DL[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_rostrum_grp_cat_length_p <- pairwise_allometry_phylo_rostrum_grp_cat_DL[["pairwise.tables"]][["P"]]

#Make list to change tables faster
pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_list <- list(pairwise_allometry_phylo_rostrum_grp_cat_dist, pairwise_allometry_phylo_rostrum_grp_cat_dist_p, pairwise_allometry_phylo_rostrum_grp_cat_angle, pairwise_allometry_phylo_rostrum_grp_cat_angle_p, 
                                                        pairwise_allometry_phylo_rostrum_grp_cat_length, pairwise_allometry_phylo_rostrum_grp_cat_length_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  rownames(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_list[[l]]) <- disparity_rostrum_vars
  colnames(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_list[[l]]) <- disparity_rostrum_vars
}

#Save only lower triangle for each
pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_lower_tri_list[[l]] <- get_upper_tri(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_list[[l]])
}

#Melt to make table in the format needed for heatmap
pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[l]] <- melt(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
pairwise_allometry_phylo_rostrum_grp_cat_dist_melt <- data.frame(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[1]], p = pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[2]][[3]])
pairwise_allometry_phylo_rostrum_grp_cat_angle_melt <- data.frame(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[3]], p = pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[4]][[3]])
pairwise_allometry_phylo_rostrum_grp_cat_length_melt <- data.frame(pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[5]], p = pairwise_allometry_phylo_rostrum_grp_cat_grp_cat_melt[[6]][[3]])

#Create columns where only significant values are shown
pairwise_allometry_phylo_rostrum_grp_cat_dist_melt <- pairwise_allometry_phylo_rostrum_grp_cat_dist_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))
pairwise_allometry_phylo_rostrum_grp_cat_angle_melt <- pairwise_allometry_phylo_rostrum_grp_cat_angle_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                                                          value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 1)))
pairwise_allometry_phylo_rostrum_grp_cat_length_melt <- pairwise_allometry_phylo_rostrum_grp_cat_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                            p_if_sig = ifelse(sig_p, p, NA),
                                                                                                            value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(pairwise_allometry_phylo_rostrum_grp_cat_dist_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_rostrum_grp_cat_angle_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_rostrum_grp_cat_length_melt$p_if_sig))

#Nice heatmap plot for each variable
pairwise_allometry_phylo_rostrum_grp_cat_dist_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_rostrum_grp_cat_dist_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8, hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1.2,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_rostrum_grp_cat_dist_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_rostrum_grp_cat_angle_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_rostrum_grp_cat_angle_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8, hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_phylo_rostrum_grp_cat_angle_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_rostrum_grp_cat_length_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_rostrum_grp_cat_length_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope length difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8, hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_phylo_rostrum_grp_cat_length_heatmap_ggplot

plotR3<-ggarrange(pairwise_allometry_phylo_rostrum_grp_cat_dist_heatmap_ggplot, pairwise_allometry_phylo_rostrum_grp_cat_angle_heatmap_ggplot,
                  pairwise_allometry_phylo_rostrum_grp_cat_length_heatmap_ggplot, ncol = 3, nrow = 1, common.legend = F)
plotR3<-annotate_figure(plotR3, top = text_grob("Rostrum", face = "bold", size = 17))
plotR3

###Braincase ----
##Regression shape on logCS size braincase
allometry_phylo_braincase <- procD.lm(as.matrix(phylo_pcscores_braincase_df[,c(1:194)]) ~ phylo_pcscores_braincase_df$size, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_phylo_braincase)

#Make interaction factor between group and category
phylo_pcscores_braincase_df$grp_cat <- interaction(phylo_pcscores_braincase_df$group, phylo_pcscores_braincase_df$category)

#Allometry by category by group models
allometry_phylo_braincase_grp_cat_comb <-  procD.lm(as.matrix(phylo_pcscores_braincase_df[,c(1:194)]) ~ phylo_pcscores_braincase_df$size + phylo_pcscores_braincase_df$grp_cat, iter=999, print.progress = TRUE) 
allometry_phylo_braincase_grp_cat_int <-  procD.lm(as.matrix(phylo_pcscores_braincase_df[,c(1:194)]) ~ phylo_pcscores_braincase_df$size * phylo_pcscores_braincase_df$grp_cat, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_phylo_braincase_grp_cat_comb)
summary(allometry_phylo_braincase_grp_cat_int) 

#Save results of significant regression to file
sink("Output/8a-Allometry phylo/allometry_phylo_shape_size_grp_cat_braincase.txt")
print("Null")
summary(allometry_phylo_braincase)

print("Combination +")
summary(allometry_phylo_braincase_grp_cat_comb) 

print("Interaction *")
summary(allometry_phylo_braincase_grp_cat_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_phylo_models_grp_cat_braincase <- anova(allometry_phylo_braincase, allometry_phylo_braincase_grp_cat_comb, allometry_phylo_braincase_grp_cat_int)
anova_allometry_phylo_models_grp_cat_braincase
#If models equivalent, compare between each other
anova_allometry_phylo_models_grp_cat_braincase1 <- anova(allometry_phylo_braincase_grp_cat_comb, allometry_phylo_braincase_grp_cat_int)
anova_allometry_phylo_models_grp_cat_braincase1

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the allometry trajectory on top of difference in intercept (comb model)
pairwise_allometry_phylo_braincase_grp_cat <- pairwise(allometry_phylo_braincase_grp_cat_int, fit.null = allometry_phylo_braincase_grp_cat_comb,
                                                 groups = phylo_pcscores_braincase_df$grp_cat, 
                                                 covariate =  phylo_pcscores_braincase_df$size, print.progress = FALSE) 
pairwise_allometry_phylo_braincase_grp_cat

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_phylo_braincase_grp_cat_dis <- summary(pairwise_allometry_phylo_braincase_grp_cat, confidence = 0.95, test.type = "dist") 
pairwise_allometry_phylo_braincase_grp_cat_dis

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_phylo_braincase_grp_cat_VC <- summary(pairwise_allometry_phylo_braincase_grp_cat, confidence = 0.95, test.type = "VC",
                                                   angle.type = "deg") 
pairwise_allometry_phylo_braincase_grp_cat_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_phylo_braincase_grp_cat_DL <-summary(pairwise_allometry_phylo_braincase_grp_cat, confidence = 0.95, test.type = "DL") 
pairwise_allometry_phylo_braincase_grp_cat_DL 

#Compare the dispersion around group slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
pairwise_allometry_phylo_braincase_grp_cat_var <-summary(pairwise_allometry_phylo_braincase_grp_cat, confidence = 0.95, test.type = "var")
pairwise_allometry_phylo_braincase_grp_cat_var

#Save results to file
sink("Output/8a-Allometry phylo/pairwise_allometry_phylo_braincase_grp_cat.txt")
print("ANOVA models")
print(anova_allometry_phylo_models_grp_cat_braincase)

print("1-Pairwise absolute distances slopes")
pairwise_allometry_phylo_braincase_grp_cat_dis

print("2-Distance between angles (slope directions)")
pairwise_allometry_phylo_braincase_grp_cat_VC

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
pairwise_allometry_phylo_braincase_grp_cat_DL

print("4-Difference in dispersion around mean slope")
pairwise_allometry_phylo_braincase_grp_cat_var
sink()

#Regression score of shape vs logCS and comb or int (best model)- regression method with "RegScore" plotting
allometry_phylo_braincase_grp_cat_plot_regscore <- plot(allometry_phylo_braincase_grp_cat_int, type = "regression",predictor = gdf$size, reg.type = "RegScore",
                                                  main = "Shape vs logCS * grp_cat phylo corrected",xlab = "logCS", pch = 21, col = mypalette_paired[1], 
                                                  bg = mypalette_paired[1], cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf$size, y = allometry_phylo_braincase_grp_cat_plot_regscore$RegScore, labels = Ids,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

####Plot allometry braincase by category by group ----
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_phylo_braincase_grp_cat_plot <- data.frame(logCS = allometry_phylo_braincase_grp_cat_plot_regscore[["plot_args"]][["x"]], 
                                               RegScores = allometry_phylo_braincase_grp_cat_plot_regscore[["plot_args"]][["y"]])
#Convert data frame to tibble
allometry_phylo_braincase_grp_cat_plot <- as_tibble(allometry_phylo_braincase_grp_cat_plot)
#Add labels and other attributes to tibble as columns
allometry_phylo_braincase_grp_cat_plot <- allometry_phylo_braincase_grp_cat_plot %>% 
  mutate(specimens = gdf_braincase$Id,  family = gdf_braincase$family, category = gdf_braincase$category, group = gdf_braincase$group,
         grp_cat = phylo_pcscores_braincase_df$grp_cat)
glimpse(allometry_phylo_braincase_grp_cat_plot)

#Plot allometry regression by category by group 
allometry_phylo_braincase_grp_cat_ggplot <- ggplot(allometry_phylo_braincase_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(aes(colour = category,fill = category), size = 0, alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group,colour = category,fill = category, group = grp_cat), inherit.aes = F,        
              se = F, linewidth = 1.5, alpha = 1)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  ggtitle("Braincase")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_braincase_grp_cat_ggplot

#Add phylopic
allometry_phylo_braincase_grp_cat_ggplot <- 
  allometry_phylo_braincase_grp_cat_ggplot +
  add_phylopic(myst, alpha = 1, x = 3.6, y = -0.3, ysize = 0.09, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2.5, y = -0.1, ysize = 0.085, fill = "gray50")
allometry_phylo_braincase_grp_cat_ggplot

ggarrange(allometry_phylo_rostrum_grp_cat_ggplot, allometry_phylo_braincase_grp_cat_ggplot,
          ncol = 2, nrow = 1, common.legend = T, legend = "bottom")

####Heatmaps plots for significant differences in pairwise ----

#Save p-values as object
pairwise_allometry_phylo_braincase_grp_cat_dist <- pairwise_allometry_phylo_braincase_grp_cat_dis[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_braincase_grp_cat_dist_p <- pairwise_allometry_phylo_braincase_grp_cat_dis[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_braincase_grp_cat_angle <- pairwise_allometry_phylo_braincase_grp_cat_VC[["pairwise.tables"]][["angle"]]
pairwise_allometry_phylo_braincase_grp_cat_angle_p <- pairwise_allometry_phylo_braincase_grp_cat_VC[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_braincase_grp_cat_length <- pairwise_allometry_phylo_braincase_grp_cat_DL[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_braincase_grp_cat_length_p <- pairwise_allometry_phylo_braincase_grp_cat_DL[["pairwise.tables"]][["P"]]

#Make list to change tables faster
pairwise_allometry_phylo_braincase_grp_cat_grp_cat_list <- list(pairwise_allometry_phylo_braincase_grp_cat_dist, pairwise_allometry_phylo_braincase_grp_cat_dist_p, pairwise_allometry_phylo_braincase_grp_cat_angle, pairwise_allometry_phylo_braincase_grp_cat_angle_p, 
                                                          pairwise_allometry_phylo_braincase_grp_cat_length, pairwise_allometry_phylo_braincase_grp_cat_length_p)

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_list[[l]]) <- disparity_braincase_vars
  colnames(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_list[[l]]) <- disparity_braincase_vars
}

#Save only lower triangle for each
pairwise_allometry_phylo_braincase_grp_cat_grp_cat_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  pairwise_allometry_phylo_braincase_grp_cat_grp_cat_lower_tri_list[[l]] <- get_upper_tri(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_list[[l]])
}

#Melt to make table in the format needed for heatmap
pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[l]] <- melt(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
pairwise_allometry_phylo_braincase_grp_cat_dist_melt <- data.frame(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[1]], p = pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[2]][[3]])
pairwise_allometry_phylo_braincase_grp_cat_angle_melt <- data.frame(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[3]], p = pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[4]][[3]])
pairwise_allometry_phylo_braincase_grp_cat_length_melt <- data.frame(pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[5]], p = pairwise_allometry_phylo_braincase_grp_cat_grp_cat_melt[[6]][[3]])

#Create columns where only significant values are shown
pairwise_allometry_phylo_braincase_grp_cat_dist_melt <- pairwise_allometry_phylo_braincase_grp_cat_dist_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                            p_if_sig = ifelse(sig_p, p, NA),
                                                                                                            value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))
pairwise_allometry_phylo_braincase_grp_cat_angle_melt <- pairwise_allometry_phylo_braincase_grp_cat_angle_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                                                              value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 1)))
pairwise_allometry_phylo_braincase_grp_cat_length_melt <- pairwise_allometry_phylo_braincase_grp_cat_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(pairwise_allometry_phylo_braincase_grp_cat_dist_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_braincase_grp_cat_angle_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_braincase_grp_cat_length_melt$p_if_sig))

#Nice heatmap plot for each variable
pairwise_allometry_phylo_braincase_grp_cat_dist_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_braincase_grp_cat_dist_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8, hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 9))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1.2,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_braincase_grp_cat_dist_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_braincase_grp_cat_angle_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_braincase_grp_cat_angle_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8,hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_phylo_braincase_grp_cat_angle_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_braincase_grp_cat_length_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_braincase_grp_cat_length_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope length difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13, vjust = 0.8,hjust = 0.7),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_phylo_braincase_grp_cat_length_heatmap_ggplot

plotB3<-ggarrange(pairwise_allometry_phylo_braincase_grp_cat_dist_heatmap_ggplot, pairwise_allometry_phylo_braincase_grp_cat_angle_heatmap_ggplot,pairwise_allometry_phylo_braincase_grp_cat_length_heatmap_ggplot,
                  ncol = 3, nrow = 1, common.legend = F)
plotB3<-annotate_figure(plotB3, top = text_grob("Braincase", face = "bold", size = 17))
plotB3

ggarrange(plotR3, plotB3, ncol = 1, nrow = 2,common.legend = T)

#TEST ALLOMETRY REGRESSION MULTIPLE CHANGE POINTS (Morris et al. 2021) ----
#Test if allometry should be considered separate for each growth category in rostrum and braincase
#Compare models dividing by groups, check if 3 break points = growth stages are preferred model  for both groups

#Rostrum ----

#Null model - no break points and single slope/intercept
null_model_allometry_phylo_r <- list(RegScores ~ logCS)

null_mcp_phylo_r <- mcp(null_model_allometry_phylo_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(null_mcp_phylo_r)
summary(null_mcp_phylo_r)

#Null model - no break point slope/intercept by group
model_allometry_phylo_null_group_r = list(RegScores ~ 1 + (1|group) + logCS)

mcp_phylo_null_group_r <- mcp(model_allometry_phylo_null_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_null_group_r)
summary(mcp_phylo_null_group_r)

#Model for single break point with multiple intercepts and slopes
model_allometry_phylo_1bp_r = list(RegScores ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_1bp_r <- mcp(model_allometry_phylo_1bp_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_1bp_r)
plot_pars(mcp_phylo_1bp_r, pars = "cp_1", type = "dens") + theme_bw(10)
summary(mcp_phylo_1bp_r)

#Model for 2 break points with multiple intercepts and slopes
model_allometry_phylo_2bp_r = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_2bp_r <- mcp(model_allometry_phylo_2bp_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_2bp_r)
plot_pars(mcp_phylo_2bp_r, pars = c("cp_1", "cp_2"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_2bp_r)

#Model for 3 break points with multiple intercepts and slopes
model_allometry_phylo_3bp_r = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_3bp_r <- mcp(model_allometry_phylo_3bp_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_3bp_r)
plot_pars(mcp_phylo_3bp_r, pars = c("cp_1", "cp_2", "cp_3"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_3bp_r)

#Model for 4 break points with multiple intercepts and slopes  - this corresponds to the growth stages
model_allometry_phylo_4bp_r = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS, ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_4bp_r <- mcp(model_allometry_phylo_4bp_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_4bp_r)
plot_pars(mcp_phylo_4bp_r, pars = c("cp_1", "cp_2", "cp_3", "cp_4"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_4bp_r)

#Model for single break point with variance among genera and multiple intercepts and slopes in both lines
model_allometry_phylo_1bp_group_r = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_1bp_group_r <- mcp(model_allometry_phylo_1bp_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_1bp_group_r, facet_by="group")
summary(mcp_phylo_1bp_group_r)
plot_pars(mcp_phylo_1bp_group_r, pars = "cp_1", type = "dens_overlay")

#Model for single break point with variance among genera and multiple intercepts and slopes in both lines
model_allometry_phylo_1bp_group_r = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_1bp_group_r <- mcp(model_allometry_phylo_1bp_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_1bp_group_r, facet_by="group")
summary(mcp_phylo_1bp_group_r)
plot_pars(mcp_phylo_1bp_group_r, pars = "cp_1", type = "dens_overlay")

#Model for 2 break points with variance among genera and multiple intercepts and slopes in both lines
model_allometry_phylo_2bp_group_r = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_2bp_group_r <- mcp(model_allometry_phylo_2bp_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_2bp_group_r, facet_by="group")
summary(mcp_phylo_2bp_group_r)
ranef(mcp_phylo_2bp_group_r)
plot_pars(mcp_phylo_2bp_group_r, pars = c("cp_1", "cp_2"), type = "dens_overlay")

#Model for 3 break points with variance among genera and multiple intercepts and slopes in both lines  - this corresponds to the growth stages
model_allometry_phylo_3bp_group_r = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_3bp_group_r <- mcp(model_allometry_phylo_3bp_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_3bp_group_r, facet_by="group")
summary(mcp_phylo_3bp_group_r)
ranef(mcp_phylo_3bp_group_r)
plot_pars(mcp_phylo_3bp_group_r, pars = c("cp_1", "cp_2","cp_3"), type = "dens_overlay")

#Model for 4 break points with variance among genera and multiple intercepts and slopes in both lines  - this corresponds to the growth stages
model_allometry_phylo_4bp_group_r = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_4bp_group_r <- mcp(model_allometry_phylo_4bp_group_r, allometry_phylo_rostrum_grp_cat_plot) 
plot(mcp_phylo_4bp_group_r, facet_by="group")
summary(mcp_phylo_4bp_group_r)
ranef(mcp_phylo_4bp_group_r)
plot_pars(mcp_phylo_4bp_group_r, pars = c("cp_1", "cp_2","cp_3","cp_4"), type = "dens_overlay")

#Find best fitting model - loo_compare
null_mcp_phylo_r$loo <- loo(null_mcp_phylo_r)
mcp_phylo_null_group_r$loo <- loo(mcp_phylo_null_group_r)
mcp_phylo_1bp_r$loo <- loo(mcp_phylo_1bp_r)
mcp_phylo_2bp_r$loo <- loo(mcp_phylo_2bp_r)
mcp_phylo_3bp_r$loo <- loo(mcp_phylo_3bp_r)
mcp_phylo_4bp_r$loo <- loo(mcp_phylo_4bp_r)
mcp_phylo_1bp_group_r$loo <- loo(mcp_phylo_1bp_group_r)
mcp_phylo_2bp_group_r$loo <- loo(mcp_phylo_2bp_group_r)
mcp_phylo_3bp_group_r$loo <- loo(mcp_phylo_3bp_group_r)
mcp_phylo_4bp_group_r$loo <- loo(mcp_phylo_4bp_group_r)

loo::loo_compare(null_mcp_phylo_r$loo, mcp_phylo_null_group_r$loo, #null models no bp 1,2
                 mcp_phylo_1bp_r$loo, mcp_phylo_2bp_r$loo, mcp_phylo_3bp_r$loo, mcp_phylo_4bp_r$loo, #bp same for all data 3,4,5,6
                 mcp_phylo_1bp_group_r$loo, mcp_phylo_2bp_group_r$loo, mcp_phylo_3bp_group_r$loo, mcp_phylo_4bp_group_r$loo) #bp by group 7,8,9,10

#Model with 3 bp by groups preferred, but 4 and 2 bp by groups next and not sign different (elpd_diff <4)
#Preference not very strong overall, large confidence intervals around cps

#Save results to file
sink("Output/8a-Allometry phylo/compare_mcps_phylo_rostrum_group.txt")
print("group and category reg scores used")

print("1-null model no bp and 1 slope, 2-null model slope and intercept different per group, 3-model 1 bp and 1 slope,\n
       4-model 2 bp and 1 slope, 5-model 3 bp and 1 slope, 6-model 4 bp and 1 slope,\n
       7-model 1 bp different slope and intercept in each group, 8-model 2 bp different slope and intercept in each group,\n 
      9-model 3 bp different slope and intercept in each group, 10-model 4 bp different slope and intercept in each group")
loo::loo_compare(null_mcp_phylo_r$loo, mcp_phylo_null_group_r$loo, mcp_phylo_1bp_r$loo, mcp_phylo_2bp_r$loo, mcp_phylo_3bp_r$loo,mcp_phylo_4bp_r$loo,
                 mcp_phylo_1bp_group_r$loo, mcp_phylo_2bp_group_r$loo, mcp_phylo_3bp_group_r$loo, mcp_phylo_4bp_group_r$loo)

print("Model with 3 bp by groups preferred, 2 bp and 4 bp next and not significantly different")

print("summary best model 3bp group")
summary(mcp_phylo_3bp_group_r)

print("summary 2nd best model 2bp group")
summary(mcp_phylo_2bp_group_r)

print("summary 3rd best model 4bp group")
summary(mcp_phylo_4bp_group_r)
sink()

#Get mean size min each growth stage fro each group - useful for plots
#New dataframe with mean by genus
mean_phylo_rostrum_group_category <- allometry_phylo_rostrum_grp_cat_plot %>% 
  group_by(group, category) %>% slice_min(logCS, n = 3) %>% summarize(mean = mean(logCS))

#Select cp means from summary
summary_mcp_phylo_rostrum <- summary(mcp_phylo_4bp_group_r)
#Create cp table
cp_phylo_mean_rostrum <- data.frame(cp = summary_mcp_phylo_rostrum$name[c(1,3,5,7)], mean = summary_mcp_phylo_rostrum$mean[c(1,3,5,7)])

#Get deltas for groups
cp_phylo_groups_rostrum <- as.data.frame(ranef(mcp_phylo_4bp_group_r))

#Create data frame and calculate cps for each group
cp_phylo_mean_groups_rostrum <- cp_phylo_mean_rostrum %>% bind_rows(replicate(1, cp_phylo_mean_rostrum, simplify = FALSE)) %>% arrange(mean)
cp_phylo_mean_groups_rostrum$group <- cp_phylo_groups_rostrum$name
cp_phylo_mean_groups_rostrum$delta <- cp_phylo_groups_rostrum$mean

cp_phylo_mean_groups_rostrum$group_mean <- rowSums(cp_phylo_mean_groups_rostrum[ , c("mean", "delta")])

#Match order to mean size data frame
cp_phylo_mean_groups_rostrum$group_name <- rep(c("mysticeti", "odontoceti"), times = 4)
cp_phylo_mean_groups_rostrum <- arrange(cp_phylo_mean_groups_rostrum, group_name)
cp_phylo_mean_groups_rostrum

#Make line graphs to compare cps and real data with error bars
#Create plot dataframe with mean, lower and upper
cp_phylo_all_rostrum  <- data.frame(cp = summary_mcp_phylo_rostrum$name[c(1,3,5,7)], 
                              mean = summary_mcp_phylo_rostrum$mean[c(1,3,5,7)],
                              lower = summary_mcp_phylo_rostrum$lower[c(1,3,5,7)],
                              upper = summary_mcp_phylo_rostrum$upper[c(1,3,5,7)])
#Create data frame and calculate cps for each group
cp_phylo_all_groups_rostrum <- cp_phylo_all_rostrum %>% bind_rows(replicate(1, cp_phylo_all_rostrum, simplify = FALSE)) %>% arrange(mean)
cp_phylo_all_groups_rostrum$group <- cp_phylo_groups_rostrum$name
cp_phylo_all_groups_rostrum$mean_delta <- cp_phylo_groups_rostrum$mean
cp_phylo_all_groups_rostrum$lower_delta <- cp_phylo_groups_rostrum$lower
cp_phylo_all_groups_rostrum$upper_delta <- cp_phylo_groups_rostrum$upper

cp_phylo_all_groups_rostrum$group_mean <- rowSums(cp_phylo_all_groups_rostrum[ , c("mean", "mean_delta")])
cp_phylo_all_groups_rostrum$group_lower <- rowSums(cp_phylo_all_groups_rostrum[ , c("lower", "lower_delta")])
cp_phylo_all_groups_rostrum$group_upper <- rowSums(cp_phylo_all_groups_rostrum[ , c("upper", "upper_delta")])

#Match order to mean size data frame
cp_phylo_all_groups_rostrum$group_name <- rep(c("mysticeti", "odontoceti"), times = 4)
cp_phylo_all_groups_rostrum <- arrange(cp_phylo_all_groups_rostrum, group_name)
cp_phylo_all_groups_rostrum 

summary_cp_phylo_all_groups_rostrum <- as_tibble(data.frame(group = cp_phylo_all_groups_rostrum$group_name, mean = cp_phylo_all_groups_rostrum$group_mean, 
                                                      lower = cp_phylo_all_groups_rostrum$group_lower,
                                                      upper = cp_phylo_all_groups_rostrum$group_upper))
summary_cp_phylo_all_groups_rostrum$type <- "estimates_mcp"
summary_cp_phylo_all_groups_rostrum$stage <- as.factor(rep(c("1","2","3","4"), times = 2))

#Get min and max values per group and category real data
min_max_rostrum_phylo_group_category <- allometry_phylo_rostrum_grp_cat_plot %>% 
  group_by(group, category) %>% summarize(lower = min(logCS), upper = max(logCS))

summary_phylo_rostrum_group_category <- cbind(mean_phylo_rostrum_group_category[,c(1,3)], min_max_rostrum_phylo_group_category[,3:4])
summary_phylo_rostrum_group_category$type <- "real_data"
summary_phylo_rostrum_group_category$stage <- as.factor(rep(c("1","2","3","4"), times = 2))

#Data frame for line plot
mcp_phylo_plot_rostrum <- bind_rows(summary_cp_phylo_all_groups_rostrum, summary_phylo_rostrum_group_category)

#Separate groups - easier plotting
mcp_phylo_plot_rostrum_mysticeti <- mcp_phylo_plot_rostrum[c(1:4,9:12),]
mcp_phylo_plot_rostrum_odontoceti <- mcp_phylo_plot_rostrum[-c(1:4,9:12),]

#Avoid error bars overlap, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

#Line plot comparing real estimated break values by group
real_mcp_phylo_breaks_ggplot_rostrum_myst <- ggplot(mcp_phylo_plot_rostrum_mysticeti, aes(x=stage, y=mean, colour=type, group=type)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="gray50", width=.2, position=pd) +
  geom_line(position=pd, linewidth = 1, show.legend = F) +
  geom_point(aes(shape=type, fill = type),position=pd, size=4, stroke = 1.5) + # 21 is filled circle
  xlab("Growth Stage") +
  ylab("Break value") +
  scale_colour_manual(name="Data type",    
                      breaks= levels(as.factor(mcp_phylo_plot_rostrum$type)),
                      labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                      values = c("darkorange3", "gray10"))+ 
  scale_shape_manual(name="Data type",   
                     labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                     values = c(23,24))+
  scale_fill_manual(name="Data type",    
                    labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                    values = alpha(c("darkorange3", "gray10"),0.7))+ 
  theme_bw(base_size = 12)+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11),  legend.position = "none",  legend.direction = "horizontal")
real_mcp_phylo_breaks_ggplot_rostrum_myst 

#Add phylopic
real_mcp_phylo_breaks_ggplot_rostrum_myst  <- 
  real_mcp_phylo_breaks_ggplot_rostrum_myst  +
  add_phylopic(myst, alpha = 1, x = 1.2, y = 4.2, ysize = 0.1, fill = "gray30")
real_mcp_phylo_breaks_ggplot_rostrum_myst 

#Line plot comparing real estimated break values by group
real_mcp_phylo_breaks_ggplot_rostrum_odont <- ggplot(mcp_phylo_plot_rostrum_odontoceti, aes(x=stage, y=mean, colour=type, group=type)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="gray50", width=.2, position=pd) +
  geom_line(position=pd, linewidth = 1, show.legend = F) +
  geom_point(aes(shape=type, fill = type),position=pd, size=4, stroke = 1.5) + # 21 is filled circle
  xlab("Growth Stage") +
  ylab("Break value") +
  scale_colour_manual(name="Data type",    
                      breaks= levels(as.factor(mcp_phylo_plot_rostrum$type)),
                      labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                      values = c("darkorange3", "gray10"))+ 
  scale_shape_manual(name="Data type",   
                     labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                     values = c(23,24))+
  scale_fill_manual(name="Data type",    
                    labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                    values = alpha(c("darkorange3", "gray10"),0.7))+ 
  theme_bw(base_size = 12)+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11),  legend.position = "none",  legend.direction = "horizontal")
real_mcp_phylo_breaks_ggplot_rostrum_odont 

#Add phylopic
real_mcp_phylo_breaks_ggplot_rostrum_odont  <- 
  real_mcp_phylo_breaks_ggplot_rostrum_odont  +
  add_phylopic(odont, alpha = 1, x = 1.2, y = 4.1, ysize = 0.1, fill = "gray50")
real_mcp_phylo_breaks_ggplot_rostrum_odont 

#Braincase ----

#Null model - no break points and single slope/intercept
null_model_allometry_phylo_b <- list(RegScores ~ logCS)

null_mcp_phylo_b <- mcp(null_model_allometry_phylo_b, allometry_phylo_braincase_grp_cat_plot) 
plot(null_mcp_phylo_b)
summary(null_mcp_phylo_b)

#Null model - no break point slope/intercept by group
model_allometry_phylo_null_group_b = list(RegScores ~ 1 + (1|group) + logCS)

mcp_phylo_null_group_b <- mcp(model_allometry_phylo_null_group_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_null_group_b)
summary(mcp_phylo_null_group_b)

#Model for single break point with multiple intercepts and slopes
model_allometry_phylo_1bp_b = list(RegScores ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_1bp_b <- mcp(model_allometry_phylo_1bp_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_1bp_b)
plot_pars(mcp_phylo_1bp_b, pars = "cp_1", type = "dens") + theme_bw(10)
summary(mcp_phylo_1bp_b)

#Model for 2 break points with multiple intercepts and slopes
model_allometry_phylo_2bp_b = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_2bp_b <- mcp(model_allometry_phylo_2bp_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_2bp_b)
plot_pars(mcp_phylo_2bp_b, pars = c("cp_1", "cp_2"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_2bp_b)

#Model for 3 break points with multiple intercepts and slopes
model_allometry_phylo_3bp_b = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_3bp_b <- mcp(model_allometry_phylo_3bp_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_3bp_b)
plot_pars(mcp_phylo_3bp_b, pars = c("cp_1", "cp_2", "cp_3"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_3bp_b)

#Model for 4 break points with multiple intercepts and slopes  - this corresponds to the growth stages
model_allometry_phylo_4bp_b = list(RegScores ~ 1 + logCS, ~ 1 + logCS,  ~ 1 + logCS, ~ 1 + logCS, ~ 1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_4bp_b <- mcp(model_allometry_phylo_4bp_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_4bp_b)
plot_pars(mcp_phylo_4bp_b, pars = c("cp_1", "cp_2", "cp_3", "cp_4"), type = "dens_overlay") + theme_bw(10)
summary(mcp_phylo_4bp_b)

#Model for single break point with variance among genera and multiple intercepts and slopes in both lines
model_allometry_phylo_1bp_group_b = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_1bp_group_b <- mcp(model_allometry_phylo_1bp_group_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_1bp_group_b, facet_by="group")
summary(mcp_phylo_1bp_group_b)
plot_pars(mcp_phylo_1bp_group_b, pars = "cp_1", type = "dens_overlay")

#Model for 2 break points with variance among genera and multiple intercepts and slopes in both lines
model_allometry_phylo_2bp_group_b = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_2bp_group_b <- mcp(model_allometry_phylo_2bp_group_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_2bp_group_b, facet_by="group")
summary(mcp_phylo_2bp_group_b)
ranef(mcp_phylo_2bp_group_b)
plot_pars(mcp_phylo_2bp_group_b, pars = c("cp_1", "cp_2"), type = "dens_overlay")

#Model for 3 break points with variance among genera and multiple intercepts and slopes in both lines  - this corresponds to the growth stages
model_allometry_phylo_3bp_group_b = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_3bp_group_b <- mcp(model_allometry_phylo_3bp_group_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_3bp_group_b, facet_by="group")
summary(mcp_phylo_3bp_group_b)
ranef(mcp_phylo_3bp_group_b)
plot_pars(mcp_phylo_3bp_group_b, pars = c("cp_1", "cp_2","cp_3"), type = "dens_overlay")

#Model for 4 break points with variance among genera and multiple intercepts and slopes in both lines  - this corresponds to the growth stages
model_allometry_phylo_4bp_group_b = list(RegScores ~ 1 + (1|group) + logCS, 
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS,
                                   RegScores ~ 1 + (1|group) ~  1 + logCS) #break point with different intercept, assigned based on data

mcp_phylo_4bp_group_b <- mcp(model_allometry_phylo_4bp_group_b, allometry_phylo_braincase_grp_cat_plot) 
plot(mcp_phylo_4bp_group_b, facet_by="group")
summary(mcp_phylo_4bp_group_b)
ranef(mcp_phylo_4bp_group_b)
plot_pars(mcp_phylo_4bp_group_b, pars = c("cp_1", "cp_2","cp_3","cp_4"), type = "dens_overlay")

#Find best fitting model - loo_compare
null_mcp_phylo_b$loo <- loo(null_mcp_phylo_b)
mcp_phylo_null_group_b$loo <- loo(mcp_phylo_null_group_b)
mcp_phylo_1bp_b$loo <- loo(mcp_phylo_1bp_b)
mcp_phylo_2bp_b$loo <- loo(mcp_phylo_2bp_b)
mcp_phylo_3bp_b$loo <- loo(mcp_phylo_3bp_b)
mcp_phylo_4bp_b$loo <- loo(mcp_phylo_4bp_b)
mcp_phylo_1bp_group_b$loo <- loo(mcp_phylo_1bp_group_b)
mcp_phylo_2bp_group_b$loo <- loo(mcp_phylo_2bp_group_b)
mcp_phylo_3bp_group_b$loo <- loo(mcp_phylo_3bp_group_b)
mcp_phylo_4bp_group_b$loo <- loo(mcp_phylo_4bp_group_b)

loo::loo_compare(null_mcp_phylo_b$loo, mcp_phylo_null_group_b$loo, #null models no bp 1,2
                 mcp_phylo_1bp_b$loo, mcp_phylo_2bp_b$loo, mcp_phylo_3bp_b$loo, mcp_phylo_4bp_b$loo, #bp same for all data 3,4,5,6
                 mcp_phylo_1bp_group_b$loo, mcp_phylo_2bp_group_b$loo, mcp_phylo_3bp_group_b$loo, mcp_phylo_4bp_group_b$loo) #bp by group 7,8,9,10

#Model with 3 bp by groups preferred, 4 bp second not sign difference. Strong preference on other models (elpd_diff >4)
#Large confidence intervals around cps

#Save results to file
sink("Output/8a-Allometry phylo/compare_mcps_phylo_braincase_group.txt")
print("group and category reg scores used")

print("1-null model no bp and 1 slope, 2-null model slope and intercept different per group, 3-model 1 bp and 1 slope,\n
       4-model 2 bp and 1 slope, 5-model 3 bp and 1 slope, 6-model 4 bp and 1 slope,\n
       7-model 1 bp different slope and intercept in each group, 8-model 2 bp different slope and intercept in each group,\n 
      9-model 3 bp different slope and intercept in each group, 10-model 4 bp different slope and intercept in each group")
loo::loo_compare(null_mcp_phylo_b$loo, mcp_phylo_null_group_b$loo, mcp_phylo_1bp_b$loo, mcp_phylo_2bp_b$loo, mcp_phylo_3bp_b$loo,mcp_phylo_4bp_b$loo,
                 mcp_phylo_1bp_group_b$loo, mcp_phylo_2bp_group_b$loo, mcp_phylo_3bp_group_b$loo, mcp_phylo_4bp_group_b$loo)

print("Model with 3 bp by groups preferred, but 4 bp next and not sign different (elpd_diff <4)")

print("summary best model 3bp group")
summary(mcp_phylo_3bp_group_b)

print("summary 2nd best model 4bp group")
summary(mcp_phylo_4bp_group_b)

sink()

#Get mean size min each growth stage fro each group - useful for plots
#New dataframe with mean by genus
mean_phylo_braincase_group_category <- allometry_phylo_braincase_grp_cat_plot %>% 
  group_by(group, category) %>% slice_min(logCS, n = 3) %>% summarize(mean = mean(logCS))

#Select cp means from summary
summary_mcp_phylo_braincase <- summary(mcp_phylo_4bp_group_b)
#Create cp table
cp_phylo_mean_braincase <- data.frame(cp = summary_mcp_phylo_braincase$name[c(1,3,5,7)], mean = summary_mcp_phylo_braincase$mean[c(1,3,5,7)])

#Get deltas for groups
cp_phylo_groups_braincase <- as.data.frame(ranef(mcp_phylo_4bp_group_r))

#Create data frame and calculate cps for each group
cp_phylo_mean_groups_braincase <- cp_phylo_mean_braincase %>% bind_rows(replicate(1, cp_phylo_mean_braincase, simplify = FALSE)) %>% arrange(mean)
cp_phylo_mean_groups_braincase$group <- cp_phylo_groups_braincase$name
cp_phylo_mean_groups_braincase$delta <- cp_phylo_groups_braincase$mean

cp_phylo_mean_groups_braincase$group_mean <- rowSums(cp_phylo_mean_groups_braincase[ , c("mean", "delta")])

#Match order to mean size data frame
cp_phylo_mean_groups_braincase$group_name <- rep(c("mysticeti", "odontoceti"), times = 4)
cp_phylo_mean_groups_braincase <- arrange(cp_phylo_mean_groups_braincase, group_name)
cp_phylo_mean_groups_braincase

#Make line graphs to compare cps and real data with error bars
#Create plot dataframe with mean, lower and upper
cp_phylo_all_braincase  <- data.frame(cp = summary_mcp_phylo_braincase$name[c(1,3,5,7)], 
                                mean = summary_mcp_phylo_braincase$mean[c(1,3,5,7)],
                                lower = summary_mcp_phylo_braincase$lower[c(1,3,5,7)],
                                upper = summary_mcp_phylo_braincase$upper[c(1,3,5,7)])
#Create data frame and calculate cps for each group
cp_phylo_all_groups_braincase <- cp_phylo_all_braincase %>% bind_rows(replicate(1, cp_phylo_all_braincase, simplify = FALSE)) %>% arrange(mean)
cp_phylo_all_groups_braincase$group <- cp_phylo_groups_braincase$name
cp_phylo_all_groups_braincase$mean_delta <- cp_phylo_groups_braincase$mean
cp_phylo_all_groups_braincase$lower_delta <- cp_phylo_groups_braincase$lower
cp_phylo_all_groups_braincase$upper_delta <- cp_phylo_groups_braincase$upper

cp_phylo_all_groups_braincase$group_mean <- rowSums(cp_phylo_all_groups_braincase[ , c("mean", "mean_delta")])
cp_phylo_all_groups_braincase$group_lower <- rowSums(cp_phylo_all_groups_braincase[ , c("lower", "lower_delta")])
cp_phylo_all_groups_braincase$group_upper <- rowSums(cp_phylo_all_groups_braincase[ , c("upper", "upper_delta")])

#Match order to mean size data frame
cp_phylo_all_groups_braincase$group_name <- rep(c("mysticeti", "odontoceti"), times = 4)
cp_phylo_all_groups_braincase <- arrange(cp_phylo_all_groups_braincase, group_name)
cp_phylo_all_groups_braincase 

summary_cp_phylo_all_groups_braincase <- as_tibble(data.frame(group = cp_phylo_all_groups_braincase$group_name, mean = cp_phylo_all_groups_braincase$group_mean, 
                                                        lower = cp_phylo_all_groups_braincase$group_lower,
                                                        upper = cp_phylo_all_groups_braincase$group_upper))
summary_cp_phylo_all_groups_braincase$type <- "estimates_mcp"
summary_cp_phylo_all_groups_braincase$stage <- as.factor(rep(c("1","2","3","4"), times = 2))

#Get min and max values per group and category real data
min_max_braincase_phylo_group_category <- allometry_phylo_braincase_grp_cat_plot %>% 
  group_by(group, category) %>% summarize(lower = min(logCS), upper = max(logCS))

summary_phylo_braincase_group_category <- cbind(mean_phylo_braincase_group_category[,c(1,3)], min_max_braincase_phylo_group_category[,3:4])
summary_phylo_braincase_group_category$type <- "real_data"
summary_phylo_braincase_group_category$stage <- as.factor(rep(c("1","2","3","4"), times = 2))

#Data frame for line plot
mcp_phylo_plot_braincase <- bind_rows(summary_cp_phylo_all_groups_braincase, summary_phylo_braincase_group_category)

#Separate groups - easier plotting
mcp_phylo_plot_braincase_mysticeti <- mcp_phylo_plot_braincase[c(1:4,9:12),]
mcp_phylo_plot_braincase_odontoceti <- mcp_phylo_plot_braincase[-c(1:4,9:12),]

#Line plot comparing real estimated break values by group
real_mcp_phylo_breaks_ggplot_braincase_myst <- ggplot(mcp_phylo_plot_braincase_mysticeti, aes(x=stage, y=mean, colour=type, group=type)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="gray50", width=.2, position=pd) +
  geom_line(position=pd, linewidth = 1, show.legend = F) +
  geom_point(aes(shape=type, fill = type),position=pd, size=4, stroke = 1.5) + # 21 is filled circle
  xlab("Growth Stage") +
  ylab("Break value") +
  scale_colour_manual(name="Data type",    
                      breaks= levels(as.factor(mcp_phylo_plot_braincase$type)),
                      labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                      values = c("darkorange3", "gray10"))+ 
  scale_shape_manual(name="Data type",   
                     labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                     values = c(23,24))+
  scale_fill_manual(name="Data type",    
                    labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                    values = alpha(c("darkorange3", "gray10"),0.7))+ 
  theme_bw(base_size = 12)+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11),  legend.position = "bottom",  legend.direction = "horizontal")
real_mcp_phylo_breaks_ggplot_braincase_myst 

#Add phylopic
real_mcp_phylo_breaks_ggplot_braincase_myst  <- 
  real_mcp_phylo_breaks_ggplot_braincase_myst  +
  add_phylopic(myst, alpha = 1, x = 1.2, y = 4, ysize = 0.1, fill = "gray30")
real_mcp_phylo_breaks_ggplot_braincase_myst 

#Line plot comparing real estimated break values by group
real_mcp_phylo_breaks_ggplot_braincase_odont <- ggplot(mcp_phylo_plot_braincase_odontoceti, aes(x=stage, y=mean, colour=type, group=type)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="gray50", width=.2, position=pd) +
  geom_line(position=pd, linewidth = 1, show.legend = F) +
  geom_point(aes(shape=type, fill = type),position=pd, size=4, stroke = 1.5) + # 21 is filled circle
  xlab("Growth Stage") +
  ylab("Break value") +
  scale_colour_manual(name="Data type",    
                      breaks= levels(as.factor(mcp_phylo_plot_braincase$type)),
                      labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                      values = c("darkorange3", "gray10"))+ 
  scale_shape_manual(name="Data type",   
                     labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                     values = c(23,24))+
  scale_fill_manual(name="Data type",    
                    labels=c("MCP estimated breaks", "Real data (average minimum size for each growth stage)"),
                    values = alpha(c("darkorange3", "gray10"),0.7))+ 
  theme_bw(base_size = 12)+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 11),  legend.position = "bottom",  legend.direction = "horizontal")
real_mcp_phylo_breaks_ggplot_braincase_odont 

#Add phylopic
real_mcp_phylo_breaks_ggplot_braincase_odont  <- 
  real_mcp_phylo_breaks_ggplot_braincase_odont  +
  add_phylopic(odont, alpha = 1, x = 1.2, y = 3.7, ysize = 0.1, fill = "gray50")
real_mcp_phylo_breaks_ggplot_braincase_odont 

plotMCP1<- ggarrange(real_mcp_phylo_breaks_ggplot_braincase_myst , real_mcp_phylo_breaks_ggplot_rostrum_myst,
                     nrow= 1, ncol = 2, common.legend = T, legend = "bottom")
plotMCP1  


plotMCP2<- ggarrange(real_mcp_phylo_breaks_ggplot_braincase_odont , real_mcp_phylo_breaks_ggplot_rostrum_odont,
                     nrow= 1, ncol = 2, common.legend = T, legend = "bottom")
plotMCP2 

plotMCP3 <- ggarrange(plotMCP1, plotMCP2, nrow=2, ncol = 1, common.legend = T, legend = "bottom")

plotMCP3<- annotate_figure(plotMCP3, top = text_grob("Braincase and Rostrum", 
                                                     face = "bold", size = 17))
plotMCP3

#ALLOMETRY DIFFERENCE ROSTRUM AND BRAINCASE COMPARE----
##Evaluate allometry using logCS for entire skull as size
##Use pcscores instead of coords - all same length of variables and same results as coords

##By group and module ----
#Test slope differences between whole skull, rostrum and braincase allometry per stage
allometry_phylo_group_module_2_null <- lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group, iter=999, print.progress = TRUE) 
allometry_phylo_group_module_2_comb <-  lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group + phylo_pcscores_R_B$module, iter=999, print.progress = TRUE) 
allometry_phylo_group_module_2_int <-  lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group * phylo_pcscores_R_B$module, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
anova(allometry_phylo_group_module_2_null)
anova(allometry_phylo_group_module_2_comb)
anova(allometry_phylo_group_module_2_int)

#Save results of significant regression to file
sink("Output/8a-Allometry phylo/allometry_phylo_group_module_2_models.txt")
print("Null")
anova(allometry_phylo_group_module_2_null)

print("Combination +")
anova(allometry_phylo_group_module_2_comb) 

print("Interaction *")
anova(allometry_phylo_group_module_2_int)
sink() 

#ANOVA models
anova_allometry_phylo_group_module_2_models <- anova(allometry_phylo_group_module_2_null, allometry_phylo_group_module_2_comb, allometry_phylo_group_module_2_int)
anova_allometry_phylo_group_module_2_models

#Create interaction object
module_phylo_2_group <- interaction(phylo_pcscores_R_B$group,phylo_pcscores_R_B$module,sep = "_")

#Pairwise test slopes
pairwise_allometry_phylo_group_module_2 <- pairwise(allometry_phylo_group_module_2_int, fit.null = allometry_phylo_group_module_2_comb,
                                              groups = module_phylo_2_group, 
                                              covariate = phylo_pcscores_R_B$size, print.progress = FALSE) 
pairwise_allometry_phylo_group_module_2

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_phylo_group_module_2_distance <- summary(pairwise_allometry_phylo_group_module_2, confidence = 0.95, test.type = "dist") 
pairwise_allometry_phylo_group_module_2_distance

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_phylo_group_module_2_VC <- summary(pairwise_allometry_phylo_group_module_2, confidence = 0.95, test.type = "VC",
                                                angle.type = "deg") 
pairwise_allometry_phylo_group_module_2_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_phylo_group_module_2_DL <-summary(pairwise_allometry_phylo_group_module_2, confidence = 0.95, test.type = "DL") 
pairwise_allometry_phylo_group_module_2_DL 

#Compare the dispersion around group slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
pairwise_allometry_phylo_group_module_2_var <-summary(pairwise_allometry_phylo_group_module_2, confidence = 0.95, test.type = "var")
pairwise_allometry_phylo_group_module_2_var

#Save results to file
sink("Output/8a-Allometry phylo/pairwise_allometry_phylo_group_module_2.txt")
print("ANOVA models")
print(anova_allometry_phylo_group_module_2_models)

print("1-Pairwise absolute distances slopes")
print(pairwise_allometry_phylo_group_module_2_distance)

print("2-Distance between angles (slope directions)")
print(pairwise_allometry_phylo_group_module_2_VC)

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
print(pairwise_allometry_phylo_group_module_2_DL)

print("4-Difference in dispersion around mean slope")
print(pairwise_allometry_phylo_group_module_2_var)
sink()

###Heatmaps plots for significant differences in pairwise ----

#Save p-values as object
pairwise_allometry_phylo_group_module_2_dist <- pairwise_allometry_phylo_group_module_2_distance[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_group_module_2_dist_p <- pairwise_allometry_phylo_group_module_2_distance[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_group_module_2_angle <- pairwise_allometry_phylo_group_module_2_VC[["pairwise.tables"]][["angle"]]
pairwise_allometry_phylo_group_module_2_angle_p <- pairwise_allometry_phylo_group_module_2_VC[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_group_module_2_length <- pairwise_allometry_phylo_group_module_2_DL[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_group_module_2_length_p <- pairwise_allometry_phylo_group_module_2_DL[["pairwise.tables"]][["P"]]

#Make list to change tables faster
pairwise_allometry_phylo_group_module_2_list <- list(pairwise_allometry_phylo_group_module_2_dist, pairwise_allometry_phylo_group_module_2_dist_p, pairwise_allometry_phylo_group_module_2_angle, pairwise_allometry_phylo_group_module_2_angle_p, 
                                               pairwise_allometry_phylo_group_module_2_length, pairwise_allometry_phylo_group_module_2_length_p)

#Make new list of variable names including modules and genera
#Save row and col names as variables to change string - colnames = rownames for both
group_module_phylo_2_vars <- rownames(pairwise_allometry_phylo_group_module_2_dist_p)

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  group_module_phylo_2_vars <- str_replace_all(group_module_phylo_2_vars, modules_2_list[m], modules_2_list_short[m])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  group_module_phylo_2_vars <- str_replace_all(group_module_phylo_2_vars, groups_list[t], groups_list_short[t])
}

#Loop change . to _
group_module_phylo_2_vars <- str_replace(group_module_phylo_2_vars, "[.]", "_")

#Check it worked
group_module_phylo_2_vars

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  rownames(pairwise_allometry_phylo_group_module_2_list[[l]]) <- group_module_phylo_2_vars
  colnames(pairwise_allometry_phylo_group_module_2_list[[l]]) <- group_module_phylo_2_vars
}

#Save only lower triangle for each
pairwise_allometry_phylo_group_module_2_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_group_module_2_lower_tri_list[[l]] <- get_upper_tri(pairwise_allometry_phylo_group_module_2_list[[l]])
}

#Melt to make table in the format needed for heatmap
pairwise_allometry_phylo_group_module_2_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_group_module_2_melt[[l]] <- melt(pairwise_allometry_phylo_group_module_2_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
pairwise_allometry_phylo_group_module_2_dist_melt <- data.frame(pairwise_allometry_phylo_group_module_2_melt[[1]], p = pairwise_allometry_phylo_group_module_2_melt[[2]][[3]])
pairwise_allometry_phylo_group_module_2_angle_melt <- data.frame(pairwise_allometry_phylo_group_module_2_melt[[3]], p = pairwise_allometry_phylo_group_module_2_melt[[4]][[3]])
pairwise_allometry_phylo_group_module_2_length_melt <- data.frame(pairwise_allometry_phylo_group_module_2_melt[[5]], p = pairwise_allometry_phylo_group_module_2_melt[[6]][[3]])

#Create columns where only significant values are shown
pairwise_allometry_phylo_group_module_2_dist_melt <- pairwise_allometry_phylo_group_module_2_dist_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                      p_if_sig = ifelse(sig_p, p, NA),
                                                                                                      value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))
pairwise_allometry_phylo_group_module_2_angle_melt <- pairwise_allometry_phylo_group_module_2_angle_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 1)))
pairwise_allometry_phylo_group_module_2_length_melt <- pairwise_allometry_phylo_group_module_2_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                                                          value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(pairwise_allometry_phylo_group_module_2_dist_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_group_module_2_angle_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_group_module_2_length_melt$p_if_sig))

#Nice heatmap plot for each variable
pairwise_allometry_phylo_group_module_2_dist_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_group_module_2_dist_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[2], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.0499), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13,  hjust = 0.75, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.2,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_group_module_2_dist_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_group_module_2_angle_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_group_module_2_angle_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[2], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.0499), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13,  hjust = 0.75, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
pairwise_allometry_phylo_group_module_2_angle_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_group_module_2_length_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_group_module_2_length_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[2], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.0499), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope length difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13,  hjust = 0.75, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
pairwise_allometry_phylo_group_module_2_length_heatmap_ggplot

plotA <- ggarrange(pairwise_allometry_phylo_group_module_2_dist_heatmap_ggplot,  pairwise_allometry_phylo_group_module_2_angle_heatmap_ggplot, pairwise_allometry_phylo_group_module_2_length_heatmap_ggplot,
                   nrow = 1, ncol = 3, common.legend = F)
plotA <- annotate_figure(plotA, top = text_grob("Allometry phylogentically corrected by group and module", face = "bold", size = 17))
plotA

##By group, category and module ----
#Test slope differences between whole skull, rostrum and braincase allometry per stage
allometry_phylo_group_cat_module_2_null <- lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group_cat, iter=999, print.progress = TRUE) 
allometry_phylo_group_cat_module_2_comb <-  lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group_cat + phylo_pcscores_R_B$module, iter=999, print.progress = TRUE) 
allometry_phylo_group_cat_module_2_int <-  lm.rrpp(phylo_pcscores_R_B_matrix ~ phylo_pcscores_R_B$size * phylo_pcscores_R_B$group_cat * phylo_pcscores_R_B$module, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
anova(allometry_phylo_group_cat_module_2_null)
anova(allometry_phylo_group_cat_module_2_comb)
anova(allometry_phylo_group_cat_module_2_int)

#Save results of significant regression to file
sink("Output/8a-Allometry phylo/allometry_phylo_group_cat_module_2_models.txt")
print("Null")
anova(allometry_phylo_group_cat_module_2_null)

print("Combination +")
anova(allometry_phylo_group_cat_module_2_comb) 

print("Interaction *")
anova(allometry_phylo_group_cat_module_2_int)
sink() 

#ANOVA models
anova_allometry_phylo_group_cat_module_2_models <- anova(allometry_phylo_group_cat_module_2_null, allometry_phylo_group_cat_module_2_comb, allometry_phylo_group_cat_module_2_int)
anova_allometry_phylo_group_cat_module_2_models

#Create interaction object
module_phylo_2_group_cat <- interaction(phylo_pcscores_R_B$group_cat,phylo_pcscores_R_B$module,sep = "_")


#Pairwise test slopes
pairwise_allometry_phylo_group_cat_module_2 <- pairwise(allometry_phylo_group_cat_module_2_int, fit.null = allometry_phylo_group_cat_module_2_comb,
                                                  groups = module_phylo_2_group_cat, 
                                                  covariate = phylo_pcscores_R_B$size, print.progress = FALSE) 
pairwise_allometry_phylo_group_cat_module_2

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_phylo_group_cat_module_2_distance <- summary(pairwise_allometry_phylo_group_cat_module_2, confidence = 0.95, test.type = "dist") 
pairwise_allometry_phylo_group_cat_module_2_distance

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_phylo_group_cat_module_2_VC <- summary(pairwise_allometry_phylo_group_cat_module_2, confidence = 0.95, test.type = "VC",
                                                    angle.type = "deg") 
pairwise_allometry_phylo_group_cat_module_2_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_phylo_group_cat_module_2_DL <-summary(pairwise_allometry_phylo_group_cat_module_2, confidence = 0.95, test.type = "DL") 
pairwise_allometry_phylo_group_cat_module_2_DL 

#Compare the dispersion around group slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
pairwise_allometry_phylo_group_cat_module_2_var <-summary(pairwise_allometry_phylo_group_cat_module_2, confidence = 0.95, test.type = "var")
pairwise_allometry_phylo_group_cat_module_2_var

#Save results to file
sink("Output/8a-Allometry phylo/pairwise_allometry_phylo_group_cat_module_2.txt")
print("ANOVA models")
print(anova_allometry_phylo_group_cat_module_2_models)

print("1-Pairwise absolute distances slopes")
print(pairwise_allometry_phylo_group_cat_module_2_distance)

print("2-Distance between angles (slope directions)")
print(pairwise_allometry_phylo_group_cat_module_2_VC)

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
print(pairwise_allometry_phylo_group_cat_module_2_DL)

print("4-Difference in dispersion around mean slope")
print(pairwise_allometry_phylo_group_cat_module_2_var)
sink()

###Heatmaps plots for significant differences in pairwise ----

#Save p-values as object
pairwise_allometry_phylo_group_cat_module_2_dist <- pairwise_allometry_phylo_group_cat_module_2_distance[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_group_cat_module_2_dist_p <- pairwise_allometry_phylo_group_cat_module_2_distance[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_group_cat_module_2_angle <- pairwise_allometry_phylo_group_cat_module_2_VC[["pairwise.tables"]][["angle"]]
pairwise_allometry_phylo_group_cat_module_2_angle_p <- pairwise_allometry_phylo_group_cat_module_2_VC[["pairwise.tables"]][["P"]]
pairwise_allometry_phylo_group_cat_module_2_length <- pairwise_allometry_phylo_group_cat_module_2_DL[["pairwise.tables"]][["D"]]
pairwise_allometry_phylo_group_cat_module_2_length_p <- pairwise_allometry_phylo_group_cat_module_2_DL[["pairwise.tables"]][["P"]]

#Make list to change tables faster
pairwise_allometry_phylo_group_cat_module_2_list <- list(pairwise_allometry_phylo_group_cat_module_2_dist, pairwise_allometry_phylo_group_cat_module_2_dist_p, pairwise_allometry_phylo_group_cat_module_2_angle, pairwise_allometry_phylo_group_cat_module_2_angle_p, 
                                                   pairwise_allometry_phylo_group_cat_module_2_length, pairwise_allometry_phylo_group_cat_module_2_length_p)

#Make new list of variable names including modules and genera
#Save row and col names as variables to change string - colnames = rownames for both
group_cat_module_phylo_2_vars <- rownames(pairwise_allometry_phylo_group_cat_module_2_dist_p)

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  group_cat_module_phylo_2_vars <- str_replace_all(group_cat_module_phylo_2_vars, modules_2_list[m], modules_2_list_short[m])
}

#Loop replacements categories
for (s in 1:length(categories_list)){
  group_cat_module_phylo_2_vars <- str_replace_all(group_cat_module_phylo_2_vars, categories_list[s], categories_list_short[s])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  group_cat_module_phylo_2_vars <- str_replace_all(group_cat_module_phylo_2_vars, groups_list[t], groups_list_short[t])
}

#Loop change . to _
group_cat_module_phylo_2_vars <- str_replace(group_cat_module_phylo_2_vars, "[.]", "_")

#Check it worked
group_cat_module_phylo_2_vars

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  rownames(pairwise_allometry_phylo_group_cat_module_2_list[[l]]) <- group_cat_module_phylo_2_vars
  colnames(pairwise_allometry_phylo_group_cat_module_2_list[[l]]) <- group_cat_module_phylo_2_vars
}

#Save only lower triangle for each
pairwise_allometry_phylo_group_cat_module_2_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_group_cat_module_2_lower_tri_list[[l]] <- get_upper_tri(pairwise_allometry_phylo_group_cat_module_2_list[[l]])
}

#Melt to make table in the format needed for heatmap
pairwise_allometry_phylo_group_cat_module_2_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by pairwise results
  pairwise_allometry_phylo_group_cat_module_2_melt[[l]] <- melt(pairwise_allometry_phylo_group_cat_module_2_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
pairwise_allometry_phylo_group_cat_module_2_dist_melt <- data.frame(pairwise_allometry_phylo_group_cat_module_2_melt[[1]], p = pairwise_allometry_phylo_group_cat_module_2_melt[[2]][[3]])
pairwise_allometry_phylo_group_cat_module_2_angle_melt <- data.frame(pairwise_allometry_phylo_group_cat_module_2_melt[[3]], p = pairwise_allometry_phylo_group_cat_module_2_melt[[4]][[3]])
pairwise_allometry_phylo_group_cat_module_2_length_melt <- data.frame(pairwise_allometry_phylo_group_cat_module_2_melt[[5]], p = pairwise_allometry_phylo_group_cat_module_2_melt[[6]][[3]])

#Create columns where only significant values are shown
pairwise_allometry_phylo_group_cat_module_2_dist_melt <- pairwise_allometry_phylo_group_cat_module_2_dist_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                                                              value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 1)))
pairwise_allometry_phylo_group_cat_module_2_angle_melt <- pairwise_allometry_phylo_group_cat_module_2_angle_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 0)))
pairwise_allometry_phylo_group_cat_module_2_length_melt <- pairwise_allometry_phylo_group_cat_module_2_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                                                                  value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 1)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(pairwise_allometry_phylo_group_cat_module_2_dist_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_group_cat_module_2_angle_melt$p_if_sig))

all(is.na(pairwise_allometry_phylo_group_cat_module_2_length_melt$p_if_sig))

#Nice heatmap plot for each variable
pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_dist_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 3.5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[2], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.0499), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13,  hjust = 0.75, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.2,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_angle_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 3.5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[2], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.0499), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x =  element_text(angle = 45, size = 13,  hjust = 0.75, vjust = 0.8),
        axis.text.y =  element_text(size = 13, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot

ggarrange(pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot,  pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot, 
          nrow = 1, ncol = 2, common.legend = F)

####Plots by group ----
#Make data frame for each stage
pairwise_allometry_phylo_group_cat_module_2_dist_melt_mysticeti <- pairwise_allometry_phylo_group_cat_module_2_dist_melt %>% filter(str_detect(Var1, "Myst")) %>% 
  filter(str_detect(Var2, "Myst"))
pairwise_allometry_phylo_group_cat_module_2_dist_melt_odontoceti <- pairwise_allometry_phylo_group_cat_module_2_dist_melt %>% filter(str_detect(Var1, "Odont")) %>% 
  filter(str_detect(Var2, "Odont"))

pairwise_allometry_phylo_group_cat_module_2_angle_melt_mysticeti <- pairwise_allometry_phylo_group_cat_module_2_angle_melt %>% filter(str_detect(Var1, "Myst")) %>% 
  filter(str_detect(Var2, "Myst"))
pairwise_allometry_phylo_group_cat_module_2_angle_melt_odontoceti <- pairwise_allometry_phylo_group_cat_module_2_angle_melt %>% filter(str_detect(Var1, "Odont")) %>% 
  filter(str_detect(Var2, "Odont"))

#Nice heatmap plot
pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_mysticeti  <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_dist_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4)+
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  scale_x_discrete(labels = category_modules_2_list)+
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.3,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_mysticeti 

pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_odontoceti  <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_dist_melt_odontoceti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4)+
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+   
  scale_x_discrete(labels = category_modules_2_list)+   
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.3,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_odontoceti

pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_mysticeti  <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_angle_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4)+
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+   
  scale_x_discrete(labels = category_modules_2_list)+   
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position.inside = c(0.18,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_mysticeti 

pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_odontoceti  <- ggplot(data = pairwise_allometry_phylo_group_cat_module_2_angle_melt_odontoceti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4)+
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+   
  scale_x_discrete(labels = category_modules_2_list)+   
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position.inside = c(0.18,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_odontoceti

plotM2<-ggarrange(pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_mysticeti,pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_mysticeti,
                  ncol = 2, nrow = 1, common.legend = F)
plotM2<-annotate_figure(plotM2, top = text_grob("Mysticeti", face = "bold", size = 17))
plotM2

plotO2<-ggarrange(pairwise_allometry_phylo_group_cat_module_2_dist_heatmap_ggplot_odontoceti,pairwise_allometry_phylo_group_cat_module_2_angle_heatmap_ggplot_odontoceti,
                  ncol = 2, nrow = 1, common.legend = F)
plotO2<-annotate_figure(plotO2, top = text_grob("Odontoceti", face = "bold", size = 17))
plotO2

ggarrange(plotM2, plotO2, ncol = 1, nrow = 2, common.legend = T)


###Plot allometry rostrum, braincase divided by group and category ----

#Regression score of shape vs logCS and comb or int (best model)- regression method with "RegScore" plotting
allometry_phylo_group_cat_module_2_plot_regscore <- plot(allometry_phylo_group_cat_module_2_int, type = "regression",predictor = phylo_pcscores_R_B$size, reg.type = "RegScore",
                                                   main = "Shape vs logCS phylo corrected * group_cat * module", xlab = "logCS", pch = 21, col = "violet", 
                                                   bg = "violet", cex = 1.2, font.main = 2)   #improve graphics
text(x = phylo_pcscores_R_B$size, y = allometry_phylo_group_cat_module_2_plot_regscore$RegScore, labels = phylo_pcscores_R_B$specimens,
     pos = 3, offset = 0.5, cex = 0.5)    #improve appearance of labels

#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_phylo_group_cat_module_2_plot <- data.frame(logCS = allometry_phylo_group_cat_module_2_plot_regscore[["plot_args"]][["x"]], 
                                                RegScores = allometry_phylo_group_cat_module_2_plot_regscore[["plot_args"]][["y"]])
#Convert data frame to tibble
allometry_phylo_group_cat_module_2_plot <- as_tibble(allometry_phylo_group_cat_module_2_plot)
#Add labels and other attributes to tibble as columns
allometry_phylo_group_cat_module_2_plot <- allometry_phylo_group_cat_module_2_plot %>% 
  mutate(specimens = phylo_pcscores_R_B$specimens,  family = phylo_pcscores_R_B$family, category = phylo_pcscores_R_B$category, group = phylo_pcscores_R_B$group,
         group_cat = phylo_pcscores_R_B$group_cat, module =  str_to_sentence(phylo_pcscores_R_B$module))
allometry_phylo_group_cat_module_2_plot$group_cat_module <- module_phylo_2_group_cat
glimpse(allometry_phylo_group_cat_module_2_plot)

####By group ----
#Divide by groups - too many lines all together
allometry_phylo_group_cat_module_2_plot_mysticeti <- allometry_phylo_group_cat_module_2_plot %>% filter(group == "mysticeti")
allometry_phylo_group_cat_module_2_plot_odontoceti <- allometry_phylo_group_cat_module_2_plot %>% filter(group == "odontoceti")

#Plot allometry regression by category mysticeti
allometry_phylo_group_cat_module_2_ggplot_mysticeti <- ggplot(allometry_phylo_group_cat_module_2_plot_mysticeti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = category, color = module), inherit.aes = F,       
              linewidth = 1.2, alpha =1, se = F, show.legend = T)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(shape = 21, linetype = 0, alpha =1, size = 4)), linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_group_cat_module_2_ggplot_mysticeti

#Add phylopic
allometry_phylo_group_cat_module_2_ggplot_mysticeti <- 
  allometry_phylo_group_cat_module_2_ggplot_mysticeti +
  add_phylopic(myst, alpha = 1, x = 3, y = 0.75, ysize = 0.1, fill = "gray30")
allometry_phylo_group_cat_module_2_ggplot_mysticeti
#Fix silhouette size

allometry_phylo_group_cat_module_2_ggplot_odontoceti <- ggplot(allometry_phylo_group_cat_module_2_plot_odontoceti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module), alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = category, color = module), inherit.aes = F,       
              linewidth = 1.2, alpha =1, se = F, show.legend = T)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module), scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(shape = 21, linetype = 0, alpha =1, size = 4)), linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_group_cat_module_2_ggplot_odontoceti

#Add phylopic
allometry_phylo_group_cat_module_2_ggplot_odontoceti <- 
  allometry_phylo_group_cat_module_2_ggplot_odontoceti +
  add_phylopic(odont, alpha = 1, x = 2.75, y = 0.25, ysize = 0.09, fill = "gray50")
allometry_phylo_group_cat_module_2_ggplot_odontoceti
#Fix silhouette size

ggarrange(allometry_phylo_group_cat_module_2_ggplot_mysticeti,allometry_phylo_group_cat_module_2_ggplot_odontoceti, 
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

###Line plots ----
#geom_abline plot by group faceted by module

#Make separate dfs for each module - easier to extract coefficients correctly
allometry_phylo_group_cat_module_2_plot_braincase <- allometry_phylo_group_cat_module_2_plot %>% filter(module == "Braincase")
allometry_phylo_group_cat_module_2_plot_rostrum <- allometry_phylo_group_cat_module_2_plot %>% filter(module == "Rostrum")

####By group and module ----
#Braincase
#Linear model for line by group and category
allometry_phylo_group_regline_braincase <- lm(RegScores ~ logCS * group, data = allometry_phylo_group_cat_module_2_plot_braincase)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_group_regline_coeffs_braincase <- as.matrix(allometry_phylo_group_regline_braincase$coefficients)

#Save intercepts and slopes separately
allometry_phylo_group_regline_intercepts_braincase <- as.matrix(allometry_phylo_group_regline_coeffs_braincase[c(1, 3:(length(levels(groups))+1)),])
allometry_phylo_group_regline_slopes_braincase <- as.matrix(allometry_phylo_group_regline_coeffs_braincase[c(2, length(levels(groups))+2:(length(levels(groups)))),])

#Calculate real intercepts and slopes
allometry_phylo_group_regline_intercepts_ok_braincase <- as.matrix(c(allometry_phylo_group_regline_intercepts_braincase[1,], allometry_phylo_group_regline_intercepts_braincase[1,]+
                                                                 allometry_phylo_group_regline_intercepts_braincase[2:length(allometry_phylo_group_regline_intercepts_braincase),]))

allometry_phylo_group_regline_slopes_ok_braincase <- as.matrix(c(allometry_phylo_group_regline_slopes_braincase[1,], allometry_phylo_group_regline_slopes_braincase[1,]+
                                                             allometry_phylo_group_regline_slopes_braincase[2:length(allometry_phylo_group_regline_slopes_braincase),]))

#Save as data frame with grouping variables
allometry_phylo_group_coeffs_braincase <- data.frame(Slope = allometry_phylo_group_regline_slopes_ok_braincase, Intercept = allometry_phylo_group_regline_intercepts_ok_braincase, 
                                               row.names = levels(groups))
#Check for NA and other issues
allometry_phylo_group_coeffs_braincase 

#Add classifiers
allometry_phylo_group_coeffs_braincase <- allometry_phylo_group_coeffs_braincase %>% mutate(group = levels(groups),
                                                                                module = 'Braincase')
allometry_phylo_group_coeffs_braincase

#Rostrum
#Linear model for line by group and category
allometry_phylo_group_regline_rostrum <- lm(RegScores ~ logCS * group, data = allometry_phylo_group_cat_module_2_plot_rostrum)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_group_regline_coeffs_rostrum <- as.matrix(allometry_phylo_group_regline_rostrum$coefficients)

#Save intercepts and slopes separately
allometry_phylo_group_regline_intercepts_rostrum <- as.matrix(allometry_phylo_group_regline_coeffs_rostrum[c(1, 3:(length(levels(groups))+1)),])
allometry_phylo_group_regline_slopes_rostrum <- as.matrix(allometry_phylo_group_regline_coeffs_rostrum[c(2, length(levels(groups))+2:(length(levels(groups)))),])

#Calculate real intercepts and slopes
allometry_phylo_group_regline_intercepts_ok_rostrum <- as.matrix(c(allometry_phylo_group_regline_intercepts_rostrum[1,], allometry_phylo_group_regline_intercepts_rostrum[1,]+
                                                               allometry_phylo_group_regline_intercepts_rostrum[2:length(allometry_phylo_group_regline_intercepts_rostrum),]))

allometry_phylo_group_regline_slopes_ok_rostrum <- as.matrix(c(allometry_phylo_group_regline_slopes_rostrum[1,], allometry_phylo_group_regline_slopes_rostrum[1,]+
                                                           allometry_phylo_group_regline_slopes_rostrum[2:length(allometry_phylo_group_regline_slopes_rostrum),]))

#Save as data frame with grouping variables
allometry_phylo_group_coeffs_rostrum <- data.frame(Slope = allometry_phylo_group_regline_slopes_ok_rostrum, Intercept = allometry_phylo_group_regline_intercepts_ok_rostrum, 
                                             row.names = levels(groups))
#Check for NA and other issues
allometry_phylo_group_coeffs_rostrum 

#Add classifiers
allometry_phylo_group_coeffs_rostrum <- allometry_phylo_group_coeffs_rostrum %>% mutate(group = levels(groups),
                                                                            module = 'Rostrum')
allometry_phylo_group_coeffs_rostrum

#Put together
allometry_phylo_group_coeffs <- rbind(allometry_phylo_group_coeffs_braincase, allometry_phylo_group_coeffs_rostrum)
allometry_phylo_group_coeffs

#Make sure group correspond to use facet_wrap
allometry_phylo_group_cat_module_2_plot$group <- str_to_title(allometry_phylo_group_cat_module_2_plot$group)

##Plot by module both groups
allometry_phylo_grp_line_ggplot <- ggplot(allometry_phylo_group_cat_module_2_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_group_coeffs, 
              aes(intercept = Intercept, slope = Slope,  colour = group, linetype = module), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Groups", labels =  levels(groups), 
                      values = mypalette_groups, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Modules",
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0))+
  guides(colour = "none",
         linetype = guide_legend(keywidth = unit(3, "char"), override.aes = list(colour = c("gray20"))))
allometry_phylo_grp_line_ggplot

#Add silhouettes groups
allometry_phylo_grp_line_ggplot  <- 
  allometry_phylo_grp_line_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 3.2, y = -0.5, ysize = 0.15, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 2.8, y = 0.1, ysize = 0.15, fill = mypalette_groups[2])
allometry_phylo_grp_line_ggplot

####By group and category by module ----
#Braincase
#Linear model for line by group and category
allometry_phylo_group_cat_regline_braincase <- lm(RegScores ~ logCS * group_cat, data = allometry_phylo_group_cat_module_2_plot_braincase)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_group_cat_regline_coeffs_braincase <- as.matrix(allometry_phylo_group_cat_regline_braincase$coefficients)

#Save intercepts and slopes separately
group_cat_vars <- levels(allometry_phylo_group_cat_module_2_plot$group_cat)

allometry_phylo_group_cat_regline_intercepts_braincase <- as.matrix(allometry_phylo_group_cat_regline_coeffs_braincase[c(1, 3:(length(group_cat_vars)+1)),])
allometry_phylo_group_cat_regline_slopes_braincase <- as.matrix(allometry_phylo_group_cat_regline_coeffs_braincase[c(2, length(group_cat_vars)+2:(length(group_cat_vars))),])

#Calculate real intercepts and slopes
allometry_phylo_group_cat_regline_intercepts_ok_braincase <- as.matrix(c(allometry_phylo_group_cat_regline_intercepts_braincase[1,], allometry_phylo_group_cat_regline_intercepts_braincase[1,]+
                                                                     allometry_phylo_group_cat_regline_intercepts_braincase[2:length(allometry_phylo_group_cat_regline_intercepts_braincase),]))

allometry_phylo_group_cat_regline_slopes_ok_braincase <- as.matrix(c(allometry_phylo_group_cat_regline_slopes_braincase[1,], allometry_phylo_group_cat_regline_slopes_braincase[1,]+
                                                                 allometry_phylo_group_cat_regline_slopes_braincase[2:length(allometry_phylo_group_cat_regline_slopes_braincase),]))

#Save as data frame with grouping variables
allometry_phylo_group_cat_coeffs_braincase <- data.frame(Slope = allometry_phylo_group_cat_regline_slopes_ok_braincase, Intercept = allometry_phylo_group_cat_regline_intercepts_ok_braincase, 
                                                   row.names = group_cat_vars )
#Check for NA and other issues
allometry_phylo_group_cat_coeffs_braincase 

#Add classifiers
allometry_phylo_group_cat_coeffs_braincase <- allometry_phylo_group_cat_coeffs_braincase %>% mutate(category = rep(categories_list, each = 2), 
                                                                                        group = rep(str_to_title(groups_list), times = 4),
                                                                                        module = 'Braincase')
allometry_phylo_group_cat_coeffs_braincase

#Rostrum
#Linear model for line by group and category
allometry_phylo_group_cat_regline_rostrum <- lm(RegScores ~ logCS * group_cat, data = allometry_phylo_group_cat_module_2_plot_rostrum)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_group_cat_regline_coeffs_rostrum <- as.matrix(allometry_phylo_group_cat_regline_rostrum$coefficients)

#Save intercepts and slopes separately
group_cat_vars <- levels(allometry_phylo_group_cat_module_2_plot$group_cat)

allometry_phylo_group_cat_regline_intercepts_rostrum <- as.matrix(allometry_phylo_group_cat_regline_coeffs_rostrum[c(1, 3:(length(group_cat_vars)+1)),])
allometry_phylo_group_cat_regline_slopes_rostrum <- as.matrix(allometry_phylo_group_cat_regline_coeffs_rostrum[c(2, length(group_cat_vars)+2:(length(group_cat_vars))),])

#Calculate real intercepts and slopes
allometry_phylo_group_cat_regline_intercepts_ok_rostrum <- as.matrix(c(allometry_phylo_group_cat_regline_intercepts_rostrum[1,], allometry_phylo_group_cat_regline_intercepts_rostrum[1,]+
                                                                   allometry_phylo_group_cat_regline_intercepts_rostrum[2:length(allometry_phylo_group_cat_regline_intercepts_rostrum),]))

allometry_phylo_group_cat_regline_slopes_ok_rostrum <- as.matrix(c(allometry_phylo_group_cat_regline_slopes_rostrum[1,], allometry_phylo_group_cat_regline_slopes_rostrum[1,]+
                                                               allometry_phylo_group_cat_regline_slopes_rostrum[2:length(allometry_phylo_group_cat_regline_slopes_rostrum),]))

#Save as data frame with grouping variables
allometry_phylo_group_cat_coeffs_rostrum <- data.frame(Slope = allometry_phylo_group_cat_regline_slopes_ok_rostrum, Intercept = allometry_phylo_group_cat_regline_intercepts_ok_rostrum, 
                                                 row.names = group_cat_vars )
#Check for NA and other issues
allometry_phylo_group_cat_coeffs_rostrum 

#Add classifiers
allometry_phylo_group_cat_coeffs_rostrum <- allometry_phylo_group_cat_coeffs_rostrum %>% mutate(category = rep(categories_list, each = 2), 
                                                                                    group = rep(str_to_title(groups_list), times = 4),
                                                                                    module = 'Rostrum')
allometry_phylo_group_cat_coeffs_rostrum

#Put together
allometry_phylo_group_cat_coeffs <- rbind(allometry_phylo_group_cat_coeffs_braincase, allometry_phylo_group_cat_coeffs_rostrum)
allometry_phylo_group_cat_coeffs

#Make sure group correspond to use facet_wrap
allometry_phylo_braincase_grp_cat_plot$group <- str_to_title(allometry_phylo_braincase_grp_cat_plot$group)
allometry_phylo_rostrum_grp_cat_plot$group <- str_to_title(allometry_phylo_rostrum_grp_cat_plot$group)

#Divide by groups - too many lines all together
allometry_phylo_group_coeffs_mysticeti <- allometry_phylo_group_cat_coeffs %>% filter(group == "Mysticeti")
allometry_phylo_group_coeffs_odontoceti <- allometry_phylo_group_cat_coeffs %>% filter(group == "Odontoceti")

#Plot allometry regression by category mysticeti
allometry_phylo_group_module_2_ggplot_line_mysticeti <- ggplot(allometry_phylo_group_cat_module_2_plot_mysticeti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_abline(data = allometry_phylo_group_coeffs_mysticeti, 
              aes(intercept = Intercept, slope = Slope,  colour = module, linetype = category), linewidth = 1.2)+
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(shape = 21, linetype = 0, alpha =1, size = 4)), linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_group_module_2_ggplot_line_mysticeti

#Add phylopic
allometry_phylo_group_module_2_ggplot_line_mysticeti <- 
  allometry_phylo_group_module_2_ggplot_line_mysticeti +
  add_phylopic(myst, alpha = 1, x = 3, y = 1, ysize = 0.1, fill = "gray30")
allometry_phylo_group_module_2_ggplot_line_mysticeti
#Fix silhouette size

#Plot allometry regression by category odontoceti
allometry_phylo_group_module_2_ggplot_line_odontoceti <- ggplot(allometry_phylo_group_cat_module_2_plot_odontoceti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_abline(data = allometry_phylo_group_coeffs_odontoceti, 
              aes(intercept = Intercept, slope = Slope,  colour = module, linetype = category), linewidth = 1.2)+
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(override.aes = list(shape = 21, linetype = 0, alpha =1, size = 4)), linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_group_module_2_ggplot_line_odontoceti

#Add phylopic
allometry_phylo_group_module_2_ggplot_line_odontoceti <- 
  allometry_phylo_group_module_2_ggplot_line_odontoceti +
  add_phylopic(odont, alpha = 1, x = 2.75, y = 0.75, ysize = 0.1, fill = "gray50")
allometry_phylo_group_module_2_ggplot_line_odontoceti
#Fix silhouette size

ggarrange(allometry_phylo_group_module_2_ggplot_line_mysticeti,allometry_phylo_group_module_2_ggplot_line_odontoceti, 
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

####By module and category by group ----
##Plot by module both groups
allometry_phylo_braincase_grp_cat_line_ggplot <- ggplot(allometry_phylo_braincase_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_group_cat_coeffs[c(1:8),], 
              aes(intercept = Intercept, slope = Slope,  colour = category, linetype = group), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  facet_wrap(vars(group))+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  ggtitle("Braincase")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_braincase_grp_cat_line_ggplot

allometry_phylo_rostrum_grp_cat_line_ggplot <- ggplot(allometry_phylo_rostrum_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_group_cat_coeffs[c(9:16),], 
              aes(intercept = Intercept, slope = Slope,  colour = category, linetype = group), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  facet_wrap(vars(group))+
  theme_bw(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  ggtitle("Rostrum")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_rostrum_grp_cat_line_ggplot

ggarrange(allometry_phylo_rostrum_grp_cat_line_ggplot, allometry_phylo_braincase_grp_cat_line_ggplot,
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")


