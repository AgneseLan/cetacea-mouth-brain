
#===========================================================#
#                                                           #
#     SKULL MODULARITY - MYSTICETI & ODONTOCETI             #
#                                                           #
#===========================================================#

#CH. 6 - Disparity analyses whole, skull rostrum and braincase and rostrum vs braincase

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

#apropos("x") lists objects with matching part of name

 
#MORPHOLOGICAL DISPARITY WHOLE SKULL, ROSTRUM AND BRAINCASE ----
#Check if whole skull should be considered

#Create data frame with pc scores
##Use pcscores instead of coords - all same length of variables

#Make sure same columns and same order in all 3 tibbles
sort(names(pcscores_all_df))
sort(names(pcscores_rostrum_df))
sort(names(pcscores_braincase_df))

#Make common data frame
pcscores_R_B_W <- bind_rows(pcscores_all_df, pcscores_rostrum_df, pcscores_braincase_df)
glimpse(pcscores_R_B_W)

#Make column with interaction groups and category
pcscores_R_B_W$group_cat <- interaction(pcscores_R_B_W$group, pcscores_R_B_W$category)

#Create 2D matrix pc scores
pcscores_R_B_W_matrix <- as.matrix(pcscores_R_B_W[,c(1:194)])

#Eliminate repeated rownames
names(pcscores_R_B_W$size) <- NULL

#Create objects for disparity formula
modules_disparity<-pcscores_R_B_W$module
group_disparity <- pcscores_R_B_W$group

#Disparity between modules overall
disparity_modules <- morphol.disparity(pcscores_R_B_W_matrix ~ 1, groups= ~modules_disparity, 
                                       iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_modules)

#Save results to file
sink("Output/disparity_modules.txt")
print(summary(disparity_modules))
sink() 

#Disparity between modules in each group 
disparity_group_modules <- morphol.disparity(pcscores_R_B_W_matrix ~ 1, groups= ~group_disparity*modules_disparity, 
                                             iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_group_modules)

#Save results to file
sink("Output/disparity_group_modules.txt")
print(summary(disparity_group_modules))
sink() 

##Heatmaps plots for significant differences in disparity ----
#Create palette for heatmap trajectory plot
mypalette_disp <- brewer.pal(9,"PuRd")
image(1:9,1, as.matrix(1:9), col = mypalette_disp,xlab="PuRd (sequential)",
      ylab = "", yaxt = "n")

###By module ----
#Save p-values as object
disp_modules_corr <- disparity_modules[["PV.dist"]]
disp_modules_pvals <- disparity_modules[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_modules_vars <- str_to_title(rownames(disp_modules_corr))
disparity_modules_vars

#Set correct row and col names for both
rownames(disp_modules_corr) <- disparity_modules_vars
rownames(disp_modules_pvals) <- disparity_modules_vars
colnames(disp_modules_corr) <- disparity_modules_vars
colnames(disp_modules_pvals) <- disparity_modules_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_modules_corr_upper_tri <- get_upper_tri(disp_modules_corr)
disp_modules_pvals_upper_tri <- get_upper_tri(disp_modules_pvals)

#Melt to make table in the format needed for heatmap
disp_modules_corr_melt <- melt(disp_modules_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_modules_pvals_melt <- melt(disp_modules_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_modules_pvals_melt$corr <- disp_modules_corr_melt$corr

#Create columns where only significant values are shown
disp_modules_pvals_melt <- disp_modules_pvals_melt %>% mutate(sig_p = ifelse(p <= .05, T, F),
                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                              corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
disp_modules_pvals_melt

#Nice heatmap plot
disparity_modules_heatmap_ggplot <- ggplot(data = disp_modules_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.05), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values <= 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by module")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.5, 0.8),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_modules_heatmap_ggplot

###Sign difference between rostrum and skull, marginally sign rostrum and braincase
###Braincase and skull results similar vs rostrum

###By group and module ----
#Save p-values as object
disp_group_modules_corr <- disparity_group_modules[["PV.dist"]]
disp_group_modules_pvals <- disparity_group_modules[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_group_modules_vars <- rownames(disp_group_modules_corr)

#Replace string names to make them shorter
disparity_group_modules_vars <- str_replace_all(disparity_group_modules_vars, "\\.", "_")

#Make lists for labels
groups_list <- str_to_lower(levels(groups))
groups_list_short <- c("Myst", "Odont")

modules_list <- levels(as.factor(modules_disparity))
modules_list_short <- c("B","R","S")

#Loop replacements groups
for (t in 1:length(groups_list)){
  disparity_group_modules_vars <- str_replace_all(disparity_group_modules_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_list)){
  disparity_group_modules_vars <- str_replace_all(disparity_group_modules_vars, modules_list[m], modules_list_short[m])
}

#Check it worked
disparity_group_modules_vars

#Set correct row and col names for both
rownames(disp_group_modules_corr) <- disparity_group_modules_vars
rownames(disp_group_modules_pvals) <- disparity_group_modules_vars
colnames(disp_group_modules_corr) <- disparity_group_modules_vars
colnames(disp_group_modules_pvals) <- disparity_group_modules_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_group_modules_corr_upper_tri <- get_upper_tri(disp_group_modules_corr)
disp_group_modules_pvals_upper_tri <- get_upper_tri(disp_group_modules_pvals)

#Melt to make table in the format needed for heatmap
disp_group_modules_corr_melt <- melt(disp_group_modules_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_group_modules_pvals_melt <- melt(disp_group_modules_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_group_modules_pvals_melt$corr <- disp_group_modules_corr_melt$corr

#Create columns where only significant values are shown
disp_group_modules_pvals_melt <- disp_group_modules_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                              corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
disp_group_modules_pvals_melt

#Nice heatmap plot
disparity_group_modules_heatmap_ggplot <- ggplot(data = disp_group_modules_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by module and group")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_group_modules_heatmap_ggplot

###Strong diff in disparity between groups across modules
###In Odontoceti, no diff between modules
###In Mysticeti, diff between braincase and rostrum and rostrum and skull

##Lack of pairwise disparity between whole skull and both rostrum and braincase
##Whole skull not needed to be considered in analysis

#MORPHOLOGICAL DISPARITY ROSTRUM AND BRAINCASE SEPARATE----

##Rostrum ----
#Disparity between modules in each group and growth stage
disparity_rostrum_group_cat <- morphol.disparity(coords ~ 1, groups= ~group*category, data = gdf_rostrum,
                                                 iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_rostrum_group_cat)

#Save results to file
sink("Output/disparity_rostrum_group_cat.txt")
print(summary(disparity_rostrum_group_cat))
sink() 

###Heatmaps plots for significant differences in disparity ----

#Create palette for comparison between modules
mypalette_seq_purple <-as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Purple"]][["value"]])
image(1:20,1, as.matrix(1:20), col = mypalette_seq_purple ,xlab="Purples (sequential)",
      ylab = "", yaxt = "n")
mypalette_seq_purple[1,] <- "#f9eef7" #change to make first color lighter
mypalette_seq_purple <- as.vector(mypalette_seq_purple)

mypalette_seq_modules <- colorRampPalette(c(mypalette_seq_purple[1], mypalette_seq_purple[10], mypalette_seq_purple[20]))
plot(rep(1,10),col=mypalette_seq_modules(10),pch=19,cex=3)
mypalette_seq_modules <- mypalette_seq_modules(10)

#Save p-values as object
disp_rostrum_corr <- disparity_rostrum_group_cat[["PV.dist"]]
disp_rostrum_pvals <- disparity_rostrum_group_cat[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_rostrum_vars <- rownames(disp_rostrum_corr)

#Replace string names to make them shorter
disparity_rostrum_vars <- str_replace_all(disparity_rostrum_vars, "\\.", "_")

#Make list for labels
categories_list <- levels(gdf$category)
categories_list_short <- c("1","2","3","4")

#Loop replacements categories
for (u in 1:length(categories_list)){
  disparity_rostrum_vars <- str_replace_all(disparity_rostrum_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  disparity_rostrum_vars <- str_replace_all(disparity_rostrum_vars, groups_list[t], groups_list_short[t])
}

#Check it worked
disparity_rostrum_vars

#Set correct row and col names for both
rownames(disp_rostrum_corr) <- disparity_rostrum_vars
rownames(disp_rostrum_pvals) <- disparity_rostrum_vars
colnames(disp_rostrum_corr) <- disparity_rostrum_vars
colnames(disp_rostrum_pvals) <- disparity_rostrum_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_rostrum_corr_upper_tri <- get_upper_tri(disp_rostrum_corr)
disp_rostrum_pvals_upper_tri <- get_upper_tri(disp_rostrum_pvals)

#Melt to make table in the format needed for heatmap
disp_rostrum_corr_melt <- melt(disp_rostrum_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_rostrum_pvals_melt <- melt(disp_rostrum_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_rostrum_pvals_melt$corr <- disp_rostrum_corr_melt$corr

#Create columns where only significant values are shown
disp_rostrum_pvals_melt <- disp_rostrum_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                              corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))

disp_rostrum_pvals_melt

#Nice heatmap plot
disparity_rostrum_group_cat_heatmap_ggplot <- ggplot(data = disp_rostrum_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (label = "Rostrum")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14), plot.subtitle = element_text(face = 3, hjust = 0.5, size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_rostrum_group_cat_heatmap_ggplot

##Braincase ----
#Disparity between modules in each group and growth stage
disparity_braincase_group_cat <- morphol.disparity(coords ~ 1, groups= ~group*category, data = gdf_braincase,
                                                   iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_braincase_group_cat)

#Save results to file
sink("Output/disparity_braincase_group_cat.txt")
print(summary(disparity_braincase_group_cat))
sink() 

###Heatmaps plots for significant differences in disparity ----
#Save p-values as object
disp_braincase_corr <- disparity_braincase_group_cat[["PV.dist"]]
disp_braincase_pvals <- disparity_braincase_group_cat[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_braincase_vars <- rownames(disp_braincase_corr)

#Replace string names to make them shorter
disparity_braincase_vars <- str_replace_all(disparity_braincase_vars, "\\.", "_")

#Loop replacements categories
for (u in 1:length(categories_list)){
  disparity_braincase_vars <- str_replace_all(disparity_braincase_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  disparity_braincase_vars <- str_replace_all(disparity_braincase_vars, groups_list[t], groups_list_short[t])
}

#Check it worked
disparity_braincase_vars

#Set correct row and col names for both
rownames(disp_braincase_corr) <- disparity_braincase_vars
rownames(disp_braincase_pvals) <- disparity_braincase_vars
colnames(disp_braincase_corr) <- disparity_braincase_vars
colnames(disp_braincase_pvals) <- disparity_braincase_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_braincase_corr_upper_tri <- get_upper_tri(disp_braincase_corr)
disp_braincase_pvals_upper_tri <- get_upper_tri(disp_braincase_pvals)

#Melt to make table in the format needed for heatmap
disp_braincase_corr_melt <- melt(disp_braincase_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_braincase_pvals_melt <- melt(disp_braincase_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_braincase_pvals_melt$corr <- disp_braincase_corr_melt$corr

#Create columns where only significant values are shown
disp_braincase_pvals_melt <- disp_braincase_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                  corr_if_sig = ifelse(sig_p, corr, NA)) %>%
mutate_at(vars(starts_with("corr")), list(~ round(., 3)))

disp_braincase_pvals_melt

#Nice heatmap plot
disparity_braincase_group_cat_heatmap_ggplot <- ggplot(data = disp_braincase_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_modules[9], high = mypalette_seq_modules[2], mid = mypalette_seq_modules[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modules[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (label = "Braincase")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14), plot.subtitle = element_text(face = 3, hjust = 0.5, size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
disparity_braincase_group_cat_heatmap_ggplot

plotD1<-ggarrange(disparity_rostrum_group_cat_heatmap_ggplot,disparity_braincase_group_cat_heatmap_ggplot,
          nrow = 1, ncol = 2, common.legend = F)
plotD1<-annotate_figure(plotD1, top = text_grob("Morphological disparity between group and growth stage by module", face = "bold", size = 16))
plotD1

#MORPHOLOGICAL DISPARITY ROSTRUM AND BRAINCASE COMPARE ----
#Create data frame with pc scores
##Use pcscores instead of coords - all same length of variables

#Make common data frame
pcscores_R_B <- bind_rows(pcscores_rostrum_df, pcscores_braincase_df)
glimpse(pcscores_R_B)

#Make column with interaction groups and category
pcscores_R_B$group_cat <- interaction(pcscores_R_B$group, pcscores_R_B$category)

#Create 2D matrix pc scores
pcscores_R_B_matrix <- as.matrix(pcscores_R_B[,c(1:194)])

#Eliminate repeated rownames
names(pcscores_R_B$size) <- NULL

##Modules by group
#Create objects for disparity formula
group_cat_2_disparity <- pcscores_R_B$group_cat
modules_2_disparity <- pcscores_R_B$module
group_2_disparity <- pcscores_R_B$group

#Disparity between modules in each group
disparity_group_modules_2 <- morphol.disparity(pcscores_R_B_matrix ~ 1, groups= ~group_2_disparity*modules_2_disparity, 
                                                   iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_group_modules_2)

#Save results to file
sink("Output/disparity_group_modules_2.txt")
print(summary(disparity_group_modules_2))
sink() 

##Modules by group and category
#Disparity between modules in each group and growth stage
disparity_group_cat_modules_2 <- morphol.disparity(pcscores_R_B_matrix ~ 1, groups= ~group_cat_2_disparity*modules_2_disparity, 
                                                   iter = 999, print.progress = FALSE)

#Results and significance
summary(disparity_group_cat_modules_2)

#Save results to file
sink("Output/disparity_group_cat_modules_2.txt")
print(summary(disparity_group_cat_modules_2))
sink() 

###Heatmaps plots for significant differences in disparity ----
#Create palette for comparison between modules and groups
mypalette_seq_mod_grp <- brewer.pal(9,"Reds")
image(1:9,1, as.matrix(1:9), col = mypalette_seq_mod_grp,xlab="Reds (sequential)",
      ylab = "", yaxt = "n")

#Create palette for comparison between groups
mypalette_seq_groups <- brewer.pal(9,"YlGn")
image(1:9,1, as.matrix(1:9), col = mypalette_seq_groups,xlab="YlGn (sequential)",
      ylab = "", yaxt = "n")

####Modules by group ----
#Save p-values as object
disp_modules_2_corr <- disparity_group_modules_2[["PV.dist"]]
disp_modules_2_pvals <- disparity_group_modules_2[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_modules_2_vars <- rownames(disp_modules_2_corr)

#Replace string names to make them shorter
disparity_modules_2_vars <- str_replace_all(disparity_modules_2_vars, "\\.", "_")

#List 2 modules
modules_2_list <- c("braincase","rostrum")
modules_2_list_short <- c("B","R")

#Loop replacements groups
for (t in 1:length(groups_list)){
  disparity_modules_2_vars <- str_replace_all(disparity_modules_2_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  disparity_modules_2_vars <- str_replace_all(disparity_modules_2_vars, modules_2_list[m], modules_2_list_short[m])
}

#Check it worked
disparity_modules_2_vars

#Set correct row and col names for both
rownames(disp_modules_2_corr) <- disparity_modules_2_vars
rownames(disp_modules_2_pvals) <- disparity_modules_2_vars
colnames(disp_modules_2_corr) <- disparity_modules_2_vars
colnames(disp_modules_2_pvals) <- disparity_modules_2_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_modules_2_corr_upper_tri <- get_upper_tri(disp_modules_2_corr)
disp_modules_2_pvals_upper_tri <- get_upper_tri(disp_modules_2_pvals)

#Melt to make table in the format needed for heatmap
disp_modules_2_corr_melt <- melt(disp_modules_2_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_modules_2_pvals_melt <- melt(disp_modules_2_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_modules_2_pvals_melt$corr <- disp_modules_2_corr_melt$corr

#Create columns where only significant values are shown
disp_modules_2_pvals_melt <- disp_modules_2_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          corr_if_sig = ifelse(sig_p, corr, NA)) %>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
disp_modules_2_pvals_melt

#Nice heatmap plot
disparity_group_modules_2_heatmap_ggplot <- ggplot(data = disp_modules_2_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[1], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by module and group")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_group_modules_2_heatmap_ggplot


####Modules by group and category ----
#Save p-values as object
disp_modules_2_cat_corr <- disparity_group_cat_modules_2[["PV.dist"]]
disp_modules_2_cat_pvals <- disparity_group_cat_modules_2[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
disparity_modules_2_cat_vars <- rownames(disp_modules_2_cat_corr)

#Replace string names to make them shorter
disparity_modules_2_cat_vars <- str_replace_all(disparity_modules_2_cat_vars, "\\.", "_")

#List 2 modules
modules_2_list <- c("braincase","rostrum")
modules_2_list_short <- c("B","R")

#Loop replacements categories
for (u in 1:length(categories_list)){
  disparity_modules_2_cat_vars <- str_replace_all(disparity_modules_2_cat_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  disparity_modules_2_cat_vars <- str_replace_all(disparity_modules_2_cat_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  disparity_modules_2_cat_vars <- str_replace_all(disparity_modules_2_cat_vars, modules_2_list[m], modules_2_list_short[m])
}

#Check it worked
disparity_modules_2_cat_vars

#Set correct row and col names for both
rownames(disp_modules_2_cat_corr) <- disparity_modules_2_cat_vars
rownames(disp_modules_2_cat_pvals) <- disparity_modules_2_cat_vars
colnames(disp_modules_2_cat_corr) <- disparity_modules_2_cat_vars
colnames(disp_modules_2_cat_pvals) <- disparity_modules_2_cat_vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_modules_2_cat_corr_upper_tri <- get_upper_tri(disp_modules_2_cat_corr)
disp_modules_2_cat_pvals_upper_tri <- get_upper_tri(disp_modules_2_cat_pvals)

#Melt to make table in the format needed for heatmap
disp_modules_2_cat_corr_melt <- melt(disp_modules_2_cat_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_modules_2_cat_pvals_melt <- melt(disp_modules_2_cat_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_modules_2_cat_pvals_melt$corr <- disp_modules_2_cat_corr_melt$corr

#Create columns where only significant values are shown
disp_modules_2_cat_pvals_melt <- disp_modules_2_cat_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          corr_if_sig = ifelse(sig_p, corr, NA)) %>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
disp_modules_2_cat_pvals_melt

#Nice heatmap plot
disparity_group_cat_modules_2_heatmap_ggplot <- ggplot(data = disp_modules_2_cat_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[1], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by module, group and growth stage")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_group_cat_modules_2_heatmap_ggplot

#Big plot, divide by group
#Make data frame for each group
disp_modules_2_cat_pvals_melt_mysticeti <- disp_modules_2_cat_pvals_melt %>% filter(str_detect(Var1, "Myst")) %>% 
  filter(str_detect(Var2, "Myst"))
disp_modules_2_cat_pvals_melt_odontoceti <- disp_modules_2_cat_pvals_melt %>% filter(str_detect(Var1, "Odont")) %>% 
  filter(str_detect(Var2, "Odont"))

#Create labels
category_modules_2_list <- as.factor(paste(rep(modules_2_list_short, times = length(categories_list_short)), 
                                           rep(categories_list_short, each = length(modules_2_list_short)), sep = "_"))

#Nice heatmap plot
disparity_group_cat_modules_2_heatmap_ggplot_myst <- ggplot(data = disp_modules_2_cat_pvals_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  scale_x_discrete(labels = category_modules_2_list)+
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Mysticeti")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_group_cat_modules_2_heatmap_ggplot_myst

disparity_group_cat_modules_2_heatmap_ggplot_odont <- ggplot(data = disp_modules_2_cat_pvals_melt_odontoceti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_groups[9], high = mypalette_seq_groups[2], mid = mypalette_seq_groups[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_groups[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  scale_x_discrete(labels = category_modules_2_list)+
  scale_y_discrete(labels = category_modules_2_list)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Odontoceti")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
disparity_group_cat_modules_2_heatmap_ggplot_odont

plotD2 <- ggarrange(disparity_group_cat_modules_2_heatmap_ggplot_myst, 
                  disparity_group_cat_modules_2_heatmap_ggplot_odont,
                  ncol = 2, nrow = 1, common.legend = F)
plotD2<-annotate_figure(plotD2, top = text_grob("Morphological disparity between module and growth stage by group", 
                                      face = "bold", size = 16, just = c(0.5,1)))
plotD2

###Line plots for Procrustes variance of groups at each stage ----
##All data
#Save variances as object
PV_group_cat_modules_2 <- data.frame(PV = disparity_group_cat_modules_2[["Procrustes.var"]])

#Look at order of variables to create vectors
PV_group_cat_modules_2

#Make vectors with correct order an number of categories and feeding
PV_modules_2 <- rep(modules_2_list, times = length(groups_list)*length(categories_list))

PV_categories_2 <- rep(categories_list, each = length(modules_2_list), times = length(groups_list))

PV_groups_2 <- rep(groups_list, each = length(PV_modules_2)/length(groups_list))

#Add labels and other attributes to tibble as columns
PV_group_cat_modules_2 <- PV_group_cat_modules_2 %>% 
  mutate(module = str_to_sentence(PV_modules_2), category= PV_categories_2, group = PV_groups_2)
PV_group_cat_modules_2

#Nice line plot by module
PV_group_cat_modules_2_ggplot <- ggplot(PV_group_cat_modules_2, aes(x = category, y = PV)) + 
  geom_line(aes(x = category, y = PV,linetype = group, colour = module, group = group), linewidth = 1, inherit.aes = F)+
  geom_point(aes(color = module,
                 shape = group),size = 4, fill = "white", stroke = 1.5)+
  facet_wrap(vars(module))+
  scale_shape_manual(name = "Groups", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Modules", labels =levels(as.factor(PV_group_cat_modules_2$module)),
                      values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour","fill"))+   
  scale_linetype_manual(name = "Groups", labels = levels(groups),
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("PV (Procrustes variances)")+
  xlab("Growth stage")+
  ggtitle ("Morphological disparity per module, group and growth stage")+ 
  scale_x_discrete(labels = levels(as.factor(categories_list_short)))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", linetype = "none", shape = "none")
PV_group_cat_modules_2_ggplot

#Add significance for each segment
PV_group_cat_modules_2_ggplot  <- 
  PV_group_cat_modules_2_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 1, y = 0.09, ysize = 0.004, color = "gray30")+
  add_phylopic(odont, alpha = 1, x = 1, y = 0.04, ysize = 0.004, color = "gray50")
PV_group_cat_modules_2_ggplot

#Nice line plot by group
group_labels <- as_labeller(c('mysticeti' = "Mysticeti", 'odontoceti' = "Odontoceti"))

PV_module_cat_modules_2_ggplot <- ggplot(PV_group_cat_modules_2, aes(x = category, y = PV)) + 
  geom_line(aes(x = category, y = PV,linetype = module, colour = group, group = module), linewidth = 1, inherit.aes = F)+
  geom_point(aes(color = group,
                 shape = module),size = 4, fill = "white", stroke = 1.5)+
  facet_wrap(vars(group), labeller = group_labels)+
  scale_shape_manual(name = "Modules", labels = levels(as.factor(PV_group_cat_modules_2$module)), values = c(23,24))+
  scale_colour_manual(name = "Groups", labels =levels(groups),
                      values = mypalette_groups, aesthetics = c("colour","fill"))+   
  scale_linetype_manual(name = "Modules", labels = levels(as.factor(PV_group_cat_modules_2$module)),
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("PV (Procrustes variances)")+
  xlab("Growth stage")+
  ggtitle ("Morphological disparity per group, module and growth stage")+ 
  scale_x_discrete(labels = levels(as.factor(categories_list_short)))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", shape = guide_legend(keywidth=unit(3, "char"),override.aes = list(size = 3)))
PV_module_cat_modules_2_ggplot

#Add silhouettes groups
PV_module_cat_modules_2_ggplot  <- 
  PV_module_cat_modules_2_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 1, y = 0.09, ysize = 0.004, color = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 1, y = 0.04, ysize = 0.004, color = mypalette_groups[2])
PV_module_cat_modules_2_ggplot
#Delete extra silhouettes and change position

###### 

#Next - ch. 7 - Trajectory analyses