#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH. 6b - Disparity analyses whole, skull rostrum and braincase and rostrum vs braincase - phylogenetically transformed components mean shapes

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
options(scipen = 30)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#remotes::install_github("r-lib/rray")

#apropos("x") lists objects with matching part of name


#MORPHOLOGICAL DISPARITY WHOLE SKULL, ROSTRUM AND BRAINCASE ----
#Check if whole skull should be considered

#Create data frame with pc scores
##Use phylo_pcscores instead of coords - all same length of variables

#Make sure same columns and same order in all 3 tibbles
sort(names(phylo_pcscores_all_means_df))
sort(names(phylo_pcscores_rostrum_means_df))
sort(names(phylo_pcscores_braincase_means_df))

#Make common data frame
phylo_pcscores_means_R_B_W <- bind_rows(phylo_pcscores_all_means_df, phylo_pcscores_rostrum_means_df, phylo_pcscores_braincase_means_df)
glimpse(phylo_pcscores_means_R_B_W)

#Make column with interaction groups and category
phylo_pcscores_means_R_B_W$group_cat <- interaction(phylo_pcscores_means_R_B_W$group, phylo_pcscores_means_R_B_W$category)

#Create 2D matrix pc scores
phylo_pcscores_means_R_B_W_matrix <- as.matrix(phylo_pcscores_means_R_B_W[,c(1:14)]) #copy number of components

#Eliminate repeated rownames
names(phylo_pcscores_means_R_B_W$size) <- NULL

#Create objects for disparity formula
modules_phylo_disparity<-phylo_pcscores_means_R_B_W$module
group_phylo_disparity <- phylo_pcscores_means_R_B_W$group

#Disparity between modules overall
phylo_disparity_means_modules <- morphol.disparity(phylo_pcscores_means_R_B_W_matrix ~ 1, groups= ~modules_phylo_disparity, 
                                             iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_modules)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_modules.txt")
print(summary(phylo_disparity_means_modules))
sink() 

#Disparity between modules in each group 
phylo_disparity_means_group_modules <- morphol.disparity(phylo_pcscores_means_R_B_W_matrix ~ 1, groups= ~group_phylo_disparity*modules_phylo_disparity, 
                                                   iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_group_modules)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_group_modules.txt")
print(summary(phylo_disparity_means_group_modules))
sink() 

##Heatmaps plots for significant differences in phylo_disparity ----

#Create palette for heatmap trajectory plot
mypalette_disp <- brewer.pal(9,"PuRd")
image(1:9,1, as.matrix(1:9), col = mypalette_disp,xlab="PuRd (sequential)",
      ylab = "", yaxt = "n")

###By module ----
#Save p-values as object
phylo_disp_means_modules_corr <- phylo_disparity_means_modules[["PV.dist"]]
phylo_disp_means_modules_pvals <- phylo_disparity_means_modules[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_modules_vars <- str_to_title(rownames(phylo_disp_means_modules_corr))
phylo_disparity_means_modules_vars

#Set correct row and col names for both
rownames(phylo_disp_means_modules_corr) <- phylo_disparity_means_modules_vars
rownames(phylo_disp_means_modules_pvals) <- phylo_disparity_means_modules_vars
colnames(phylo_disp_means_modules_corr) <- phylo_disparity_means_modules_vars
colnames(phylo_disp_means_modules_pvals) <- phylo_disparity_means_modules_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_modules_corr_upper_tri <- get_upper_tri(phylo_disp_means_modules_corr)
phylo_disp_means_modules_pvals_upper_tri <- get_upper_tri(phylo_disp_means_modules_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_modules_corr_melt <- melt(phylo_disp_means_modules_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_modules_pvals_melt <- melt(phylo_disp_means_modules_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_modules_pvals_melt$corr <- phylo_disp_means_modules_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_modules_pvals_melt <- phylo_disp_means_modules_pvals_melt %>% mutate(sig_p = ifelse(p <= .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 4)))
phylo_disp_means_modules_pvals_melt

#Nice heatmap plot
phylo_disparity_means_modules_heatmap_ggplot <- ggplot(data = phylo_disp_means_modules_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.05), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values <= 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity mean shapes phylogenetically corrected by module")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.8),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_modules_heatmap_ggplot

###Sign differences between all pairs

###By group and module ----
#Save p-values as object
phylo_disp_means_group_modules_corr <- phylo_disparity_means_group_modules[["PV.dist"]]
phylo_disp_means_group_modules_pvals <- phylo_disparity_means_group_modules[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_group_modules_vars <- rownames(phylo_disp_means_group_modules_corr)

#Replace string names to make them shorter
phylo_disparity_means_group_modules_vars <- str_replace_all(phylo_disparity_means_group_modules_vars, "\\.", "_")

#Make lists for labels
modules_list <- levels(as.factor(modules_phylo_disparity))
modules_list_short <- c("B","R","S")

#Loop replacements groups
for (t in 1:length(groups_list)){
  phylo_disparity_means_group_modules_vars <- str_replace_all(phylo_disparity_means_group_modules_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_list)){
  phylo_disparity_means_group_modules_vars <- str_replace_all(phylo_disparity_means_group_modules_vars, modules_list[m], modules_list_short[m])
}

#Check it worked
phylo_disparity_means_group_modules_vars

#Set correct row and col names for both
rownames(phylo_disp_means_group_modules_corr) <- phylo_disparity_means_group_modules_vars
rownames(phylo_disp_means_group_modules_pvals) <- phylo_disparity_means_group_modules_vars
colnames(phylo_disp_means_group_modules_corr) <- phylo_disparity_means_group_modules_vars
colnames(phylo_disp_means_group_modules_pvals) <- phylo_disparity_means_group_modules_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_group_modules_corr_upper_tri <- get_upper_tri(phylo_disp_means_group_modules_corr)
phylo_disp_means_group_modules_pvals_upper_tri <- get_upper_tri(phylo_disp_means_group_modules_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_group_modules_corr_melt <- melt(phylo_disp_means_group_modules_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_group_modules_pvals_melt <- melt(phylo_disp_means_group_modules_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_group_modules_pvals_melt$corr <- phylo_disp_means_group_modules_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_group_modules_pvals_melt <- phylo_disp_means_group_modules_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                      p_if_sig = ifelse(sig_p, p, NA),
                                                                                      corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 4)))
phylo_disp_means_group_modules_pvals_melt

#Nice heatmap plot
phylo_disparity_means_group_modules_heatmap_ggplot <- ggplot(data = phylo_disp_means_group_modules_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity mean shapes phylogenetically corrected by module and group")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_group_modules_heatmap_ggplot

###Difference in disparity between the braincase and other modules across groups but not between the Odont and Myst braincase
###In Odontoceti, diff between all modules
###In Mysticeti, diff between braincase and rostrum and braincase and skull

##Significant difference between braincase and rostrum
##Whole skull not needed to be considered in analysis

#MORPHOLOGICAL DISPARITY ROSTRUM AND BRAINCASE SEPARATE----

##Rostrum ----
#Create 2D matrix pc scores
phylo_pcscores_rostrum_means_matrix <- as.matrix(phylo_pcscores_rostrum_means_df[,1:14]) #copy number of components

#Create interaction for disparity formula
group_cat_phylo_disparity <- interaction(phylo_pcscores_rostrum_means_df$group,phylo_pcscores_rostrum_means_df$category)

#Disparity between modules in each group and growth stage
phylo_disparity_means_rostrum_group_cat <- morphol.disparity(phylo_pcscores_rostrum_means_matrix ~ 1, groups= ~group_cat_phylo_disparity, 
                                                       iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_rostrum_group_cat)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_rostrum_group_cat.txt")
print(summary(phylo_disparity_means_rostrum_group_cat))
sink() 

###Heatmaps plots for significant differences in phylo_disparity ----

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
phylo_disp_means_rostrum_corr <- phylo_disparity_means_rostrum_group_cat[["PV.dist"]]
phylo_disp_means_rostrum_pvals <- phylo_disparity_means_rostrum_group_cat[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_rostrum_vars <- rownames(phylo_disp_means_rostrum_corr)

#Replace string names to make them shorter
phylo_disparity_means_rostrum_vars <- str_replace_all(phylo_disparity_means_rostrum_vars, "\\.", "_")

#Loop replacements categories
for (u in 1:length(categories_list)){
  phylo_disparity_means_rostrum_vars <- str_replace_all(phylo_disparity_means_rostrum_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  phylo_disparity_means_rostrum_vars <- str_replace_all(phylo_disparity_means_rostrum_vars, groups_list[t], groups_list_short[t])
}

#Check it worked
phylo_disparity_means_rostrum_vars

#Set correct row and col names for both
rownames(phylo_disp_means_rostrum_corr) <- phylo_disparity_means_rostrum_vars
rownames(phylo_disp_means_rostrum_pvals) <- phylo_disparity_means_rostrum_vars
colnames(phylo_disp_means_rostrum_corr) <- phylo_disparity_means_rostrum_vars
colnames(phylo_disp_means_rostrum_pvals) <- phylo_disparity_means_rostrum_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_rostrum_corr_upper_tri <- get_upper_tri(phylo_disp_means_rostrum_corr)
phylo_disp_means_rostrum_pvals_upper_tri <- get_upper_tri(phylo_disp_means_rostrum_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_rostrum_corr_melt <- melt(phylo_disp_means_rostrum_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_rostrum_pvals_melt <- melt(phylo_disp_means_rostrum_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_rostrum_pvals_melt$corr <- phylo_disp_means_rostrum_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_rostrum_pvals_melt <- phylo_disp_means_rostrum_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          corr_if_sig = ifelse(sig_p, corr, NA))%>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))

phylo_disp_means_rostrum_pvals_melt

#Nice heatmap plot
phylo_disparity_means_rostrum_group_cat_heatmap_ggplot <- ggplot(data = phylo_disp_means_rostrum_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
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
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_rostrum_group_cat_heatmap_ggplot

##Braincase ----
#Create 2D matrix pc scores
phylo_pcscores_braincase_means_matrix <- as.matrix(phylo_pcscores_braincase_means_df[,1:14]) #copy number of components

#Create interaction for disparity formula
group_cat_phylo_disparity <- interaction(phylo_pcscores_braincase_means_df$group,phylo_pcscores_braincase_means_df$category)

#Disparity between modules in each group and growth stage
phylo_disparity_means_braincase_group_cat <- morphol.disparity(phylo_pcscores_braincase_means_matrix ~ 1, groups= ~group_cat_phylo_disparity, 
                                                         iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_braincase_group_cat)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_braincase_group_cat.txt")
print(summary(phylo_disparity_means_braincase_group_cat))
sink() 

###Heatmaps plots for significant differences in phylo_disparity ----
#Save p-values as object
phylo_disp_means_braincase_corr <- phylo_disparity_means_braincase_group_cat[["PV.dist"]]
phylo_disp_means_braincase_pvals <- phylo_disparity_means_braincase_group_cat[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_braincase_vars <- rownames(phylo_disp_means_braincase_corr)

#Replace string names to make them shorter
phylo_disparity_means_braincase_vars <- str_replace_all(phylo_disparity_means_braincase_vars, "\\.", "_")

#Loop replacements categories
for (u in 1:length(categories_list)){
  phylo_disparity_means_braincase_vars <- str_replace_all(phylo_disparity_means_braincase_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  phylo_disparity_means_braincase_vars <- str_replace_all(phylo_disparity_means_braincase_vars, groups_list[t], groups_list_short[t])
}

#Check it worked
phylo_disparity_means_braincase_vars

#Set correct row and col names for both
rownames(phylo_disp_means_braincase_corr) <- phylo_disparity_means_braincase_vars
rownames(phylo_disp_means_braincase_pvals) <- phylo_disparity_means_braincase_vars
colnames(phylo_disp_means_braincase_corr) <- phylo_disparity_means_braincase_vars
colnames(phylo_disp_means_braincase_pvals) <- phylo_disparity_means_braincase_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_braincase_corr_upper_tri <- get_upper_tri(phylo_disp_means_braincase_corr)
phylo_disp_means_braincase_pvals_upper_tri <- get_upper_tri(phylo_disp_means_braincase_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_braincase_corr_melt <- melt(phylo_disp_means_braincase_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_braincase_pvals_melt <- melt(phylo_disp_means_braincase_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_braincase_pvals_melt$corr <- phylo_disp_means_braincase_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_braincase_pvals_melt <- phylo_disp_means_braincase_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                              corr_if_sig = ifelse(sig_p, corr, NA)) %>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))

phylo_disp_means_braincase_pvals_melt

#Nice heatmap plot
phylo_disparity_means_braincase_group_cat_heatmap_ggplot <- ggplot(data = phylo_disp_means_braincase_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
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
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
phylo_disparity_means_braincase_group_cat_heatmap_ggplot

plotD1<-ggarrange(phylo_disparity_means_rostrum_group_cat_heatmap_ggplot,
                  nrow = 1, ncol = 1, common.legend = F)
plotD1<-annotate_figure(plotD1, top = text_grob("Morphological disparity means phylogenetically corrected between group and growth stage by module", face = "bold", size = 16))
plotD1

#MORPHOLOGICAL DISPARITY ROSTRUM AND BRAINCASE COMPARE ----
#Create data frame with pc scores
##Use phylo_pcscores instead of coords - all same length of variables

#Make common data frame
phylo_pcscores_means_R_B <- bind_rows(phylo_pcscores_rostrum_means_df, phylo_pcscores_braincase_means_df)
glimpse(phylo_pcscores_means_R_B)

#Make column with interaction groups and category
phylo_pcscores_means_R_B$group_cat <- interaction(phylo_pcscores_means_R_B$group, phylo_pcscores_means_R_B$category)

#Create 2D matrix pc scores
phylo_pcscores_means_R_B_matrix <- as.matrix(phylo_pcscores_means_R_B[,c(1:14)])

#Eliminate repeated rownames
names(phylo_pcscores_means_R_B$size) <- NULL

##Modules by group
#Create objects for phylo_disparity formula
group_cat_2_phylo_disparity <- phylo_pcscores_means_R_B$group_cat
modules_2_phylo_disparity <- phylo_pcscores_means_R_B$module
group_2_phylo_disparity <- phylo_pcscores_means_R_B$group

#Disparity between modules in each group
phylo_disparity_means_group_modules_2 <- morphol.disparity(phylo_pcscores_means_R_B_matrix ~ 1, groups= ~group_2_phylo_disparity*modules_2_phylo_disparity, 
                                                     iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_group_modules_2)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_group_modules_2.txt")
print(summary(phylo_disparity_means_group_modules_2))
sink() 

##Modules by group and category
#Disparity between modules in each group and growth stage
phylo_disparity_means_group_cat_modules_2 <- morphol.disparity(phylo_pcscores_means_R_B_matrix ~ 1, groups= ~group_cat_2_phylo_disparity*modules_2_phylo_disparity, 
                                                         iter = 999, print.progress = FALSE)

#Results and significance
summary(phylo_disparity_means_group_cat_modules_2)

#Save results to file
sink("Output/6b-Disparity phylo means/phylo_disparity_means_group_cat_modules_2.txt")
print(summary(phylo_disparity_means_group_cat_modules_2))
sink() 

###Heatmaps plots for significant differences in phylo_disparity ----
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
phylo_disp_means_modules_2_corr <- phylo_disparity_means_group_modules_2[["PV.dist"]]
phylo_disp_means_modules_2_pvals <- phylo_disparity_means_group_modules_2[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_modules_2_vars <- rownames(phylo_disp_means_modules_2_corr)

#Replace string names to make them shorter
phylo_disparity_means_modules_2_vars <- str_replace_all(phylo_disparity_means_modules_2_vars, "\\.", "_")

#List 2 modules
modules_2_list <- c("braincase","rostrum")
modules_2_list_short <- c("B","R")

#Loop replacements groups
for (t in 1:length(groups_list)){
  phylo_disparity_means_modules_2_vars <- str_replace_all(phylo_disparity_means_modules_2_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  phylo_disparity_means_modules_2_vars <- str_replace_all(phylo_disparity_means_modules_2_vars, modules_2_list[m], modules_2_list_short[m])
}

#Check it worked
phylo_disparity_means_modules_2_vars

#Set correct row and col names for both
rownames(phylo_disp_means_modules_2_corr) <- phylo_disparity_means_modules_2_vars
rownames(phylo_disp_means_modules_2_pvals) <- phylo_disparity_means_modules_2_vars
colnames(phylo_disp_means_modules_2_corr) <- phylo_disparity_means_modules_2_vars
colnames(phylo_disp_means_modules_2_pvals) <- phylo_disparity_means_modules_2_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_modules_2_corr_upper_tri <- get_upper_tri(phylo_disp_means_modules_2_corr)
phylo_disp_means_modules_2_pvals_upper_tri <- get_upper_tri(phylo_disp_means_modules_2_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_modules_2_corr_melt <- melt(phylo_disp_means_modules_2_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_modules_2_pvals_melt <- melt(phylo_disp_means_modules_2_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_modules_2_pvals_melt$corr <- phylo_disp_means_modules_2_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_modules_2_pvals_melt <- phylo_disp_means_modules_2_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                              p_if_sig = ifelse(sig_p, p, NA),
                                                                              corr_if_sig = ifelse(sig_p, corr, NA)) %>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
phylo_disp_means_modules_2_pvals_melt

#Nice heatmap plot
phylo_disparity_means_group_modules_2_heatmap_ggplot <- ggplot(data = phylo_disp_means_modules_2_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = corr_if_sig), color = "white", size = 5) +
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[1], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity phylogenetically corrected means by module and group")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_group_modules_2_heatmap_ggplot


####Modules by group and category ----
#Save p-values as object
phylo_disp_means_modules_2_cat_corr <- phylo_disparity_means_group_cat_modules_2[["PV.dist"]]
phylo_disp_means_modules_2_cat_pvals <- phylo_disparity_means_group_cat_modules_2[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
phylo_disparity_means_modules_2_cat_vars <- rownames(phylo_disp_means_modules_2_cat_corr)

#Replace string names to make them shorter
phylo_disparity_means_modules_2_cat_vars <- str_replace_all(phylo_disparity_means_modules_2_cat_vars, "\\.", "_")

#List 2 modules
modules_2_list <- c("braincase","rostrum")
modules_2_list_short <- c("B","R")

#Loop replacements categories
for (u in 1:length(categories_list)){
  phylo_disparity_means_modules_2_cat_vars <- str_replace_all(phylo_disparity_means_modules_2_cat_vars, categories_list[u], categories_list_short[u])
}

#Loop replacements groups
for (t in 1:length(groups_list)){
  phylo_disparity_means_modules_2_cat_vars <- str_replace_all(phylo_disparity_means_modules_2_cat_vars, groups_list[t], groups_list_short[t])
}

#Loop replacements modules
for (m in 1:length(modules_2_list)){
  phylo_disparity_means_modules_2_cat_vars <- str_replace_all(phylo_disparity_means_modules_2_cat_vars, modules_2_list[m], modules_2_list_short[m])
}

#Check it worked
phylo_disparity_means_modules_2_cat_vars

#Set correct row and col names for both
rownames(phylo_disp_means_modules_2_cat_corr) <- phylo_disparity_means_modules_2_cat_vars
rownames(phylo_disp_means_modules_2_cat_pvals) <- phylo_disparity_means_modules_2_cat_vars
colnames(phylo_disp_means_modules_2_cat_corr) <- phylo_disparity_means_modules_2_cat_vars
colnames(phylo_disp_means_modules_2_cat_pvals) <- phylo_disparity_means_modules_2_cat_vars

#Get upper triangles only - half matrix, eliminates redundant info
phylo_disp_means_modules_2_cat_corr_upper_tri <- get_upper_tri(phylo_disp_means_modules_2_cat_corr)
phylo_disp_means_modules_2_cat_pvals_upper_tri <- get_upper_tri(phylo_disp_means_modules_2_cat_pvals)

#Melt to make table in the format needed for heatmap
phylo_disp_means_modules_2_cat_corr_melt <- melt(phylo_disp_means_modules_2_cat_corr_upper_tri, value.name = "corr", na.rm = TRUE)
phylo_disp_means_modules_2_cat_pvals_melt <- melt(phylo_disp_means_modules_2_cat_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
phylo_disp_means_modules_2_cat_pvals_melt$corr <- phylo_disp_means_modules_2_cat_corr_melt$corr

#Create columns where only significant values are shown
phylo_disp_means_modules_2_cat_pvals_melt <- phylo_disp_means_modules_2_cat_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                      p_if_sig = ifelse(sig_p, p, NA),
                                                                                      corr_if_sig = ifelse(sig_p, corr, NA)) %>%
  mutate_at(vars(starts_with("corr")), list(~ round(., 3)))
phylo_disp_means_modules_2_cat_pvals_melt

#Nice heatmap plot
phylo_disparity_means_group_cat_modules_2_heatmap_ggplot <- ggplot(data = phylo_disp_means_modules_2_cat_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_mod_grp[9], high = mypalette_seq_mod_grp[1], mid = mypalette_seq_mod_grp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_mod_grp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity phylogenetically corrected means by module, group and growth stage")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_group_cat_modules_2_heatmap_ggplot

#Big plot, divide by group
#Make data frame for each group
phylo_disp_means_modules_2_cat_pvals_melt_mysticeti <- phylo_disp_means_modules_2_cat_pvals_melt %>% filter(str_detect(Var1, "Myst")) %>% 
  filter(str_detect(Var2, "Myst"))
phylo_disp_means_modules_2_cat_pvals_melt_odontoceti <- phylo_disp_means_modules_2_cat_pvals_melt %>% filter(str_detect(Var1, "Odont")) %>% 
  filter(str_detect(Var2, "Odont"))

#Create labels
category_modules_2_list <- as.factor(paste(rep(modules_2_list_short, times = length(categories_list_short)), 
                                           rep(categories_list_short, each = length(modules_2_list_short)), sep = "_"))

#Nice heatmap plot
phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_myst <- ggplot(data = phylo_disp_means_modules_2_cat_pvals_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
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
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_myst

phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_odont <- ggplot(data = phylo_disp_means_modules_2_cat_pvals_melt_odontoceti, aes(Var2, Var1, fill = p_if_sig))+
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
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_odont

plotD2 <- ggarrange(phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_myst, 
                    phylo_disparity_means_group_cat_modules_2_heatmap_ggplot_odont,
                    ncol = 2, nrow = 1, common.legend = F)
plotD2<-annotate_figure(plotD2, top = text_grob("Morphological disparity phylogenetically corrected means between module and growth stage by group", 
                                                face = "bold", size = 16, just = c(0.5,1)))
plotD2

###Line plots for Procrustes variance of groups at each stage ----
##All data
#Save variances as object
PV_group_cat_means_modules_2 <- data.frame(PV = phylo_disparity_means_group_cat_modules_2[["Procrustes.var"]])

#Look at order of variables to create vectors
PV_group_cat_means_modules_2

#Make vectors with correct order an number of categories and feeding
PV_modules_2 <- rep(modules_2_list, times = length(groups_list)*length(categories_list))

PV_categories_2 <- rep(categories_list, each = length(modules_2_list), times = length(groups_list))

PV_groups_2 <- rep(groups_list, each = length(PV_modules_2)/length(groups_list))

#Add labels and other attributes to tibble as columns
PV_group_cat_means_modules_2 <- PV_group_cat_means_modules_2 %>% 
  mutate(module = str_to_sentence(PV_modules_2), category= PV_categories_2, group = PV_groups_2)
PV_group_cat_means_modules_2

#Multiply PV for 100 as too small values crash ggplot
PV_group_cat_means_modules_2$PV_100 <- c(PV_group_cat_means_modules_2$PV * 100)
PV_group_cat_means_modules_2

#Nice line plot by module
PV_group_cat_means_modules_2_ggplot <- ggplot(PV_group_cat_means_modules_2, aes(x = category, y = PV_100)) + 
  geom_line(aes(x = category, y = PV_100,linetype = group, colour = module, group = group), linewidth = 1, inherit.aes = F)+
  geom_point(aes(color = module,
                 shape = group),size = 4, fill = "white", stroke = 1.5)+
  facet_wrap(vars(module))+
  scale_shape_manual(name = "Groups", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Modules", labels =levels(as.factor(PV_group_cat_means_modules_2$module)),
                      values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour","fill"))+   
  scale_linetype_manual(name = "Groups", labels = levels(groups),
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("PV (Procrustes variances) x 100")+
  xlab("Growth stage")+
  ggtitle ("Morphological disparity phylogenetically corrected means per module, group and growth stage")+ 
  scale_x_discrete(labels = levels(as.factor(categories_list_short)))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", linetype = "none", shape = "none")
PV_group_cat_means_modules_2_ggplot

#Add phylopics
PV_group_cat_means_modules_2_ggplot  <- 
  PV_group_cat_means_modules_2_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 1, y = 0.12, ysize = 0.005, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2, y = 0.06, ysize = 0.005, fill = "gray50")
PV_group_cat_means_modules_2_ggplot
#Delete extra silhouettes and change position

#Nice line plot by group
group_labels <- as_labeller(c('mysticeti' = "Mysticeti", 'odontoceti' = "Odontoceti"))

PV_module_cat_means_modules_2_ggplot <- ggplot(PV_group_cat_means_modules_2, aes(x = category, y = PV_100)) + 
  geom_line(aes(x = category, y = PV_100,linetype = module, colour = group, group = module), linewidth = 1, inherit.aes = F)+
  geom_point(aes(color = group,
                 shape = module),size = 4, fill = "white", stroke = 1.5)+
  facet_wrap(vars(group), labeller = group_labels)+
  scale_shape_manual(name = "Modules", labels = levels(as.factor(PV_group_cat_means_modules_2$module)), values = c(23,24))+
  scale_colour_manual(name = "Groups", labels =levels(groups),
                      values = mypalette_groups, aesthetics = c("colour","fill"))+   
  scale_linetype_manual(name = "Modules", labels = levels(as.factor(PV_group_cat_means_modules_2$module)),
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("PV (Procrustes variances) x 100")+
  xlab("Growth stage")+
  ggtitle ("Morphological disparity phylogenetically corrected means per group, module and growth stage")+ 
  scale_x_discrete(labels = levels(as.factor(categories_list_short)))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", shape = guide_legend(keywidth=unit(3, "char"),override.aes = list(size = 3)))
PV_module_cat_means_modules_2_ggplot

#Add silhouettes groups
PV_module_cat_means_modules_2_ggplot  <- 
  PV_module_cat_means_modules_2_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 2, y = 0.06, ysize = 0.005, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 3, y = 0.1, ysize = 0.005, fill = mypalette_groups[2])
PV_module_cat_means_modules_2_ggplot
#Delete extra silhouettes and change position

###### 

#Next - ch. 7b - Trajectory analyses - phylogenetically transformed components means