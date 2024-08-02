#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH. 8 - Allometry analyses rostrum vs braincase, mcp testing  - mean shapes

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
library(emmeans)

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
allometry_phylo_means_rostrum_list <- list()
  
allometry_phylo_means_rostrum_list[[1]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list[[2]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list[[3]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list[[4]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_rostrum <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_rostrum[[c]] <- summary(allometry_phylo_means_rostrum_list[[c]])
  names(summary_allometry_phylo_means_rostrum)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_rostrum

#Allometry by category by group models
#Combination
allometry_phylo_means_rostrum_list_comb <- list()

allometry_phylo_means_rostrum_list_comb[[1]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[1]]] + gdf_mean_rostrum$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_comb[[2]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[2]]] + gdf_mean_rostrum$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_comb[[3]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[3]]] + gdf_mean_rostrum$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_comb[[4]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[4]]] + gdf_mean_rostrum$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_rostrum_comb <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_rostrum_comb[[c]] <- summary(allometry_phylo_means_rostrum_list_comb[[c]])
  names(summary_allometry_phylo_means_rostrum_comb)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_rostrum_comb

#Interaction
allometry_phylo_means_rostrum_list_int <- list()

allometry_phylo_means_rostrum_list_int[[1]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[1]]] * gdf_mean_rostrum$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_int[[2]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[2]]] * gdf_mean_rostrum$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_int[[3]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[3]]] * gdf_mean_rostrum$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_rostrum_list_int[[4]] <- procD.pgls(gdf_mean_rostrum$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_rostrum$size[rows_categories_mean_all[[4]]] * gdf_mean_rostrum$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_rostrum_int <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_rostrum_int[[c]] <- summary(allometry_phylo_means_rostrum_list_int[[c]])
  names(summary_allometry_phylo_means_rostrum_int)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_rostrum_int

#Save results of significant regression to file
sink("Output/8-Allometry means/allometry_phylo_means_shape_size_grp_cat_rostrum.txt")
print("Null")
summary_allometry_phylo_means_rostrum

print("Combination +")
summary_allometry_phylo_means_rostrum_comb

print("Interaction *")
summary_allometry_phylo_means_rostrum_int
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_phylo_means_models_grp_cat_rostrum_list <- list()
  
#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_rostrum_list[[c]] <- anova(allometry_phylo_means_rostrum_list[[c]], allometry_phylo_means_rostrum_list_comb[[c]], allometry_phylo_means_rostrum_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_list_rostrum_list)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_rostrum_list

#ANOVAs - is a model significantly better than the others? - check between comb and int only
anova_allometry_phylo_means_models_grp_cat_rostrum_list_1 <- list()

#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_rostrum_list_1[[c]] <- anova(allometry_phylo_means_rostrum_list_comb[[c]], allometry_phylo_means_rostrum_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_rostrum_list_1)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_rostrum_list_1

#In all categories, comb better than null but interaction not difference - no sign diff in slope between groups

##Extract RegScores best model for comparison
#Regression score of shape vs logCS - regression method with "RegScore" plotting

allometry_phylo_means_rostrum_plot_regscore <- list()

for (c in 1:length(categories_list)){  
allometry_phylo_means_rostrum_plot_regscore[[c]] <- plot(allometry_phylo_means_rostrum_list_comb[[c]], type = "regression",predictor = gdf_mean_rostrum$size[rows_categories_mean_all[[c]]], reg.type = "RegScore", 
                                      main = paste("Shape vs logCS by group", "-", levels(categories)[[c]]), xlab = "logCS", pch = 21, cex = 1.2, font.main = 2,
                                      col = mypalette_category[c], bg = mypalette_category[c])
text(x = gdf_mean_rostrum$size[rows_categories_mean_all[[c]]], y = allometry_phylo_means_rostrum_plot_regscore[[c]]$RegScore, labels = gdf_mean_rostrum$genus[rows_categories_mean_all[[c]]],
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels  
}

#Create object to use for linear model
allometry_phylo_means_rostrum_regscores <- rbind(allometry_phylo_means_rostrum_plot_regscore[[1]][["RegScore"]], allometry_phylo_means_rostrum_plot_regscore[[2]][["RegScore"]] ,
                                                 allometry_phylo_means_rostrum_plot_regscore[[3]][["RegScore"]] ,allometry_phylo_means_rostrum_plot_regscore[[4]][["RegScore"]] )

#Linear model for line by group and category 
allometry_phylo_means_rostrum_regscores_df <- data.frame(RegScores = allometry_phylo_means_rostrum_regscores, logCS = gdf_mean_rostrum$size, genus = gdf_mean_rostrum$genus, group = gdf_mean_rostrum$group,
                                                       category = gdf_mean_rostrum$category)

allometry_phylo_means_rostrum_regscores_df$grp_cat <- interaction(allometry_phylo_means_rostrum_regscores_df$group, allometry_phylo_means_rostrum_regscores_df$category)



###Pairwise comparison of regression model between groups and category ----
#Create models, with different slopes and int or just int
allometry_phylo_means_rostrum_null_cat <- lm.rrpp(RegScores ~ logCS * category,   #null model only differences by category
                                     data = allometry_phylo_means_rostrum_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_rostrum_null_grp <- lm.rrpp(RegScores ~ logCS * group,   #null model only differences by group
                                                  data = allometry_phylo_means_rostrum_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_rostrum_comb <- lm.rrpp(RegScores ~ logCS * category + group,
                                     data = allometry_phylo_means_rostrum_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_rostrum_int <- lm.rrpp(RegScores ~ logCS * category * group,
                                    data = allometry_phylo_means_rostrum_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_rostrum_null_cat)
summary(allometry_phylo_means_rostrum_null_grp)
summary(allometry_phylo_means_rostrum_comb)
summary(allometry_phylo_means_rostrum_int)

#Anova for difference between models
anova_allometry_phylo_means_models_grp_cat_rostrum <- anova(allometry_phylo_means_rostrum_null_cat,allometry_phylo_means_rostrum_null_grp, allometry_phylo_means_rostrum_comb, allometry_phylo_means_rostrum_int)
anova_allometry_phylo_means_models_grp_cat_rostrum

#Check interactions in models
anova(allometry_phylo_means_rostrum_null_cat)
anova(allometry_phylo_means_rostrum_null_grp)
anova(allometry_phylo_means_rostrum_comb)
anova(allometry_phylo_means_rostrum_int)

#No sign interaction between category and group. Group and category don't influence regression

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_rostrum_int1 <- lm(RegScores ~ logCS * grp_cat,
                                data = allometry_phylo_means_rostrum_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_rostrum_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_rostrum_emms <- emmeans(allometry_phylo_means_rostrum_int1, "grp_cat")

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_rostrum_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_rostrum_ems_table <- pwpm(allometry_phylo_means_rostrum_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_rostrum_ems_table

#Nothing significant

#Save results to file
sink("Output/8-Allometry means/pairwise_allometry_phylo_means_rostrum.txt")
print("ANOVA models")
print(anova_allometry_phylo_means_models_grp_cat_rostrum)

print("summary models - lmrpp")
anova(allometry_phylo_means_rostrum_null_cat)
anova(allometry_phylo_means_rostrum_null_grp)
anova(allometry_phylo_means_rostrum_comb)
anova(allometry_phylo_means_rostrum_int)

print("summary model used for comparisons - lm")
summary(allometry_phylo_means_rostrum_int1)
anova(allometry_phylo_means_rostrum_int1)

print("Pairwise comparison using emmeans")
summary(allometry_phylo_means_rostrum_emms)

print("Full results table emmeans pairwise comparions")
pwpm(allometry_phylo_means_rostrum_emms)
sink()

####Plot allometry rostrum by category by group ----
#Create data frame object that ggplot can read - use data from plot object you want to improve
#Convert data frame to tibble
allometry_phylo_means_rostrum_grp_cat_plot <- as_tibble(allometry_phylo_means_rostrum_regscores_df)

#Plot allometry regression by category by group 
allometry_phylo_means_rostrum_grp_cat_ggplot <- ggplot(allometry_phylo_means_rostrum_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(aes(colour = category,fill = category), size = 0, alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group,colour = category,fill = category, group = grp_cat), inherit.aes = F,        
              se = F, linewidth = 1.5, alpha = 1)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle("Rostrum")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0), legend.box = "vertical")+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_means_rostrum_grp_cat_ggplot

#Add phylopic
allometry_phylo_means_rostrum_grp_cat_ggplot <- 
  allometry_phylo_means_rostrum_grp_cat_ggplot +
  add_phylopic(myst, alpha = 1, x = 3.1, y = 0.035, ysize = 0.007, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2.6, y = 0.005, ysize = 0.005, fill = "gray50")
allometry_phylo_means_rostrum_grp_cat_ggplot

###Braincase ----

##Regression shape on logCS size braincase
allometry_phylo_means_braincase_list <- list()

allometry_phylo_means_braincase_list[[1]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list[[2]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list[[3]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list[[4]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_braincase <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_braincase[[c]] <- summary(allometry_phylo_means_braincase_list[[c]])
  names(summary_allometry_phylo_means_braincase)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_braincase

#Allometry by category by group models
#Combination
allometry_phylo_means_braincase_list_comb <- list()

allometry_phylo_means_braincase_list_comb[[1]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[1]]] + gdf_mean_braincase$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_comb[[2]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[2]]] + gdf_mean_braincase$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_comb[[3]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[3]]] + gdf_mean_braincase$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_comb[[4]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[4]]] + gdf_mean_braincase$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_braincase_comb <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_braincase_comb[[c]] <- summary(allometry_phylo_means_braincase_list_comb[[c]])
  names(summary_allometry_phylo_means_braincase_comb)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_braincase_comb

#Interaction
allometry_phylo_means_braincase_list_int <- list()

allometry_phylo_means_braincase_list_int[[1]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[1]]] * gdf_mean_braincase$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_int[[2]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[2]]] * gdf_mean_braincase$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_int[[3]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[3]]] * gdf_mean_braincase$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_braincase_list_int[[4]] <- procD.pgls(gdf_mean_braincase$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_braincase$size[rows_categories_mean_all[[4]]] * gdf_mean_braincase$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_braincase_int <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_braincase_int[[c]] <- summary(allometry_phylo_means_braincase_list_int[[c]])
  names(summary_allometry_phylo_means_braincase_int)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_braincase_int

#Save results of significant regression to file
sink("Output/8-Allometry means/allometry_phylo_means_shape_size_grp_cat_braincase.txt")
print("Null")
summary_allometry_phylo_means_braincase

print("Combination +")
summary_allometry_phylo_means_braincase_comb

print("Interaction *")
summary_allometry_phylo_means_braincase_int
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_phylo_means_models_grp_cat_braincase_list <- list()

#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_braincase_list[[c]] <- anova(allometry_phylo_means_braincase_list[[c]], allometry_phylo_means_braincase_list_comb[[c]], allometry_phylo_means_braincase_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_braincase_list)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_braincase_list

#ANOVAs - is a model significantly better than the others? - check between comb and int only
anova_allometry_phylo_means_models_grp_cat_braincase_list_1 <- list()

#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_braincase_list_1[[c]] <- anova(allometry_phylo_means_braincase_list_comb[[c]], allometry_phylo_means_braincase_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_braincase_list_1)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_braincase_list_1

#In all categories, comb better than null but interaction not difference - no sign diff in slope between groups

##Extract RegScores best model for comparison
#Regression score of shape vs logCS - regression method with "RegScore" plotting

allometry_phylo_means_braincase_plot_regscore <- list()

for (c in 1:length(categories_list)){  
  allometry_phylo_means_braincase_plot_regscore[[c]] <- plot(allometry_phylo_means_braincase_list_comb[[c]], type = "regression",predictor = gdf_mean_braincase$size[rows_categories_mean_all[[c]]], reg.type = "RegScore", 
                                                           main = paste("Shape vs logCS by group", "-", levels(categories)[[c]]), xlab = "logCS", pch = 21, cex = 1.2, font.main = 2,
                                                           col = mypalette_category[c], bg = mypalette_category[c])
  text(x = gdf_mean_braincase$size[rows_categories_mean_all[[c]]], y = allometry_phylo_means_braincase_plot_regscore[[c]]$RegScore, labels = gdf_mean_braincase$genus[rows_categories_mean_all[[c]]],
       pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels  
}

#Create object to use for linear model
allometry_phylo_means_braincase_regscores <- rbind(allometry_phylo_means_braincase_plot_regscore[[1]][["RegScore"]], allometry_phylo_means_braincase_plot_regscore[[2]][["RegScore"]] ,
                                                 allometry_phylo_means_braincase_plot_regscore[[3]][["RegScore"]] ,allometry_phylo_means_braincase_plot_regscore[[4]][["RegScore"]] )

#Linear model for line by group and category 
allometry_phylo_means_braincase_regscores_df <- data.frame(RegScores = allometry_phylo_means_braincase_regscores, logCS = gdf_mean_braincase$size, genus = gdf_mean_braincase$genus, group = gdf_mean_braincase$group,
                                                         category = gdf_mean_braincase$category)

allometry_phylo_means_braincase_regscores_df$grp_cat <- interaction(allometry_phylo_means_braincase_regscores_df$group, allometry_phylo_means_braincase_regscores_df$category)



###Pairwise comparison of regression model between groups and category ----
#Create models, with different slopes and int or just int
allometry_phylo_means_braincase_null_cat <- lm.rrpp(RegScores ~ logCS * category,   #null model only differences by category
                                                  data = allometry_phylo_means_braincase_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_braincase_null_grp <- lm.rrpp(RegScores ~ logCS * group,   #null model only differences by group
                                                  data = allometry_phylo_means_braincase_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_braincase_comb <- lm.rrpp(RegScores ~ logCS * category + group,
                                              data = allometry_phylo_means_braincase_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_braincase_int <- lm.rrpp(RegScores ~ logCS * category * group,
                                             data = allometry_phylo_means_braincase_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_braincase_null_cat)
summary(allometry_phylo_means_braincase_null_grp)
summary(allometry_phylo_means_braincase_comb)
summary(allometry_phylo_means_braincase_int)

#Anova for difference between models
anova_allometry_phylo_means_models_grp_cat_braincase <- anova(allometry_phylo_means_braincase_null_cat,allometry_phylo_means_braincase_null_grp, allometry_phylo_means_braincase_comb, allometry_phylo_means_braincase_int)
anova_allometry_phylo_means_models_grp_cat_braincase

#Check interactions in models
anova(allometry_phylo_means_braincase_null_cat)
anova(allometry_phylo_means_braincase_null_grp)
anova(allometry_phylo_means_braincase_comb)
anova(allometry_phylo_means_braincase_int)

#Sign interaction between size and category (diff slopes). Groups different in intercepts

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_braincase_int1 <- lm(RegScores ~ logCS * grp_cat,
                                         data = allometry_phylo_means_braincase_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_braincase_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_braincase_emms <- emmeans(allometry_phylo_means_braincase_int1, "grp_cat")

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_braincase_ems_table <- pwpm(allometry_phylo_means_braincase_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_braincase_ems_table

#Nothing significant

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_braincase_comb1 <- lm(RegScores ~ logCS + grp_cat,
                                           data = allometry_phylo_means_braincase_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_braincase_comb1)

#Get pairwise comparisons of slopes
allometry_phylo_means_braincase_comb_emms <- emmeans(allometry_phylo_means_braincase_comb1, "grp_cat")

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_braincase_comb_ems_table <- pwpm(allometry_phylo_means_braincase_comb_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_braincase_comb_ems_table

#Only sign differences in Mysticeti early vs Odontoceti early and adult

#Save results to file
sink("Output/8-Allometry means/pairwise_allometry_phylo_means_braincase.txt")
print("ANOVA models")
print(anova_allometry_phylo_means_models_grp_cat_braincase)

print("summary models - lmrpp")
anova(allometry_phylo_means_braincase_null_cat)
anova(allometry_phylo_means_braincase_null_grp)
anova(allometry_phylo_means_braincase_comb)
anova(allometry_phylo_means_braincase_int)

print("summary model used for comparisons - lm")
summary(allometry_phylo_means_braincase_int1)
anova(allometry_phylo_means_braincase_int1)

print("Pairwise comparison using emmeans")
summary(allometry_phylo_means_braincase_emms)

print("Full results table emmeans pairwise comparions")
pwpm(allometry_phylo_means_braincase_emms)

print("summary model used for comparisons comb - lm")
summary(allometry_phylo_means_braincase_comb1)
anova(allometry_phylo_means_braincase_comb1)

print("Pairwise comparison using emmeans")
summary(allometry_phylo_means_braincase_comb_emms)

print("Full results table emmeans pairwise comparions")
pwpm(allometry_phylo_means_braincase_comb_emms)
sink()

####Plot allometry braincase by category by group ----
#Create data frame object that ggplot can read - use data from plot object you want to improve
#Convert data frame to tibble
allometry_phylo_means_braincase_grp_cat_plot <- as_tibble(allometry_phylo_means_braincase_regscores_df)

#Plot allometry regression by category by group 
allometry_phylo_means_braincase_grp_cat_ggplot <- ggplot(allometry_phylo_means_braincase_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(aes(colour = category,fill = category), size = 0, alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group,colour = category,fill = category, group = grp_cat), inherit.aes = F,        
              se = F, linewidth = 1.5, alpha = 1)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle("Braincase")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0), legend.box = "vertical")+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_means_braincase_grp_cat_ggplot

#Add phylopic
allometry_phylo_means_braincase_grp_cat_ggplot <- 
  allometry_phylo_means_braincase_grp_cat_ggplot +
  add_phylopic(myst, alpha = 1, x = 3.05, y = -0.07, ysize = 0.009, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2.8, y = 0.02, ysize = 0.007, fill = "gray50")
allometry_phylo_means_braincase_grp_cat_ggplot

ggarrange(allometry_phylo_means_rostrum_grp_cat_ggplot, allometry_phylo_means_braincase_grp_cat_ggplot, common.legend = T, legend = "bottom")

#ALLOMETRY DIFFERENCE ROSTRUM AND BRAINCASE COMPARE----
##Evaluate allometry using logCS for entire skull as size
##Use pcscores instead of coords - all same length of variables and same results as coords

##By group, category and module ----
#Make common df regscores both modules
#Add module column first
allometry_phylo_means_rostrum_regscores_df$module <- "rostrum"
allometry_phylo_means_braincase_regscores_df$module <- "braincase"

allometry_phylo_means_regscores_df <- rbind(allometry_phylo_means_rostrum_regscores_df, allometry_phylo_means_braincase_regscores_df)

#Create models, with different slopes and int or just int
allometry_phylo_means_mod_null <- lm.rrpp(RegScores ~ logCS * module,   #null model only differences by category
                                                    data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_comb <- lm.rrpp(RegScores ~ logCS * module + grp_cat,
                                                data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_int <- lm.rrpp(RegScores ~ logCS * module * grp_cat,
                                               data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_mod_null)
summary(allometry_phylo_means_mod_comb)
summary(allometry_phylo_means_mod_int)

#Anova for difference between models
anova_allometry_phylo_means_models_mod_grp_cat <- anova(allometry_phylo_means_mod_null, allometry_phylo_means_mod_comb, allometry_phylo_means_mod_int)
anova_allometry_phylo_means_models_mod_grp_cat

#Check interactions in models
anova(allometry_phylo_means_mod_null)
anova(allometry_phylo_means_mod_comb)
anova(allometry_phylo_means_mod_int)

#Sign difference in slopes between modules, sign diff in intercept between group and category in the 2 modules

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_mod_int1 <- lm(RegScores ~ logCS * module * grp_cat,
                                           data = allometry_phylo_means_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_mod_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_mod_emms <- emmeans(allometry_phylo_means_mod_int1, ~ module * grp_cat)

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_mod_ems_table <- pwpm(allometry_phylo_means_mod_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_mod_ems_table

#Nothing sign

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_mod_null1 <- lm(RegScores ~ logCS * module,
                                     data = allometry_phylo_means_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_mod_null1)

#Get pairwise comparisons of slopes
allometry_phylo_means_mod_null_emms <- emmeans(allometry_phylo_means_mod_null1, ~ module)

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_mod_null_ems_table <- pwpm(allometry_phylo_means_mod_null_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_mod_null_ems_table

#Nothing sign

##By category and module ----
#Make common df regscores both modules

#Create models, with different slopes and int or just int
allometry_phylo_means_mod_cat_null <- lm.rrpp(RegScores ~ logCS * module,   #null model only differences by category
                                          data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_cat_comb <- lm.rrpp(RegScores ~ logCS * module + category,
                                          data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_cat_int <- lm.rrpp(RegScores ~ logCS * module * category,
                                         data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_mod_cat_null)
summary(allometry_phylo_means_mod_cat_comb)
summary(allometry_phylo_means_mod_cat_int)

#Anova for difference between models
anova_allometry_phylo_means_models_mod_cat_grp_cat <- anova(allometry_phylo_means_mod_cat_null, allometry_phylo_means_mod_cat_comb, allometry_phylo_means_mod_cat_int)
anova_allometry_phylo_means_models_mod_cat_grp_cat

#Check interactions in models
anova(allometry_phylo_means_mod_cat_null)
anova(allometry_phylo_means_mod_cat_comb)
anova(allometry_phylo_means_mod_cat_int)

#Sign difference in slopes between modules, sign diff in intercept between group and category in the 2 modules

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_mod_cat_int1 <- lm(RegScores ~ logCS * module * category,
                                     data = allometry_phylo_means_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_mod_cat_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_mod_cat_emms <- emmeans(allometry_phylo_means_mod_cat_int1, ~ module * category)

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_mod_cat_ems_table <- pwpm(allometry_phylo_means_mod_cat_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_mod_cat_ems_table

#Nothing sign

##By group and module ----
#Make common df regscores both modules

#Create models, with different slopes and int or just int
allometry_phylo_means_mod_grp_null <- lm.rrpp(RegScores ~ logCS * module,   #null model only differences by category
                                              data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_grp_comb <- lm.rrpp(RegScores ~ logCS * module + group,
                                              data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_mod_grp_int <- lm.rrpp(RegScores ~ logCS * module * group,
                                             data = allometry_phylo_means_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_mod_grp_null)
summary(allometry_phylo_means_mod_grp_comb)
summary(allometry_phylo_means_mod_grp_int)

#Anova for difference between models
anova_allometry_phylo_means_models_mod_grp_grp_grp <- anova(allometry_phylo_means_mod_grp_null, allometry_phylo_means_mod_grp_comb, allometry_phylo_means_mod_grp_int)
anova_allometry_phylo_means_models_mod_grp_grp_grp

#Check interactions in models
anova(allometry_phylo_means_mod_grp_null)
anova(allometry_phylo_means_mod_grp_comb)
anova(allometry_phylo_means_mod_grp_int)

#Sign difference in slopes between modules, sign diff in intercept between group and category in the 2 modules

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_mod_grp_int1 <- lm(RegScores ~ logCS * module * group,
                                         data = allometry_phylo_means_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_mod_grp_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_mod_grp_emms <- emmeans(allometry_phylo_means_mod_grp_int1, ~ module * group)

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_braincase_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_mod_grp_ems_table <- pwpm(allometry_phylo_means_mod_grp_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_mod_grp_ems_table

#Nothing sign

###Plot allometry rostrum, braincase divided by group and category ----

#Create common tibble for plotting
allometry_phylo_means_grp_cat_plot <- as_tibble(allometry_phylo_means_regscores_df)
allometry_phylo_means_grp_cat_plot$module <- str_to_title(allometry_phylo_means_grp_cat_plot$module)
glimpse(allometry_phylo_means_grp_cat_plot)

####By group ----
#Divide by groups - too many lines all together
allometry_phylo_means_grp_cat_plot_mysticeti <- allometry_phylo_means_grp_cat_plot %>% filter(group == "mysticeti")
allometry_phylo_means_grp_cat_plot_odontoceti <- allometry_phylo_means_grp_cat_plot %>% filter(group == "odontoceti")

#Plot allometry regression by category mysticeti
allometry_phylo_means_group_cat_module_ggplot_mysticeti <- ggplot(allometry_phylo_means_grp_cat_plot_mysticeti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = category, color = module), inherit.aes = F,       
              linewidth = 1.2, alpha =1, se = F, show.legend = T)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", fill = "none", linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_means_group_cat_module_ggplot_mysticeti

#Add phylopic
allometry_phylo_means_group_cat_module_ggplot_mysticeti <- 
  allometry_phylo_means_group_cat_module_ggplot_mysticeti +
  add_phylopic(myst, alpha = 1, x = 2.75, y = -0.02, ysize = 0.003, fill = "gray30")
allometry_phylo_means_group_cat_module_ggplot_mysticeti
#Fix silhouette size

#Plot allometry regression by category odontoceti
allometry_phylo_means_group_cat_module_ggplot_odontoceti <- ggplot(allometry_phylo_means_grp_cat_plot_odontoceti, aes(x = logCS, y = RegScores))+
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
  guides(colour = "none", fill = "none", linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_means_group_cat_module_ggplot_odontoceti

#Add phylopic
allometry_phylo_means_group_cat_module_ggplot_odontoceti <- 
  allometry_phylo_means_group_cat_module_ggplot_odontoceti +
  add_phylopic(odont, alpha = 1, x = 2.6, y = -0.03, ysize = 0.004, fill = "gray50")
allometry_phylo_means_group_cat_module_ggplot_odontoceti
#Fix silhouette size

ggarrange(allometry_phylo_means_group_cat_module_ggplot_mysticeti,allometry_phylo_means_group_cat_module_ggplot_odontoceti, 
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

###Line plots ----
#geom_abline plot by group faceted by module

#Make separate dfs for each module - easier to extract coefficients correctly
allometry_phylo_means_grp_cat_plot_braincase <- allometry_phylo_means_grp_cat_plot %>% filter(module == "Braincase")
allometry_phylo_means_grp_cat_plot_rostrum <- allometry_phylo_means_grp_cat_plot %>% filter(module == "Rostrum")

####By group and module ----
#Braincase
#Linear model for line by group and category
allometry_phylo_means_group_regline_braincase <- lm(RegScores ~ logCS * group, data = allometry_phylo_means_grp_cat_plot_braincase)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_means_group_regline_coeffs_braincase <- as.matrix(allometry_phylo_means_group_regline_braincase$coefficients)

#Save intercepts and slopes separately
allometry_phylo_means_group_regline_intercepts_braincase <- as.matrix(allometry_phylo_means_group_regline_coeffs_braincase[c(1, 3:(length(levels(groups))+1)),])
allometry_phylo_means_group_regline_slopes_braincase <- as.matrix(allometry_phylo_means_group_regline_coeffs_braincase[c(2, length(levels(groups))+2:(length(levels(groups)))),])

#Calculate real intercepts and slopes
allometry_phylo_means_group_regline_intercepts_ok_braincase <- as.matrix(c(allometry_phylo_means_group_regline_intercepts_braincase[1,], allometry_phylo_means_group_regline_intercepts_braincase[1,]+
                                                                       allometry_phylo_means_group_regline_intercepts_braincase[2:length(allometry_phylo_means_group_regline_intercepts_braincase),]))

allometry_phylo_means_group_regline_slopes_ok_braincase <- as.matrix(c(allometry_phylo_means_group_regline_slopes_braincase[1,], allometry_phylo_means_group_regline_slopes_braincase[1,]+
                                                                   allometry_phylo_means_group_regline_slopes_braincase[2:length(allometry_phylo_means_group_regline_slopes_braincase),]))

#Save as data frame with grouping variables
allometry_phylo_means_group_coeffs_braincase <- data.frame(Slope = allometry_phylo_means_group_regline_slopes_ok_braincase, Intercept = allometry_phylo_means_group_regline_intercepts_ok_braincase, 
                                                     row.names = levels(groups))
#Check for NA and other issues
allometry_phylo_means_group_coeffs_braincase 

#Add classifiers
allometry_phylo_means_group_coeffs_braincase <- allometry_phylo_means_group_coeffs_braincase %>% mutate(group = levels(groups),
                                                                                            module = 'Braincase')
allometry_phylo_means_group_coeffs_braincase

#Rostrum
#Linear model for line by group and category
allometry_phylo_means_group_regline_rostrum <- lm(RegScores ~ logCS * group, data = allometry_phylo_means_grp_cat_plot_rostrum)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_means_group_regline_coeffs_rostrum <- as.matrix(allometry_phylo_means_group_regline_rostrum$coefficients)

#Save intercepts and slopes separately
allometry_phylo_means_group_regline_intercepts_rostrum <- as.matrix(allometry_phylo_means_group_regline_coeffs_rostrum[c(1, 3:(length(levels(groups))+1)),])
allometry_phylo_means_group_regline_slopes_rostrum <- as.matrix(allometry_phylo_means_group_regline_coeffs_rostrum[c(2, length(levels(groups))+2:(length(levels(groups)))),])

#Calculate real intercepts and slopes
allometry_phylo_means_group_regline_intercepts_ok_rostrum <- as.matrix(c(allometry_phylo_means_group_regline_intercepts_rostrum[1,], allometry_phylo_means_group_regline_intercepts_rostrum[1,]+
                                                                     allometry_phylo_means_group_regline_intercepts_rostrum[2:length(allometry_phylo_means_group_regline_intercepts_rostrum),]))

allometry_phylo_means_group_regline_slopes_ok_rostrum <- as.matrix(c(allometry_phylo_means_group_regline_slopes_rostrum[1,], allometry_phylo_means_group_regline_slopes_rostrum[1,]+
                                                                 allometry_phylo_means_group_regline_slopes_rostrum[2:length(allometry_phylo_means_group_regline_slopes_rostrum),]))

#Save as data frame with grouping variables
allometry_phylo_means_group_coeffs_rostrum <- data.frame(Slope = allometry_phylo_means_group_regline_slopes_ok_rostrum, Intercept = allometry_phylo_means_group_regline_intercepts_ok_rostrum, 
                                                   row.names = levels(groups))
#Check for NA and other issues
allometry_phylo_means_group_coeffs_rostrum 

#Add classifiers
allometry_phylo_means_group_coeffs_rostrum <- allometry_phylo_means_group_coeffs_rostrum %>% mutate(group = levels(groups),
                                                                                        module = 'Rostrum')
allometry_phylo_means_group_coeffs_rostrum

#Put together
allometry_phylo_means_group_coeffs <- rbind(allometry_phylo_means_group_coeffs_braincase, allometry_phylo_means_group_coeffs_rostrum)
allometry_phylo_means_group_coeffs

#Make sure group correspond to use facet_wrap
allometry_phylo_means_grp_cat_plot$group <- str_to_title(allometry_phylo_means_grp_cat_plot$group)

##Plot by module both groups
allometry_phylo_means_grp_line_ggplot <- ggplot(allometry_phylo_means_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_means_group_coeffs, 
              aes(intercept = Intercept, slope = Slope,  colour = group, linetype = module), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Groups", labels =  levels(groups), 
                      values = mypalette_groups, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Modules",
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  theme(legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0))+
  guides(colour = "none",
         linetype = guide_legend(keywidth = unit(3, "char"), override.aes = list(colour = c("gray20"))))
allometry_phylo_means_grp_line_ggplot

#Add silhouettes groups
allometry_phylo_means_grp_line_ggplot  <- 
  allometry_phylo_means_grp_line_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 2.45, y = 0.002, ysize = 0.007, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 3.2, y = 0, ysize = 0.006, fill = mypalette_groups[2])
allometry_phylo_means_grp_line_ggplot

####By group and category by module ----
#Braincase
#Linear model for line by group and category
allometry_phylo_means_group_cat_regline_braincase <- lm(RegScores ~ logCS * grp_cat, data = allometry_phylo_means_grp_cat_plot_braincase)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_means_group_cat_regline_coeffs_braincase <- as.matrix(allometry_phylo_means_group_cat_regline_braincase$coefficients)

#Save intercepts and slopes separately
group_cat_vars <- levels(allometry_phylo_means_grp_cat_plot$grp_cat)

allometry_phylo_means_group_cat_regline_intercepts_braincase <- as.matrix(allometry_phylo_means_group_cat_regline_coeffs_braincase[c(1, 3:(length(group_cat_vars)+1)),])
allometry_phylo_means_group_cat_regline_slopes_braincase <- as.matrix(allometry_phylo_means_group_cat_regline_coeffs_braincase[c(2, length(group_cat_vars)+2:(length(group_cat_vars))),])

#Calculate real intercepts and slopes
allometry_phylo_means_group_cat_regline_intercepts_ok_braincase <- as.matrix(c(allometry_phylo_means_group_cat_regline_intercepts_braincase[1,], allometry_phylo_means_group_cat_regline_intercepts_braincase[1,]+
                                                                           allometry_phylo_means_group_cat_regline_intercepts_braincase[2:length(allometry_phylo_means_group_cat_regline_intercepts_braincase),]))

allometry_phylo_means_group_cat_regline_slopes_ok_braincase <- as.matrix(c(allometry_phylo_means_group_cat_regline_slopes_braincase[1,], allometry_phylo_means_group_cat_regline_slopes_braincase[1,]+
                                                                       allometry_phylo_means_group_cat_regline_slopes_braincase[2:length(allometry_phylo_means_group_cat_regline_slopes_braincase),]))

#Save as data frame with grouping variables
allometry_phylo_means_group_cat_coeffs_braincase <- data.frame(Slope = allometry_phylo_means_group_cat_regline_slopes_ok_braincase, Intercept = allometry_phylo_means_group_cat_regline_intercepts_ok_braincase, 
                                                         row.names = group_cat_vars )
#Check for NA and other issues
allometry_phylo_means_group_cat_coeffs_braincase 

#Add classifiers
allometry_phylo_means_group_cat_coeffs_braincase <- allometry_phylo_means_group_cat_coeffs_braincase %>% mutate(category = rep(categories_list, each = 2), 
                                                                                                    group = rep(str_to_title(groups_list), times = 4),
                                                                                                    module = 'Braincase')
allometry_phylo_means_group_cat_coeffs_braincase

#Rostrum
#Linear model for line by group and category
allometry_phylo_means_group_cat_regline_rostrum <- lm(RegScores ~ logCS * grp_cat, data = allometry_phylo_means_grp_cat_plot_rostrum)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_means_group_cat_regline_coeffs_rostrum <- as.matrix(allometry_phylo_means_group_cat_regline_rostrum$coefficients)

#Save intercepts and slopes separately
group_cat_vars <- levels(allometry_phylo_means_grp_cat_plot$grp_cat)

allometry_phylo_means_group_cat_regline_intercepts_rostrum <- as.matrix(allometry_phylo_means_group_cat_regline_coeffs_rostrum[c(1, 3:(length(group_cat_vars)+1)),])
allometry_phylo_means_group_cat_regline_slopes_rostrum <- as.matrix(allometry_phylo_means_group_cat_regline_coeffs_rostrum[c(2, length(group_cat_vars)+2:(length(group_cat_vars))),])

#Calculate real intercepts and slopes
allometry_phylo_means_group_cat_regline_intercepts_ok_rostrum <- as.matrix(c(allometry_phylo_means_group_cat_regline_intercepts_rostrum[1,], allometry_phylo_means_group_cat_regline_intercepts_rostrum[1,]+
                                                                         allometry_phylo_means_group_cat_regline_intercepts_rostrum[2:length(allometry_phylo_means_group_cat_regline_intercepts_rostrum),]))

allometry_phylo_means_group_cat_regline_slopes_ok_rostrum <- as.matrix(c(allometry_phylo_means_group_cat_regline_slopes_rostrum[1,], allometry_phylo_means_group_cat_regline_slopes_rostrum[1,]+
                                                                     allometry_phylo_means_group_cat_regline_slopes_rostrum[2:length(allometry_phylo_means_group_cat_regline_slopes_rostrum),]))

#Save as data frame with grouping variables
allometry_phylo_means_group_cat_coeffs_rostrum <- data.frame(Slope = allometry_phylo_means_group_cat_regline_slopes_ok_rostrum, Intercept = allometry_phylo_means_group_cat_regline_intercepts_ok_rostrum, 
                                                       row.names = group_cat_vars )
#Check for NA and other issues
allometry_phylo_means_group_cat_coeffs_rostrum 

#Add classifiers
allometry_phylo_means_group_cat_coeffs_rostrum <- allometry_phylo_means_group_cat_coeffs_rostrum %>% mutate(category = rep(categories_list, each = 2), 
                                                                                                group = rep(str_to_title(groups_list), times = 4),
                                                                                                module = 'Rostrum')
allometry_phylo_means_group_cat_coeffs_rostrum

#Put together
allometry_phylo_means_group_cat_coeffs <- rbind(allometry_phylo_means_group_cat_coeffs_braincase, allometry_phylo_means_group_cat_coeffs_rostrum)
allometry_phylo_means_group_cat_coeffs

#Make sure group correspond to use facet_wrap
allometry_phylo_means_braincase_grp_cat_plot$group <- str_to_title(allometry_phylo_means_braincase_grp_cat_plot$group)
allometry_phylo_means_rostrum_grp_cat_plot$group <- str_to_title(allometry_phylo_means_rostrum_grp_cat_plot$group)

#Divide by groups - too many lines all together
allometry_phylo_means_group_coeffs_mysticeti <- allometry_phylo_means_group_cat_coeffs %>% filter(group == "Mysticeti")
allometry_phylo_means_group_coeffs_odontoceti <- allometry_phylo_means_group_cat_coeffs %>% filter(group == "Odontoceti")

#Plot allometry regression by category mysticeti
allometry_phylo_means_group_module_2_ggplot_line_mysticeti <- ggplot(allometry_phylo_means_grp_cat_plot_mysticeti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_abline(data = allometry_phylo_means_group_coeffs_mysticeti, 
              aes(intercept = Intercept, slope = Slope,  colour = module, linetype = category), linewidth = 1.2)+
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", fill = "none", linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_means_group_module_2_ggplot_line_mysticeti

#Add phylopic
allometry_phylo_means_group_module_2_ggplot_line_mysticeti <- 
  allometry_phylo_means_group_module_2_ggplot_line_mysticeti +
  add_phylopic(myst, alpha = 1, x = 2.85, y = 0.005, ysize = 0.005, fill = "gray30")
allometry_phylo_means_group_module_2_ggplot_line_mysticeti
#Fix silhouette size

#Plot allometry regression by category odontoceti
allometry_phylo_means_group_module_2_ggplot_line_odontoceti <- ggplot(allometry_phylo_means_grp_cat_plot_odontoceti, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, aes(color = module, fill = module), alpha = 0)+  
  geom_abline(data = allometry_phylo_means_group_coeffs_odontoceti, 
              aes(intercept = Intercept, slope = Slope,  colour = module, linetype = category), linewidth = 1.2)+
  #points after, so they are on top
  scale_linetype_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                        values = c(3,2,4,1))+
  scale_color_manual(name = "Module",values = c(mypalette_paired[2],mypalette_paired[5]), aesthetics = c("colour", "fill"))+
  facet_wrap(vars(module),  scales = "free")+
  theme_bw(base_size = 12)+
  ylab("Regression Score")+
  theme(legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.box = "horizontal",    legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = "none", fill = "none", linetype = guide_legend(override.aes = list(colour = "gray20"), keywidth = unit(3, "char")))
allometry_phylo_means_group_module_2_ggplot_line_odontoceti

#Add phylopic
allometry_phylo_means_group_module_2_ggplot_line_odontoceti <- 
  allometry_phylo_means_group_module_2_ggplot_line_odontoceti +
  add_phylopic(odont, alpha = 1, x = 2.65, y = -0.03, ysize = 0.005, fill = "gray50")
allometry_phylo_means_group_module_2_ggplot_line_odontoceti
#Fix silhouette size

ggarrange(allometry_phylo_means_group_module_2_ggplot_line_mysticeti,allometry_phylo_means_group_module_2_ggplot_line_odontoceti, 
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

####By module and category by group ----
##Plot by module both groups
allometry_phylo_means_braincase_grp_cat_line_ggplot <- ggplot(allometry_phylo_means_braincase_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_means_group_cat_coeffs[c(1:8),], 
              aes(intercept = Intercept, slope = Slope,  colour = category, linetype = group), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  facet_wrap(vars(group))+
  theme_bw(base_size = 12)+
  ylab("Regression Score")+
  ggtitle("Braincase")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = "none")
allometry_phylo_means_braincase_grp_cat_line_ggplot

allometry_phylo_means_rostrum_grp_cat_line_ggplot <- ggplot(allometry_phylo_means_rostrum_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_means_group_cat_coeffs[c(9:16),], 
              aes(intercept = Intercept, slope = Slope,  colour = category, linetype = group), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  facet_wrap(vars(group))+
  theme_bw(base_size = 12)+
  ylab("Regression Score")+
  ggtitle("Rostrum")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = "none")
allometry_phylo_means_rostrum_grp_cat_line_ggplot

ggarrange(allometry_phylo_means_rostrum_grp_cat_line_ggplot, allometry_phylo_means_braincase_grp_cat_line_ggplot,
          ncol = 1, nrow = 2, common.legend = T, legend = "bottom")



#ALLOMETRY ANALYSIS - WHOLE SKULL AND COMPARE PLOTS ROSTRUM AND BRAINCASE ----
##Whole skull ----
##Regression shape on logCS size whole skull
allometry_phylo_means_whole_list <- list()

allometry_phylo_means_whole_list[[1]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_all$size[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list[[2]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_all$size[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list[[3]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_all$size[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list[[4]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_all$size[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_whole <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_whole[[c]] <- summary(allometry_phylo_means_whole_list[[c]])
  names(summary_allometry_phylo_means_whole)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_whole

#Allometry by category by group models
#Combination
allometry_phylo_means_whole_list_comb <- list()

allometry_phylo_means_whole_list_comb[[1]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_all$size[rows_categories_mean_all[[1]]] + gdf_mean_all$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_comb[[2]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_all$size[rows_categories_mean_all[[2]]] + gdf_mean_all$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_comb[[3]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_all$size[rows_categories_mean_all[[3]]] + gdf_mean_all$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_comb[[4]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_all$size[rows_categories_mean_all[[4]]] + gdf_mean_all$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_whole_comb <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_whole_comb[[c]] <- summary(allometry_phylo_means_whole_list_comb[[c]])
  names(summary_allometry_phylo_means_whole_comb)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_whole_comb

#Interaction
allometry_phylo_means_whole_list_int <- list()

allometry_phylo_means_whole_list_int[[1]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[1]]] ~ gdf_mean_all$size[rows_categories_mean_all[[1]]] * gdf_mean_all$group[rows_categories_mean_all[[1]]], phy = trees_list_cat[[1]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_int[[2]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[2]]] ~ gdf_mean_all$size[rows_categories_mean_all[[2]]] * gdf_mean_all$group[rows_categories_mean_all[[2]]], phy = trees_list_cat[[2]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_int[[3]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[3]]] ~ gdf_mean_all$size[rows_categories_mean_all[[3]]] * gdf_mean_all$group[rows_categories_mean_all[[3]]], phy = trees_list_cat[[3]], iter=999, print.progress = TRUE)
allometry_phylo_means_whole_list_int[[4]] <- procD.pgls(gdf_mean_all$coords[,,rows_categories_mean_all[[4]]] ~ gdf_mean_all$size[rows_categories_mean_all[[4]]] * gdf_mean_all$group[rows_categories_mean_all[[4]]], phy = trees_list_cat[[4]], iter=999, print.progress = TRUE)

#Main results of ANOVA analysis of allometry with logCS
summary_allometry_phylo_means_whole_int <- list()

#Loop
for (c in 1:length(categories_list)){
  summary_allometry_phylo_means_whole_int[[c]] <- summary(allometry_phylo_means_whole_list_int[[c]])
  names(summary_allometry_phylo_means_whole_int)[[c]] <- categories_list[[c]]
}

summary_allometry_phylo_means_whole_int

#Save results of significant regression to file
sink("Output/8-Allometry means/allometry_phylo_means_shape_size_grp_cat_wholeskull.txt")
print("Null")
summary_allometry_phylo_means_whole

print("Combination +")
summary_allometry_phylo_means_whole_comb

print("Interaction *")
summary_allometry_phylo_means_whole_int
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_phylo_means_models_grp_cat_whole_list <- list()

#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_whole_list[[c]] <- anova(allometry_phylo_means_whole_list[[c]], allometry_phylo_means_whole_list_comb[[c]], allometry_phylo_means_whole_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_whole_list)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_whole_list

#ANOVAs - is a model significantly better than the others? - check between comb and int only
anova_allometry_phylo_means_models_grp_cat_whole_list_1 <- list()

#Loop
for (c in 1:length(categories_list)){  
  anova_allometry_phylo_means_models_grp_cat_whole_list_1[[c]] <- anova(allometry_phylo_means_whole_list_comb[[c]], allometry_phylo_means_whole_list_int[[c]])
  names(anova_allometry_phylo_means_models_grp_cat_whole_list_1)[[c]] <- categories_list[[c]]
}

anova_allometry_phylo_means_models_grp_cat_whole_list_1

#In all categories, comb better than null but interaction not difference - no sign diff in slope between groups

##Extract RegScores best model for comparison
#Regression score of shape vs logCS - regression method with "RegScore" plotting

allometry_phylo_means_whole_plot_regscore <- list()

for (c in 1:length(categories_list)){  
  allometry_phylo_means_whole_plot_regscore[[c]] <- plot(allometry_phylo_means_whole_list_comb[[c]], type = "regression",predictor = gdf_mean_all$size[rows_categories_mean_all[[c]]], reg.type = "RegScore", 
                                                           main = paste("Shape vs logCS by group", "-", levels(categories)[[c]]), xlab = "logCS", pch = 21, cex = 1.2, font.main = 2,
                                                           col = mypalette_category[c], bg = mypalette_category[c])
  text(x = gdf_mean_all$size[rows_categories_mean_all[[c]]], y = allometry_phylo_means_whole_plot_regscore[[c]]$RegScore, labels = gdf_mean_all$genus[rows_categories_mean_all[[c]]],
       pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels  
}

#Create object to use for linear model
allometry_phylo_means_whole_regscores <- rbind(allometry_phylo_means_whole_plot_regscore[[1]][["RegScore"]], allometry_phylo_means_whole_plot_regscore[[2]][["RegScore"]] ,
                                                 allometry_phylo_means_whole_plot_regscore[[3]][["RegScore"]] ,allometry_phylo_means_whole_plot_regscore[[4]][["RegScore"]] )

#Linear model for line by group and category 
allometry_phylo_means_whole_regscores_df <- data.frame(RegScores = allometry_phylo_means_whole_regscores, logCS = gdf_mean_all$size, genus = gdf_mean_all$genus, group = gdf_mean_all$group,
                                                         category = gdf_mean_all$category)

allometry_phylo_means_whole_regscores_df$grp_cat <- interaction(allometry_phylo_means_whole_regscores_df$group, allometry_phylo_means_whole_regscores_df$category)



###Pairwise comparison of regression model between groups and category ----
#Create models, with different slopes and int or just int
allometry_phylo_means_whole_null_cat <- lm.rrpp(RegScores ~ logCS * category,   #null model only differences by category
                                                  data = allometry_phylo_means_whole_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_whole_null_grp <- lm.rrpp(RegScores ~ logCS * group,   #null model only differences by group
                                                  data = allometry_phylo_means_whole_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_whole_comb <- lm.rrpp(RegScores ~ logCS * category + group,
                                              data = allometry_phylo_means_whole_regscores_df, print.progress = FALSE, iter = 999) 
allometry_phylo_means_whole_int <- lm.rrpp(RegScores ~ logCS * category * group,
                                             data = allometry_phylo_means_whole_regscores_df, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_phylo_means_whole_null_cat)
summary(allometry_phylo_means_whole_null_grp)
summary(allometry_phylo_means_whole_comb)
summary(allometry_phylo_means_whole_int)

#Anova for difference between models
anova_allometry_phylo_means_models_grp_cat_whole <- anova(allometry_phylo_means_whole_null_cat,allometry_phylo_means_whole_null_grp, allometry_phylo_means_whole_comb, allometry_phylo_means_whole_int)
anova_allometry_phylo_means_models_grp_cat_whole

#Check interactions in models
anova(allometry_phylo_means_whole_null_cat)
anova(allometry_phylo_means_whole_null_grp)
anova(allometry_phylo_means_whole_comb)
anova(allometry_phylo_means_whole_int)

#No sign interaction between category and group. Group and category don't influence regression

#Calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_phylo_means_whole_int1 <- lm(RegScores ~ logCS * grp_cat,
                                         data = allometry_phylo_means_whole_regscores_df) 
#Check anova still ok
anova(allometry_phylo_means_whole_int1)

#Get pairwise comparisons of slopes
allometry_phylo_means_whole_emms <- emmeans(allometry_phylo_means_whole_int1, "grp_cat")

#to make graph, confusing for lots of groups - pwpp(allometry_phylo_means_whole_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_phylo_means_whole_ems_table <- pwpm(allometry_phylo_means_whole_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
allometry_phylo_means_whole_ems_table

#Sign differences between Odontoceti early and other categories of Odont.

#Save results to file
sink("Output/8-Allometry means/pairwise_allometry_phylo_means_whole.txt")
print("ANOVA models")
print(anova_allometry_phylo_means_models_grp_cat_whole)

print("summary models - lmrpp")
anova(allometry_phylo_means_whole_null_cat)
anova(allometry_phylo_means_whole_null_grp)
anova(allometry_phylo_means_whole_comb)
anova(allometry_phylo_means_whole_int)

print("summary model used for comparisons - lm")
summary(allometry_phylo_means_whole_int1)
anova(allometry_phylo_means_whole_int1)

print("Pairwise comparison using emmeans")
summary(allometry_phylo_means_whole_emms)

print("Full results table emmeans pairwise comparions")
pwpm(allometry_phylo_means_whole_emms)
sink()

####Plot allometry whole skull by category by group ----
#Create data frame object that ggplot can read - use data from plot object you want to improve
#Convert data frame to tibble
allometry_phylo_means_whole_grp_cat_plot <- as_tibble(allometry_phylo_means_whole_regscores_df)

#Plot allometry regression by category by group 
allometry_phylo_means_whole_grp_cat_ggplot <- ggplot(allometry_phylo_means_whole_grp_cat_plot, aes(x = logCS, y = RegScores))+
  geom_point(aes(colour = category,fill = category), size = 0, alpha = 0)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group,colour = category,fill = category, group = grp_cat), inherit.aes = F,        
              se = F, linewidth = 1.5, alpha = 1)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult"), 
                      values = mypalette_category, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Groups", labels =  levels(groups),
                        values = c(1,6))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle("Whole skull")+
  theme(plot.title =  element_text(face = 3, hjust = 0.5, size = 16), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0), legend.box = "vertical")+
  guides(colour = guide_legend(keywidth = unit(4, "char"), nrow = 1, byrow = F, override.aes = list(alpha = 1, size = 0, linetype =5, colour = mypalette_category, fill = mypalette_category)),
         linetype = guide_legend(keywidth = unit(4, "char"), override.aes = list(colour = c("grey30","gray50"))))
allometry_phylo_means_whole_grp_cat_ggplot

#Add phylopic
allometry_phylo_means_whole_grp_cat_ggplot <- 
  allometry_phylo_means_whole_grp_cat_ggplot +
  add_phylopic(myst, alpha = 1, x = 3.1, y = -0.025, ysize = 0.004, fill = "gray30")+
  add_phylopic(odont, alpha = 1, x = 2.5, y = 0.005, ysize = 0.003, fill = "gray50")
allometry_phylo_means_whole_grp_cat_ggplot

ggarrange(allometry_phylo_means_whole_grp_cat_ggplot, allometry_phylo_means_rostrum_grp_cat_ggplot, allometry_phylo_means_braincase_grp_cat_ggplot, 
          nrow = 1, ncol = 3, common.legend = T, legend = "bottom")

####Line plot by group and module whole skull ----
#Whole skull

allometry_phylo_means_whole_grp_cat_plot$module <- "Whole skull"

#Linear model for line by group and category
allometry_phylo_means_group_regline_whole <- lm(RegScores ~ logCS * group, data = allometry_phylo_means_whole_grp_cat_plot)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_phylo_means_group_regline_coeffs_whole <- as.matrix(allometry_phylo_means_group_regline_whole$coefficients)

#Save intercepts and slopes separately
allometry_phylo_means_group_regline_intercepts_whole <- as.matrix(allometry_phylo_means_group_regline_coeffs_whole[c(1, 3:(length(levels(groups))+1)),])
allometry_phylo_means_group_regline_slopes_whole <- as.matrix(allometry_phylo_means_group_regline_coeffs_whole[c(2, length(levels(groups))+2:(length(levels(groups)))),])

#Calculate real intercepts and slopes
allometry_phylo_means_group_regline_intercepts_ok_whole <- as.matrix(c(allometry_phylo_means_group_regline_intercepts_whole[1,], allometry_phylo_means_group_regline_intercepts_whole[1,]+
                                                                             allometry_phylo_means_group_regline_intercepts_whole[2:length(allometry_phylo_means_group_regline_intercepts_whole),]))

allometry_phylo_means_group_regline_slopes_ok_whole <- as.matrix(c(allometry_phylo_means_group_regline_slopes_whole[1,], allometry_phylo_means_group_regline_slopes_whole[1,]+
                                                                         allometry_phylo_means_group_regline_slopes_whole[2:length(allometry_phylo_means_group_regline_slopes_whole),]))

#Save as data frame with grouping variables
allometry_phylo_means_group_coeffs_whole <- data.frame(Slope = allometry_phylo_means_group_regline_slopes_ok_whole, Intercept = allometry_phylo_means_group_regline_intercepts_ok_whole, 
                                                           row.names = levels(groups))
#Check for NA and other issues
allometry_phylo_means_group_coeffs_whole 

#Add classifiers
allometry_phylo_means_group_coeffs_whole <- allometry_phylo_means_group_coeffs_whole %>% mutate(group = levels(groups),
                                                                                                        module = 'Whole skull')
allometry_phylo_means_group_coeffs_whole

#Put together
allometry_phylo_means_group_coeffs_3 <- rbind(allometry_phylo_means_group_coeffs_whole, allometry_phylo_means_group_coeffs_braincase, allometry_phylo_means_group_coeffs_rostrum)
allometry_phylo_means_group_coeffs_3

#Create common tibble for plotting
allometry_phylo_means_grp_cat_plot_3 <- rbind(allometry_phylo_means_whole_grp_cat_plot, allometry_phylo_means_grp_cat_plot)
glimpse(allometry_phylo_means_grp_cat_plot_3)

##Plot by module both groups
allometry_phylo_means_grp_line_3_ggplot <- ggplot(allometry_phylo_means_grp_cat_plot_3, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+  
  #line on plot
  geom_abline(data = allometry_phylo_means_group_coeffs_3, 
              aes(intercept = Intercept, slope = Slope,  colour = group, linetype = module), linewidth = 1.2)+
  #points after, so they are on top
  scale_colour_manual(name = "Groups", labels =  levels(groups), 
                      values = mypalette_groups, aesthetics = c("colour","fill"))+
  scale_linetype_manual(name = "Modules",
                        values = c(1,2,4))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  theme(legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom",  legend.direction = "horizontal", legend.justification = c(0,0))+
  guides(colour = "none",
         linetype = guide_legend(keywidth = unit(3, "char"), override.aes = list(colour = c("gray20"))))
allometry_phylo_means_grp_line_3_ggplot

#Add silhouettes groups
allometry_phylo_means_grp_line_3_ggplot  <- 
  allometry_phylo_means_grp_line_3_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 2.95, y = -0.04, ysize = 0.0075, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 2.7, y = 0.02, ysize = 0.0065, fill = mypalette_groups[2])
allometry_phylo_means_grp_line_3_ggplot
