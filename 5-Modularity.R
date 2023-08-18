#===========================================================#
#                                                           #
#     SKULL MODULARITY - MYSTICETI & ODONTOCETI   #
#                                                           #
#===========================================================#

#CH. 5 - Modularity analyses

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
library(paleomorph)
library(qgraph)
library(rray)
library(abind)
library(reshape2)
library(scales)
library(mcp)
library(tidytext)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

# RPushbullet::pbPost(type = "note", title = "R process done")

#apropos("x") lists objects with matching part of name

#MODULARITY TEST - RAW DATA AND ALLOMETRY CORRECTED ----

##Define modules ----

#Change modules for different hypothesis
modules_11 <- modules_all

modules_11[premaxilla]<-1 
modules_11[maxilla]<-2 
modules_11[nasals]<-3 
modules_11[orbit]<-4 
modules_11[squamosal]<-5
modules_11[interparietal]<-6
modules_11[supraoccipital]<-7 
modules_11[exoccipital]<-8
modules_11[condyles]<-9
modules_11[basioccipital]<-10 
modules_11[palatine]<-11 
modules_11

#Change modules for different hypothesis
#Reduced 5 modules Churchill et al. (2019)-extra module added for premax/rostrum for curves, best for dolphins
modules_5_C19 <- modules_all

#Anterior rostrum - premaxilla
modules_5_C19[my_curves$Curve.in$SC1]<-1 
modules_5_C19[my_curves$Curve.in$SC6]<-1 
modules_5_C19[my_curves$Curve.in$SC11]<-1 
modules_5_C19[my_curves$Curve.in$SC12]<-1 
modules_5_C19[c(1,10,20,19,11,2,5,14)]<-1
#Palate - ventral anterior rostrum and palatine
modules_5_C19[my_curves$Curve.in$SC29]<-2 
modules_5_C19[my_curves$Curve.in$SC30]<-2 
modules_5_C19[my_curves$Curve.in$SC31]<-2 
modules_5_C19[my_curves$Curve.in$SC32]<-2 
modules_5_C19[c(43,48,49,44,51,52,46,47)]<-2 
#Face - maxilla ascending, orbit, jugal, side of exoccipital
modules_5_C19[my_curves$Curve.in$SC2]<-3 
modules_5_C19[my_curves$Curve.in$SC3]<-3 
modules_5_C19[my_curves$Curve.in$SC7]<-3 
modules_5_C19[my_curves$Curve.in$SC8]<-3 
modules_5_C19[orbit]<-3
modules_5_C19[c(3,4,12,13)]<-3
modules_5_C19[my_curves$Curve.in$SC23]<-3 
modules_5_C19[my_curves$Curve.in$SC26]<-3 
modules_5_C19[c(35,27,39,34)]<-3
#Zygomatic - squamosal, basioccipital and vomer
modules_5_C19[squamosal]<-4 
modules_5_C19[basioccipital]<-4 
modules_5_C19[c(64)]<-4
#Nasals
modules_5_C19[nasals]<-5 
#Vault - interparietal, supraoccipital and condyles
modules_5_C19[interparietal]<-6
modules_5_C19[supraoccipital]<-6 
modules_5_C19[condyles]<-6
modules_5_C19

#Goswami Therian 6 modules (from Churchill et al. (2019))-best with symmetric component dolphins
modules_6_th <- modules_all

#Vault - interparietal, supraoccipital and condyles
modules_6_th[interparietal]<-1
modules_6_th[my_curves$Curve.in$SC41]<-1 
modules_6_th[my_curves$Curve.in$SC42]<-1
modules_6_th[c(60,61,62)]<-1
#Basicranium - condyle, midline supraoccipital, exoccipital
modules_6_th[my_curves$Curve.in$SC43]<-2
modules_6_th[c(63)]<-2
modules_6_th[condyles]<-2
modules_6_th[exoccipital]<-2
#Anterior oral-nasal - nasals, premaxilla and ascending maxilla back
modules_6_th[nasals]<-3
modules_6_th[my_curves$Curve.in$SC11]<-3
modules_6_th[my_curves$Curve.in$SC12]<-3
modules_6_th[my_curves$Curve.in$SC3]<-3
modules_6_th[my_curves$Curve.in$SC8]<-3
modules_6_th[c(13,4,5,14,19,20)]<-3
#Orbit - orbit and jugal
modules_6_th[orbit]<-4
#Molar - maxilla dorsal and ventral, palatine, vomer
modules_6_th[my_curves$Curve.in$SC1]<-5
modules_6_th[my_curves$Curve.in$SC6]<-5
modules_6_th[my_curves$Curve.in$SC2]<-5
modules_6_th[my_curves$Curve.in$SC7]<-5
modules_6_th[c(1,10,11,2,3,12)]<-5
modules_6_th[my_curves$Curve.in$SC29]<-5 
modules_6_th[my_curves$Curve.in$SC30]<-5 
modules_6_th[my_curves$Curve.in$SC31]<-5 
modules_6_th[my_curves$Curve.in$SC32]<-5 
modules_6_th[c(43,48,49,44,51,52,46,47,64)]<-5 
#Zygomatic-Pterygoid - basioccipital and squamosal
modules_6_th[which(modules_6_th %in% c("squamosal", "basioccipital"))]<-6
modules_6_th

#DelCastillo 3 modules (from Churchill et al. (2019))
modules_3_DCA <- modules_all

#Rostrum - maxilla, premaxilla, palatine, vomer, nasals
modules_3_DCA[nasals]<-1
modules_3_DCA[maxilla]<-1
modules_3_DCA[premaxilla]<-1
modules_3_DCA[palatine]<-1
modules_3_DCA[c(64)]<-1
#Neurocranium - supraoccipital, exoccipital, interparietal, condyle, orbit
modules_3_DCA[supraoccipital]<-2
modules_3_DCA[exoccipital]<-2
modules_3_DCA[interparietal]<-2
modules_3_DCA[condyles]<-2
modules_3_DCA[orbit]<-2
#Basicranium - basioccipital and squamosal
modules_3_DCA[basioccipital]<-3
modules_3_DCA[squamosal]<-3
modules_3_DCA

#Drake&Klingenberg 2 modules (from Churchill et al. (2019))
modules_2_DK <- modules_all

#Facial - maxilla, premaxilla, palatine, vomer, nasals
modules_2_DK[nasals]<-1
modules_2_DK[maxilla]<-1
modules_2_DK[premaxilla]<-1
modules_2_DK[palatine]<-1
modules_2_DK[orbit]<-1
modules_2_DK[c(64)]<-1
#Neurocranial - supraoccipital, exoccipital, interparietal, condyle, squamosal, basioccipital
modules_2_DK[supraoccipital]<-2
modules_2_DK[exoccipital]<-2
modules_2_DK[interparietal]<-2
modules_2_DK[condyles]<-2
modules_2_DK[squamosal]<-2
modules_2_DK[basioccipital]<-2
modules_2_DK

#Main hypothesis
#Developmental modules neural crest and mesoderm 2 modules (from Goswami et al. (2022))
modules_2_dev <- modules_all

#Neural crest - maxilla, premaxilla, palatine, vomer, nasals, squamosal - ROSTRUM
modules_2_dev[nasals]<-1
modules_2_dev[maxilla]<-1
modules_2_dev[premaxilla]<-1
modules_2_dev[palatine]<-1
modules_2_dev[orbit]<-1
modules_2_dev[c(64)]<-1
modules_2_dev[squamosal]<-1
#Mesoderm - supraoccipital, exoccipital, interparietal, condyle, basioccipital - BRAINCASE
modules_2_dev[supraoccipital]<-2
modules_2_dev[exoccipital]<-2
modules_2_dev[interparietal]<-2
modules_2_dev[condyles]<-2
modules_2_dev[basioccipital]<-2
modules_2_dev

#Save modules to file
modules_df <- data.frame(lm = c(1:dim(gdf$coords)[[1]]), modules_all = modules_all,
                             m11 = modules_11, mR5 = modules_5_C19, m6 = modules_6_th, m3 = modules_3_DCA, m2 = modules_2_DK,
                             m2_dev = modules_2_dev)
write_csv(modules_df, "Output/modules_all.csv", col_names = T)

 
#COMPARE MODULARITY HYPOTHESES - compare.CR ----

##All data ----
#Perform modularity test all modules - compare selected modules to the null assumption of random assignment of partitions (no modularity at all)
#Run CR analysis for each module hypothesis
m11_all <- modularity.test(gdf$coords, modules_df$m11, CI = F, iter = 999, print.progress = T)
mR5_all <- modularity.test(gdf$coords, modules_df$mR5, CI = F, iter = 999, print.progress = T)
m6_all <- modularity.test(gdf$coords, modules_df$m6, CI = F, iter = 999, print.progress = T)
m3_all <- modularity.test(gdf$coords, modules_df$m3, CI = F, iter = 999, print.progress = T)
m2_all <- modularity.test(gdf$coords, modules_df$m2, CI = F, iter = 999, print.progress = T)
m2dev_all <- modularity.test(gdf$coords, modules_df$m2_dev, CI = F, iter = 999, print.progress = T)

#Compare strength of modularity with CR test
CR_compare_all <- compare.CR(m11_all, mR5_all, m6_all, m3_all, m2_all,m2dev_all, CR.null = TRUE)

#The best supported hypothesis is the one with the smallest effect size in summary()
summary(CR_compare_all)

#Make graph to check which module has the lowest Z score
CR_compare_all_df <- data.frame(modules = names(CR_compare_all[["sample.z"]]), z = CR_compare_all[["sample.z"]], 
                                se = CR_compare_all[["sample.se"]])
CR_compare_all_df <- CR_compare_all_df[-1,] #remove null hypothesis

#Make factor for variable
CR_compare_all_df$modules <- factor(CR_compare_all_df$modules, 
                                    levels = c("m11_all", "mR5_all", "m6_all","m3_all","m2_all", "m2dev_all")) #copy from string printed with the code above
#Order
CR_compare_all_df <- CR_compare_all_df[order(CR_compare_all_df$modules),]
CR_compare_all_df

#Modules list names for plot
module_hyp_list <- c("All sep. (11 mod.)", "5 mod.", "6 mod.","3 mod.","2 mod.", "2 mod. dev.")

CR_compare_all_plot <- ggplot(CR_compare_all_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = 2, colour = mypalette_paired[12], stroke = 1.5)+
  ggtitle("Cetacea")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_all_plot

#Save best model data
best_model_mod_all <- CR_compare_all_df %>% slice_min(n = 1, z)
best_model_mod_all

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/Compare_CR_all.txt")
summary(CR_compare_all)

print("Best modularity model")
summary(m2dev_all)
print("Pairwise CR best model")
m2dev_all[["CR.mat"]]
sink() 

#Histogram plot best model
plot(m2dev_all)

##Mysticeti ----

#Find rows for Mysticeti
rows_mysticeti <- which(gdf$group == "mysticeti")

#Run CR analysis for each module hypothesis
m11_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$m11, CI = F, iter = 999, print.progress = T)
mR5_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$mR5, CI = F, iter = 999, print.progress = T)
m6_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$m6, CI = F, iter = 999, print.progress = T)
m3_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$m3, CI = F, iter = 999, print.progress = T)
m2_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$m2, CI = F, iter = 999, print.progress = T)
m2_dev_myst <- modularity.test(gdf$coords[,,rows_mysticeti], modules_df$m2_dev, CI = F, iter = 999, print.progress = T)

#Compare strength of modularity with CR test
CR_compare_myst <- compare.CR(m11_myst, mR5_myst, m6_myst, m3_myst, m2_myst, m2_dev_myst, CR.null = TRUE)

#The best supported hypothesis is the one with the smallest effect size in summary()
summary(CR_compare_myst)

#Make graph to check which module has the lowest Z score
CR_compare_myst_df <- data.frame(modules = names(CR_compare_myst[["sample.z"]]), z = CR_compare_myst[["sample.z"]], 
                                 se = CR_compare_myst[["sample.se"]])
CR_compare_myst_df <- CR_compare_myst_df[-1,] #remove null hypothesis

#Make factor for variable
CR_compare_myst_df$modules <- factor(CR_compare_myst_df$modules, 
                                     levels = c("m11_myst", "mR5_myst", "m6_myst","m3_myst","m2_myst", "m2_dev_myst")) #copy from string printed with the code above
#Order
CR_compare_myst_df <- CR_compare_myst_df[order(CR_compare_myst_df$modules),]
CR_compare_myst_df

CR_compare_myst_plot <- ggplot(CR_compare_myst_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = shapes[1], colour = mypalette_groups[1], stroke = 1.5)+
  ggtitle("Mysticeti")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_myst_plot

#Save best model data
best_model_mod_myst <- CR_compare_myst_df %>% slice_min(n = 1, z)
best_model_mod_myst

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/Compare_CR_myst.txt")
summary(CR_compare_myst)

print("Best modularity model")
summary(m2_myst)
print("Pairwise CR best model")
m2_myst[["CR.mat"]]
sink() 

#Histogram plot best model
plot(m2_myst)

##Odontoceti ----

#Find rows for Odontoceti
rows_odontoceti <- which(gdf$group == "odontoceti")

#Run CR analysis for each module hypothesis
m11_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$m11, CI = F, iter = 999, print.progress = T)
mR5_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$mR5, CI = F, iter = 999, print.progress = T)
m6_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$m6, CI = F, iter = 999, print.progress = T)
m3_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$m3, CI = F, iter = 999, print.progress = T)
m2_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$m2, CI = F, iter = 999, print.progress = T)
m2_dev_odont <- modularity.test(gdf$coords[,,rows_odontoceti], modules_df$m2_dev, CI = F, iter = 999, print.progress = T)

#Compare strength of modularity with CR test
CR_compare_odont <- compare.CR(m11_odont, mR5_odont, m6_odont, m3_odont, m2_odont, m2_dev_odont, CR.null = TRUE)

#The best supported hypothesis is the one with the smallest effect size in summary()
summary(CR_compare_odont)

#Make graph to check which module has the lowest Z score
CR_compare_odont_df <- data.frame(modules = names(CR_compare_odont[["sample.z"]]), z = CR_compare_odont[["sample.z"]], 
                                  se = CR_compare_odont[["sample.se"]])
CR_compare_odont_df <- CR_compare_odont_df[-1,] #remove null hypothesis

#Make factor for variable
CR_compare_odont_df$modules <- factor(CR_compare_odont_df$modules, 
                                      levels = c("m11_odont", "mR5_odont", "m6_odont","m3_odont","m2_odont","m2_dev_odont")) #copy from string printed with the code above
#Order
CR_compare_odont_df <- CR_compare_odont_df[order(CR_compare_odont_df$modules),]
CR_compare_odont_df

CR_compare_odont_plot <- ggplot(CR_compare_odont_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = shapes[2], colour = mypalette_groups[2], stroke = 1.5)+ #choose slightly darker shade of green to make sure it shows
  ggtitle("Odontoceti")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_odont_plot

ggarrange(CR_compare_all_plot, CR_compare_myst_plot, CR_compare_odont_plot, nrow=1, ncol=3)

#Save best model data
best_model_mod_odont <- CR_compare_odont_df %>% slice_min(n = 1, z)
best_model_mod_odont

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/Compare_CR_odont.txt")
summary(CR_compare_odont)

print("Best modularity model")
summary(m2_dev_odont)
print("Pairwise CR best model")
m2_dev_odont[["CR.mat"]]
sink() 

#Histogram plot best model
plot(m2_dev_odont)

###Bar plot to compare Z-scores modularity hypothesis by group ----

#Create combined data frame
CR_compare_all_df$data <- "All"
CR_compare_myst_df$data <- levels(groups)[1]
CR_compare_odont_df$data <- levels(groups)[2]

CR_compare_plot_df <- bind_rows(CR_compare_all_df, CR_compare_myst_df, CR_compare_odont_df)

#Add column with proper mod. hypothesis labels
CR_compare_plot_df$labels <- rep(module_hyp_list, times = 3)
CR_compare_plot_df

#Arrange data by z-score in each group
CR_compare_plot_df <- CR_compare_plot_df %>% 
  arrange(data, z) 
CR_compare_plot_df

#Make vector to fix labels for plot
sub_labels_string <- paste0("___", levels(as.factor(CR_compare_plot_df$data)))

#Create a vertical bar plot faceted by group
CR_compare_bar_plot <- ggplot(CR_compare_plot_df, aes(x = reorder_within(labels, z, data), y = z, fill = data)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8, color = "gray50") +
  facet_wrap(vars(data), scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c(mypalette_paired[12],mypalette_groups))+   
  scale_fill_manual(values = c(mypalette_paired[12],mypalette_groups))+   
  labs(y = "Z-Scores") +
  theme_minimal(base_size = 13)+
  scale_x_discrete(labels = function(x) {
    x <- gsub(sub_labels_string[1], "", x)
    x <- gsub(sub_labels_string[2], "", x)
    x <- gsub(sub_labels_string[3], "", x)
    return(x)
  })+ #Function to delete faceting from labels
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15), axis.line.y = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(),
        strip.text.x = element_text(size=12), legend.position = "none",
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"))
CR_compare_bar_plot

###Heatmaps plots for significant differences between modularity hypothesis by group ----
#####Functions----
{#Get lower triangle of the correlation matrix
  get_lower_tri<-function(x){
    x[upper.tri(x)] <- NA
    return(x)
  }
  #Get upper triangle of the correlation matrix
  get_upper_tri <- function(x){
    x[lower.tri(x)]<- NA
    return(x)
  }
  #Reorder table
  reorder_corr_table <- function(x){
    # Use correlation between variables as distance
    dd <- as.dist((1-x)/2)
    hc <- hclust(dd)
    x <-x[hc$order, hc$order]
  }
}


#Create palette for comparison between modules and groups
mypalette_seq_modularity <- brewer.pal(9,"YlGnBu")
image(1:9,1, as.matrix(1:9), col = mypalette_seq_modularity,xlab="Yel-Blue-Green (sequential)",
      ylab = "", yaxt = "n")

#All data
#Save p-values as object
modularity_z_all <- CR_compare_all[["pairwise.z"]]
modularity_pvals_all <- CR_compare_all[["pairwise.P"]]

#Remove null hypothesis
modularity_z_all <- modularity_z_all[-1,]
modularity_z_all <- modularity_z_all[,-1]
modularity_pvals_all <- modularity_pvals_all[-1,]
modularity_pvals_all <- modularity_pvals_all[,-1]  

#Set correct row and col names for both
rownames(modularity_z_all) <- module_hyp_list
rownames(modularity_pvals_all) <- module_hyp_list
colnames(modularity_z_all) <- module_hyp_list
colnames(modularity_pvals_all) <- module_hyp_list

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_all_upper_tri <- get_upper_tri(modularity_z_all)
modularity_pvals_all_upper_tri <- get_upper_tri(modularity_pvals_all)

#Melt to make table in the format needed for heatmap
modularity_z_all_melt <- melt(modularity_z_all_upper_tri, value.name = "z", na.rm = TRUE)
modularity_pvals_all_melt <- melt(modularity_pvals_all_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
modularity_pvals_all_melt$z <- modularity_z_all_melt$z

#Create columns where only significant values are shown
modularity_pvals_all_melt <- modularity_pvals_all_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          z_if_sig = ifelse(sig_p, z, NA)) %>%
  mutate_at(vars(starts_with("z")), list(~ round(., 2)))
modularity_pvals_all_melt 

#Nice heatmap plot
modularity_pvals_all_heatmap_ggplot <- ggplot(data = modularity_pvals_all_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_all_melt$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("All data")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
modularity_pvals_all_heatmap_ggplot 

#Mysticeti
#Save p-values as object
modularity_z_myst <- CR_compare_myst[["pairwise.z"]]
modularity_pvals_myst <- CR_compare_myst[["pairwise.P"]]

#Remove null hypothesis
modularity_z_myst <- modularity_z_myst[-1,]
modularity_z_myst <- modularity_z_myst[,-1]
modularity_pvals_myst <- modularity_pvals_myst[-1,]
modularity_pvals_myst <- modularity_pvals_myst[,-1]  

#Set correct row and col names for both
rownames(modularity_z_myst) <- module_hyp_list
rownames(modularity_pvals_myst) <- module_hyp_list
colnames(modularity_z_myst) <- module_hyp_list
colnames(modularity_pvals_myst) <- module_hyp_list

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_myst_upper_tri <- get_upper_tri(modularity_z_myst)
modularity_pvals_myst_upper_tri <- get_upper_tri(modularity_pvals_myst)

#Melt to make table in the format needed for heatmap
modularity_z_myst_melt <- melt(modularity_z_myst_upper_tri, value.name = "z", na.rm = TRUE)
modularity_pvals_myst_melt <- melt(modularity_pvals_myst_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
modularity_pvals_myst_melt$z <- modularity_z_myst_melt$z

#Create columns where only significant values are shown
modularity_pvals_myst_melt <- modularity_pvals_myst_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                  z_if_sig = ifelse(sig_p, z, NA)) %>%
  mutate_at(vars(starts_with("z")), list(~ round(., 2)))
modularity_pvals_myst_melt 

#Nice heatmap plot
modularity_pvals_myst_heatmap_ggplot <- ggplot(data = modularity_pvals_myst_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_myst_melt$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Mysticeti")+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = "none")
modularity_pvals_myst_heatmap_ggplot 

#Odontoceti
#Save p-values as object
modularity_z_odont <- CR_compare_odont[["pairwise.z"]]
modularity_pvals_odont <- CR_compare_odont[["pairwise.P"]]

#Remove null hypothesis
modularity_z_odont <- modularity_z_odont[-1,]
modularity_z_odont <- modularity_z_odont[,-1]
modularity_pvals_odont <- modularity_pvals_odont[-1,]
modularity_pvals_odont <- modularity_pvals_odont[,-1]  

#Set correct row and col names for both
rownames(modularity_z_odont) <- module_hyp_list
rownames(modularity_pvals_odont) <- module_hyp_list
colnames(modularity_z_odont) <- module_hyp_list
colnames(modularity_pvals_odont) <- module_hyp_list

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_odont_upper_tri <- get_upper_tri(modularity_z_odont)
modularity_pvals_odont_upper_tri <- get_upper_tri(modularity_pvals_odont)

#Melt to make table in the format needed for heatmap
modularity_z_odont_melt <- melt(modularity_z_odont_upper_tri, value.name = "z", na.rm = TRUE)
modularity_pvals_odont_melt <- melt(modularity_pvals_odont_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
modularity_pvals_odont_melt$z <- modularity_z_odont_melt$z

#Create columns where only significant values are shown
modularity_pvals_odont_melt <- modularity_pvals_odont_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                  z_if_sig = ifelse(sig_p, z, NA)) %>%
  mutate_at(vars(starts_with("z")), list(~ round(., 2)))
modularity_pvals_odont_melt 

#Nice heatmap plot
modularity_pvals_odont_heatmap_ggplot <- ggplot(data = modularity_pvals_odont_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_odont_melt$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
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
modularity_pvals_odont_heatmap_ggplot 

plotM1 <- ggarrange(modularity_pvals_all_heatmap_ggplot , modularity_pvals_myst_heatmap_ggplot ,
                    modularity_pvals_odont_heatmap_ggplot ,
                    ncol = 3, nrow = 1, common.legend = F)
plotM1 <- annotate_figure(plotM1, top = text_grob("Modularity hypothesis comparison (Z-scores difference)", 
                                                face = "bold", size = 16, just = c(0.5,6)))
plotM1

##Plot modules on surfaces ----

col_modules_2 <-  as.factor(modules_2_DK)

levels(col_modules_2) <- c("tomato","turquoise4" )


col_modules_2dev <-  as.factor(modules_2_dev)

levels(col_modules_2dev) <- c("orange2","darkblue" )

#All data+Odont
#Plot on surface
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[,,41], col = col_modules_2dev, type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/all_odont_modules.png") 
rgl.snapshot(filename = "Output/all_odont_modules1.png")
rgl.snapshot(filename = "Output/all_odont_modules2.png")
play3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6, movie = "all_odont_modules" ,dir = "Output/")
clear3d()

#Myst
shade3d(myst_fetus, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Ff3",Ids)], col =  col_modules_2, type = "s",
          radius = 1, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/myst_modules.png") 
rgl.snapshot(filename = "Output/myst_modules1.png") 
rgl.snapshot(filename = "Output/myst_modules2.png") 
play3d(spin3d(axis = c(0, 1,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 1,1), rpm = 10), duration = 6, movie = "myst_modules" ,dir = "Output/")

###### 

#Next - ch. 6 - Disparity analyses
