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
  geom_point(size = 6, shape = shapes[2], colour = mypalette_paired[4], stroke = 1.5)+ #choose slightly darker shade of green to make sure it shows
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

#HERE####
##Bar plot for all, Myst and Odont

CR_compare_all_df$data <- "all"
CR_compare_myst_df$data <- "myst"
CR_compare_odont_df$data <- "odont"

CR_compare_plot_df <- bind_rows(CR_compare_all_df, CR_compare_myst_df, CR_compare_odont_df)

CR_compare_plot_df <- CR_compare_plot_df %>% 
  arrange(data, z) 

# Create a vertical bar plot faceted by data with inverted axis
ggplot(CR_compare_plot_df, aes(x = reorder_within(modules, -z, data), y = z, fill = data)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_errorbar(aes(ymin = z - se, ymax = z + se), width = 0.4, position = position_dodge(width = 0.5)) +
  facet_wrap(~ data, scales = "free_x", ncol = 3) +
  labs(title = "Z Values with Confidence Intervals",
       x = "Modules",
       y = "Z Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank()) +
  scale_y_reverse()


p


###Plot modules on surfaces ----

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
