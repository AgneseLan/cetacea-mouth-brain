#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH. 4b - Modularity analyses - mean shapes phylogenetic correction

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
library(paleomorph)
library(qgraph)
library(abind)
library(reshape2)
library(scales)
library(mcp)
library(tidytext)
library(ggpubr)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)

#RPushbullet::pbPost(type = "note", title = "R process done")

#apropos("x") lists objects with matching part of name

#MODULARITY TEST - PHYLOGENETICALLY CORRECTED ----

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
write_csv(modules_df, "Output/4b-Modularity phylo means/modules_all.csv", col_names = T)

##Plot modules on surfaces ----

col_modules_11 <-  col_modules

col_modules_5 <-  as.factor(modules_5_C19)

levels(col_modules_5) <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')

col_modules_6 <-  as.factor(modules_6_th)

levels(col_modules_6) <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')

col_modules_3 <-  as.factor(modules_3_DCA)

levels(col_modules_3) <- c('#fbb4ae','#b3cde3','#ccebc5')

col_modules_2 <-  as.factor(modules_2_DK)

levels(col_modules_2) <- c("tomato","turquoise4" )

col_modules_2dev <-  as.factor(modules_2_dev)

levels(col_modules_2dev) <- c("orange2","darkblue" )

#Plot on surface - adult odontoceti
shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_11, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/11_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/11_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/11_modules2.png")
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_5, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/5_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/5_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/5_modules2.png")
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_6, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/6_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/6_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/6_modules2.png")
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_3, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/3_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/3_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/3_modules2.png")
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_2, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2_modules2.png")
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules_2dev, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2dev_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2dev_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/2dev_modules2.png")
clear3d()

#IMPORT AND PREPARE TREE AND DATA ----
#All
#Import trees in Nexus format - branch lengths needed!!
tree_1 <- "Data/tree_all_early.txt"   #tree with selected families with category data
tree_2 <- "Data/tree_all_late_new.txt"  
tree_3 <- "Data/tree_all_immature.txt"  
tree_4 <- "Data/tree_all_adult.txt"  

##Read the trees for analysis
tree_early <- read.nexus(tree_1) #tree with selected families with category data
plot(tree_early)
tree_late_new <- read.nexus(tree_2) #tree with selected families with category data
plot(tree_late_new)
tree_immature <- read.nexus(tree_3) #tree with selected families with category data
plot(tree_immature)
tree_adult <- read.nexus(tree_4) #tree with selected families with category data
plot(tree_adult)

#Make sure tip labels match taxa names in data frame
wrong_tips_all <- sort(tree_early$tip.label)
tree1 <- tree_early

tree1$tip.label <- levels(as.factor(gdf_mean_early$genus))[match(tree1$tip.label, wrong_tips_all)]
plot(tree1)
tree_early <- tree1

wrong_tips_all <- sort(tree_late_new$tip.label)
tree1 <- tree_late_new

tree1$tip.label <- levels(as.factor(gdf_mean_late_new$genus))[match(tree1$tip.label, wrong_tips_all)]
plot(tree1)
tree_late_new <- tree1

wrong_tips_all <- sort(tree_immature$tip.label)
tree1 <- tree_immature

tree1$tip.label <- levels(as.factor(gdf_mean_immature$genus))[match(tree1$tip.label, wrong_tips_all)]
plot(tree1)
tree_immature <- tree1

wrong_tips_all <- sort(tree_adult$tip.label)
tree1 <- tree_adult

tree1$tip.label <- levels(as.factor(gdf_mean_adult$genus))[match(tree1$tip.label, wrong_tips_all)]
plot(tree1)
tree_adult <- tree1

#List trees
trees_list <- list(tree_early,tree_late_new,tree_immature,tree_adult)

#Mysticeti
#Import trees in Nexus format - branch lengths needed!!
tree_myst_1 <- "Data/tree_all_early_myst.txt"   #tree with selected families with category data
tree_myst_2 <- "Data/tree_all_late_new_myst.txt"  
tree_myst_3 <- "Data/tree_all_immature_myst.txt"  
tree_myst_4 <- "Data/tree_all_adult_myst.txt"  

##Read the trees for analysis
tree_myst_early <- read.nexus(tree_myst_1) #tree with selected families with category data
plot(tree_myst_early)
tree_myst_late_new <- read.nexus(tree_myst_2) #tree with selected families with category data
plot(tree_myst_late_new)
tree_myst_immature <- read.nexus(tree_myst_3) #tree with selected families with category data
plot(tree_myst_immature)
tree_myst_adult <- read.nexus(tree_myst_4) #tree with selected families with category data
plot(tree_myst_adult)

#List trees
trees_myst_list <- list(tree_myst_early,tree_myst_late_new,tree_myst_immature,tree_myst_adult)

#Odontoceti
#Import trees in Nexus format - branch lengths needed!!
tree_odont_1 <- "Data/tree_all_early_odont.txt"   #tree with selected families with category data
tree_odont_2 <- "Data/tree_all_late_new_odont.txt"  
tree_odont_3 <- "Data/tree_all_immature_odont.txt"  
tree_odont_4 <- "Data/tree_all_adult_odont.txt"  

##Read the trees for analysis
tree_odont_early <- read.nexus(tree_odont_1) #tree with selected families with category data
plot(tree_odont_early)
tree_odont_late_new <- read.nexus(tree_odont_2) #tree with selected families with category data
plot(tree_odont_late_new)
tree_odont_immature <- read.nexus(tree_odont_3) #tree with selected families with category data
plot(tree_odont_immature)
tree_odont_adult <- read.nexus(tree_odont_4) #tree with selected families with category data
plot(tree_odont_adult)

#List trees
trees_odont_list <- list(tree_odont_early,tree_odont_late_new,tree_odont_immature,tree_odont_adult)

#COMPARE MODULARITY HYPOTHESES - compare.CR ----

##All data ----
#Perform modularity test all modules - compare selected modules to the null assumption of random assignment of partitions (no modularity at all)
#Run CR analysis for each module hypothesis

#Loop for separate calculation each category
m11_all_phylo <- list()
mR5_all_phylo <- list()
m6_all_phylo <- list()
m3_all_phylo <- list()
m2_all_phylo <- list()
m2dev_all_phylo <- list()

for (c in 1:length(categories_list)){
  m11_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$m11, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  mR5_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$mR5, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  m6_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$m6, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  m3_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$m3, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  m2_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$m2, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  m2dev_all_phylo[[c]] <- phylo.modularity(gdf_mean_shapes[[c]][[1]], modules_df$m2_dev, phy = trees_list[[c]], CI = F, iter = 999, print.progress = T)
  
}

#Compare strength of modularity with CR test
#Loop
CR_compare_all_phylo <- list()

for (c in 1:length(categories_list)){  
  CR_compare_all_phylo[[c]] <- compare.CR(m11_all_phylo[[c]], mR5_all_phylo[[c]], m6_all_phylo[[c]], m3_all_phylo[[c]], m2_all_phylo[[c]],m2dev_all_phylo[[c]], CR.null = TRUE)
}


#The best supported hypothesis is the one with the smallest effect size in summary()
for (c in 1:length(categories_list)){  
  summary(CR_compare_all_phylo[[c]])
}

#Make graph to check which module has the lowest Z score

CR_compare_all_phylo_df <- list()

for (c in 1:length(categories_list)){  
  CR_compare_all_phylo_df[[c]] <- data.frame(modules = c("null","m11_all", "mR5_all", "m6_all","m3_all","m2_all", "m2dev_all"), z = CR_compare_all_phylo[[c]][["sample.z"]], 
                                se = CR_compare_all_phylo[[c]][["sample.se"]])
}

CR_compare_all_phylo_df <- as.data.frame(rbind(CR_compare_all_phylo_df[[1]],CR_compare_all_phylo_df[[2]],CR_compare_all_phylo_df[[3]],CR_compare_all_phylo_df[[4]]))


CR_compare_all_phylo_df <- CR_compare_all_phylo_df[-c(1,8,15,22),] #remove null hypothesis

rownames(CR_compare_all_phylo_df) <- NULL #remove rownames

CR_compare_all_phylo_df$category <- as.factor(rep(levels(categories), each = 6)) #add category variable

#Make factor for variable
CR_compare_all_phylo_df$category <- factor(CR_compare_all_phylo_df$category, 
                                          levels = c("Early Fetus"  ,      "Late Fetus/Neonate", "Juvenile"   ,        "Adult" )) #copy from string printed with the code above
#Order
CR_compare_all_phylo_df <- CR_compare_all_phylo_df[order(CR_compare_all_phylo_df$category),]
CR_compare_all_phylo_df

#Make factor for variable
CR_compare_all_phylo_df$modules <- factor(CR_compare_all_phylo_df$modules, 
                                    levels = c("m11_all", "mR5_all", "m6_all","m3_all","m2_all", "m2dev_all")) #copy from string printed with the code above
#Order
CR_compare_all_phylo_df <- CR_compare_all_phylo_df[order(CR_compare_all_phylo_df$modules),]
CR_compare_all_phylo_df

#Modules list names for plot
module_hyp_list <- c("All sep. (11 mod.)", "5 mod.", "6 mod.","3 mod.","2 mod.", "2 mod. dev.")

CR_compare_all_phylo_plot <- ggplot(CR_compare_all_phylo_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = 24, colour = mypalette_paired[12], fill = adjustcolor(mypalette_paired[12],alpha.f=0.3), stroke = 1.5)+
  ggtitle("Cetacea")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  facet_wrap(vars(category), scales = "free")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_all_phylo_plot

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/4b-Modularity phylo means/Compare_CR_all.txt")
for (c in 1:length(categories_list)){  
  summary(CR_compare_all_phylo[[c]])
}

print("Best modularity model")
for (c in 1:length(categories_list)){  
  summary(m2dev_all_phylo[[c]])
}
sink() 

#Histogram plot best model
for (c in 1:length(categories_list)){  
  plot(m2dev_all_phylo[[c]])
}

##Mysticeti ----
#Run CR analysis for each module hypothesis
#Loop for separate calculation each category
m11_myst_phylo <- list()
mR5_myst_phylo <- list()
m6_myst_phylo <- list()
m3_myst_phylo <- list()
m2_myst_phylo <- list()
m2dev_myst_phylo <- list()

for (c in 1:length(categories_list)){
  m11_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$m11, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  mR5_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$mR5, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  m6_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$m6, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  m3_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$m3, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  m2_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$m2, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  m2dev_myst_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_myst[[c]][[1]], modules_df$m2_dev, phy = trees_myst_list[[c]], CI = F, iter = 999, print.progress = T)
  
}

#Compare strength of modularity with CR test
#Loop
CR_compare_myst_phylo <- list()

for (c in 1:length(categories_list)){  
  CR_compare_myst_phylo[[c]] <- compare.CR(m11_myst_phylo[[c]], mR5_myst_phylo[[c]], m6_myst_phylo[[c]], m3_myst_phylo[[c]], m2_myst_phylo[[c]],m2dev_myst_phylo[[c]], CR.null = TRUE)
}


#The best supported hypothesis is the one with the smallest effect size in summary()
for (c in 1:length(categories_list)){  
  summary(CR_compare_myst_phylo[[c]])
}

#Make graph to check which module has the lowest Z score

CR_compare_myst_phylo_df <- list()

for (c in 1:length(categories_list)){  
  CR_compare_myst_phylo_df[[c]] <- data.frame(modules = c("null","m11_myst", "mR5_myst", "m6_myst","m3_myst","m2_myst", "m2dev_myst"), z = CR_compare_myst_phylo[[c]][["sample.z"]], 
                                             se = CR_compare_myst_phylo[[c]][["sample.se"]])
}

CR_compare_myst_phylo_df <- as.data.frame(rbind(CR_compare_myst_phylo_df[[1]],CR_compare_myst_phylo_df[[2]],CR_compare_myst_phylo_df[[3]],CR_compare_myst_phylo_df[[4]]))


CR_compare_myst_phylo_df <- CR_compare_myst_phylo_df[-c(1,8,15,22),] #remove null hypothesis

rownames(CR_compare_myst_phylo_df) <- NULL #remove rownames

CR_compare_myst_phylo_df$category <- as.factor(rep(levels(categories), each = 6)) #add category variable

#Make factor for variable
CR_compare_myst_phylo_df$category <- factor(CR_compare_myst_phylo_df$category, 
                                           levels = c("Early Fetus"  ,      "Late Fetus/Neonate", "Juvenile"   ,        "Adult" )) #copy from string printed with the code above
#Order
CR_compare_myst_phylo_df <- CR_compare_myst_phylo_df[order(CR_compare_myst_phylo_df$category),]
CR_compare_myst_phylo_df

#Make factor for variable
CR_compare_myst_phylo_df$modules <- factor(CR_compare_myst_phylo_df$modules, 
                                          levels = c("m11_myst", "mR5_myst", "m6_myst","m3_myst","m2_myst", "m2dev_myst")) #copy from string printed with the code above
#Order
CR_compare_myst_phylo_df <- CR_compare_myst_phylo_df[order(CR_compare_myst_phylo_df$modules),]
CR_compare_myst_phylo_df

#Modules list names for plot
module_hyp_list <- c("All sep. (11 mod.)", "5 mod.", "6 mod.","3 mod.","2 mod.", "2 mod. dev.")

CR_compare_myst_phylo_plot <- ggplot(CR_compare_myst_phylo_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = 24, colour = mypalette_groups[1], fill = adjustcolor(mypalette_groups[1],alpha.f=0.3), stroke = 1.5)+
  ggtitle("Mysticeti")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  facet_wrap(vars(category), scales = "free")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_myst_phylo_plot

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/4b-Modularity phylo means/Compare_CR_myst.txt")
for (c in 1:length(categories_list)){  
  summary(CR_compare_myst_phylo[[c]])
}

print("Best modularity model")
for (c in 1:length(categories_list)){  
  summary(m2dev_myst_phylo[[c]])
}
sink() 

#Histogram plot best model
for (c in 1:length(categories_list)){  
  plot(m2dev_myst_phylo[[c]])
}

##Odontoceti ----
#Run CR analysis for each module hypothesis
#Loop for separate calculation each category
m11_odont_phylo <- list()
mR5_odont_phylo <- list()
m6_odont_phylo <- list()
m3_odont_phylo <- list()
m2_odont_phylo <- list()
m2dev_odont_phylo <- list()

for (c in 1:length(categories_list)){
  m11_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$m11, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  mR5_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$mR5, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  m6_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$m6, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  m3_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$m3, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  m2_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$m2, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  m2dev_odont_phylo[[c]] <- phylo.modularity(gdf_mean_shapes_odont[[c]][[1]], modules_df$m2_dev, phy = trees_odont_list[[c]], CI = F, iter = 999, print.progress = T)
  
}

#Compare strength of modularity with CR test
#Loop
CR_compare_odont_phylo <- list()

for (c in 1:length(categories_list)){  
  CR_compare_odont_phylo[[c]] <- compare.CR(m11_odont_phylo[[c]], mR5_odont_phylo[[c]], m6_odont_phylo[[c]], m3_odont_phylo[[c]], m2_odont_phylo[[c]],m2dev_odont_phylo[[c]], CR.null = TRUE)
}


#The best supported hypothesis is the one with the smallest effect size in summary()
for (c in 1:length(categories_list)){  
  summary(CR_compare_odont_phylo[[c]])
}

#Make graph to check which module has the lowest Z score

CR_compare_odont_phylo_df <- list()

for (c in 1:length(categories_list)){  
  CR_compare_odont_phylo_df[[c]] <- data.frame(modules = c("null","m11_odont", "mR5_odont", "m6_odont","m3_odont","m2_odont", "m2dev_odont"), z = CR_compare_odont_phylo[[c]][["sample.z"]], 
                                              se = CR_compare_odont_phylo[[c]][["sample.se"]])
}

CR_compare_odont_phylo_df <- as.data.frame(rbind(CR_compare_odont_phylo_df[[1]],CR_compare_odont_phylo_df[[2]],CR_compare_odont_phylo_df[[3]],CR_compare_odont_phylo_df[[4]]))


CR_compare_odont_phylo_df <- CR_compare_odont_phylo_df[-c(1,8,15,22),] #remove null hypothesis

rownames(CR_compare_odont_phylo_df) <- NULL #remove rownames

CR_compare_odont_phylo_df$category <- as.factor(rep(levels(categories), each = 6)) #add category variable

#Make factor for variable
CR_compare_odont_phylo_df$category <- factor(CR_compare_odont_phylo_df$category, 
                                            levels = c("Early Fetus"  ,      "Late Fetus/Neonate", "Juvenile"   ,        "Adult" )) #copy from string printed with the code above
#Order
CR_compare_odont_phylo_df <- CR_compare_odont_phylo_df[order(CR_compare_odont_phylo_df$category),]
CR_compare_odont_phylo_df

#Make factor for variable
CR_compare_odont_phylo_df$modules <- factor(CR_compare_odont_phylo_df$modules, 
                                           levels = c("m11_odont", "mR5_odont", "m6_odont","m3_odont","m2_odont", "m2dev_odont")) #copy from string printed with the code above
#Order
CR_compare_odont_phylo_df <- CR_compare_odont_phylo_df[order(CR_compare_odont_phylo_df$modules),]
CR_compare_odont_phylo_df

#Modules list names for plot
module_hyp_list <- c("All sep. (11 mod.)", "5 mod.", "6 mod.","3 mod.","2 mod.", "2 mod. dev.")

CR_compare_odont_phylo_plot <- ggplot(CR_compare_odont_phylo_df, aes(x=modules, y=z)) + 
  geom_errorbar(aes(ymin=z-se, ymax=z+se), width = 0.5, colour = "gray30") +
  geom_point(size = 6, shape = 24, colour = mypalette_groups[2], fill = adjustcolor(mypalette_groups[2],alpha.f=0.3), stroke = 1.5)+
  ggtitle("Odontoceti")+
  xlab("Modularity hypothesis")+
  ylab("Z-scores")+
  facet_wrap(vars(category), scales = "free")+
  scale_x_discrete(labels = module_hyp_list)+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
CR_compare_odont_phylo_plot

#Save compare CR and best modularity model results to file - include pairwise modularity scores
sink("Output/4b-Modularity phylo means/Compare_CR_odont.txt")
for (c in 1:length(categories_list)){  
  summary(CR_compare_odont_phylo[[c]])
}

print("Best modularity model")
for (c in 1:length(categories_list)){  
  summary(m2dev_odont_phylo[[c]])
}
sink() 

#Histogram plot best model
for (c in 1:length(categories_list)){  
  plot(m2dev_odont_phylo[[c]])
}

###Scatter plot and bar plot to compare Z-scores modularity hypothesis by group ----

#Create combined data frame
CR_compare_all_phylo_df$data <- "All"
CR_compare_myst_phylo_df$data <- levels(groups)[1]
CR_compare_odont_phylo_df$data <- levels(groups)[2]

CR_compare_plot_phylo_df <- bind_rows(CR_compare_all_phylo_df, CR_compare_myst_phylo_df, CR_compare_odont_phylo_df)

#Add column with proper mod. hypothesis labels
CR_compare_plot_phylo_df$labels <- rep(module_hyp_list, each = 4, times = 3)
CR_compare_plot_phylo_df

#Calculate error for each group and module hyp.
CR_compare_plot_phylo_df <- CR_compare_plot_phylo_df %>% 
  group_by(data, modules) %>% mutate(z_min = min(z), z_max = max(z))
CR_compare_plot_phylo_df

#Make vector to fix labels for plot
sub_labels_string <- paste0("___", levels(as.factor(CR_compare_plot_phylo_df$data)))

#Create a vertical bar plot faceted by group
CR_compare_phylo_bar_plot <- ggplot(CR_compare_plot_phylo_df, aes(x = reorder_within(labels, z, data), y = z, fill = data)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8, color = "gray50") +
  geom_errorbar(aes(ymin=z_min, ymax=z_max), width=.2,
                position=position_dodge(.9)) +
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
CR_compare_phylo_bar_plot

###Heatmaps plots signifcant difference in modularity ----
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
modularity_z_phylo_all <- list()

for (c in 1:length(categories_list)){  
modularity_z_phylo_all[[c]] <- CR_compare_all_phylo[[c]][["pairwise.z"]]
}

modularity_pvals_phylo_all <- list()

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_all[[c]] <- CR_compare_all_phylo[[c]][["pairwise.P"]]
}


#Remove null hypothesis
for (c in 1:length(categories_list)){  
  modularity_z_phylo_all[[c]] <- modularity_z_phylo_all[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_z_phylo_all[[c]] <- modularity_z_phylo_all[[c]][,-1]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_all[[c]] <- modularity_pvals_phylo_all[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_all[[c]] <- modularity_pvals_phylo_all[[c]][,-1]
}

#Set correct row and col names for both
for (c in 1:length(categories_list)){  
  rownames(modularity_z_phylo_all[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_z_phylo_all[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  rownames(modularity_pvals_phylo_all[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_pvals_phylo_all[[c]]) <- module_hyp_list
}

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_phylo_all_upper_tri <- list()
  
for (c in 1:length(categories_list)){    
  modularity_z_phylo_all_upper_tri[[c]] <- get_upper_tri(modularity_z_phylo_all[[c]])
}

modularity_pvals_phylo_all_upper_tri <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_all_upper_tri[[c]] <- get_upper_tri(modularity_pvals_phylo_all[[c]])
}

#Melt to make table in the format needed for heatmap
modularity_z_phylo_all_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_z_phylo_all_melt[[c]] <-  melt(modularity_z_phylo_all_upper_tri[[c]], value.name = "z", na.rm = TRUE)
}

modularity_pvals_phylo_all_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_all_melt[[c]] <-  melt(modularity_pvals_phylo_all_upper_tri[[c]], value.name = "p", na.rm = TRUE)
}


#Add column to main table

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_all_melt[[c]]$z <- modularity_z_phylo_all_melt[[c]]$z
}

#Create columns where only significant values are shown
for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_all_melt[[c]] <- modularity_pvals_phylo_all_melt[[c]] %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                  z_if_sig = ifelse(sig_p, z, NA)) %>%
  mutate_at(vars(starts_with("z")), list(~ round(., 2)))
}

modularity_pvals_phylo_all_melt 

#Check if all NAs in some categories
#Check NAs first - if TRUE do not plot

for (c in 1:length(categories_list)){    
  print(all(is.na(modularity_pvals_phylo_all_melt[[c]]$p_if_sig)))
}

#Nice heatmap plot for each category
modularity_pvals_phylo_all_early_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_all_melt[[1]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_all_melt[[1]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[1])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
modularity_pvals_phylo_all_early_heatmap_ggplot 

modularity_pvals_phylo_all_immature_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_all_melt[[3]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_all_melt[[3]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[3])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_all_immature_heatmap_ggplot 

modularity_pvals_phylo_all_adult_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_all_melt[[4]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_all_melt[[4]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[4])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_all_adult_heatmap_ggplot 

plotMP <- ggarrange(modularity_pvals_phylo_all_early_heatmap_ggplot, modularity_pvals_phylo_all_immature_heatmap_ggplot, modularity_pvals_phylo_all_adult_heatmap_ggplot,
                   nrow = 1, ncol = 3, common.legend = F)
plotMP <- annotate_figure(plotMP, top = text_grob("All data", face = "bold", size = 17))
plotMP

#Mysticeti
#Save p-values as object
modularity_z_phylo_myst <- list()

for (c in 1:length(categories_list)){  
  modularity_z_phylo_myst[[c]] <- CR_compare_myst_phylo[[c]][["pairwise.z"]]
}

modularity_pvals_phylo_myst <- list()

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_myst[[c]] <- CR_compare_myst_phylo[[c]][["pairwise.P"]]
}


#Remove null hypothesis
for (c in 1:length(categories_list)){  
  modularity_z_phylo_myst[[c]] <- modularity_z_phylo_myst[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_z_phylo_myst[[c]] <- modularity_z_phylo_myst[[c]][,-1]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_myst[[c]] <- modularity_pvals_phylo_myst[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_myst[[c]] <- modularity_pvals_phylo_myst[[c]][,-1]
}

#Set correct row and col names for both
for (c in 1:length(categories_list)){  
  rownames(modularity_z_phylo_myst[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_z_phylo_myst[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  rownames(modularity_pvals_phylo_myst[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_pvals_phylo_myst[[c]]) <- module_hyp_list
}

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_phylo_myst_upper_tri <- list()

for (c in 1:length(categories_list)){    
  modularity_z_phylo_myst_upper_tri[[c]] <- get_upper_tri(modularity_z_phylo_myst[[c]])
}

modularity_pvals_phylo_myst_upper_tri <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_myst_upper_tri[[c]] <- get_upper_tri(modularity_pvals_phylo_myst[[c]])
}

#Melt to make table in the format needed for heatmap
modularity_z_phylo_myst_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_z_phylo_myst_melt[[c]] <-  melt(modularity_z_phylo_myst_upper_tri[[c]], value.name = "z", na.rm = TRUE)
}

modularity_pvals_phylo_myst_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_myst_melt[[c]] <-  melt(modularity_pvals_phylo_myst_upper_tri[[c]], value.name = "p", na.rm = TRUE)
}


#Add column to main table

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_myst_melt[[c]]$z <- modularity_z_phylo_myst_melt[[c]]$z
}

#Create columns where only significant values are shown
for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_myst_melt[[c]] <- modularity_pvals_phylo_myst_melt[[c]] %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                                          z_if_sig = ifelse(sig_p, z, NA)) %>%
    mutate_at(vars(starts_with("z")), list(~ round(., 2)))
}

modularity_pvals_phylo_myst_melt 

#Check if all NAs in some categories
#Check NAs first - if TRUE do not plot

for (c in 1:length(categories_list)){    
  print(all(is.na(modularity_pvals_phylo_myst_melt[[c]]$p_if_sig)))
}

#Nice heatmap plot for each category
modularity_pvals_phylo_myst_early_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_myst_melt[[1]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_myst_melt[[1]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[1])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
modularity_pvals_phylo_myst_early_heatmap_ggplot 

modularity_pvals_phylo_myst_late_new_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_myst_melt[[2]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_myst_melt[[2]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[2])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_myst_late_new_heatmap_ggplot 

modularity_pvals_phylo_myst_immature_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_myst_melt[[3]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_myst_melt[[3]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[3])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_myst_immature_heatmap_ggplot 

modularity_pvals_phylo_myst_adult_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_myst_melt[[4]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_myst_melt[[4]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[4])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_myst_adult_heatmap_ggplot 

plotMP1 <- ggarrange(modularity_pvals_phylo_myst_early_heatmap_ggplot, modularity_pvals_phylo_myst_late_new_heatmap_ggplot, modularity_pvals_phylo_myst_immature_heatmap_ggplot, modularity_pvals_phylo_myst_adult_heatmap_ggplot,
                    nrow = 2, ncol = 2, common.legend = F)
plotMP1 <- annotate_figure(plotMP1, top = text_grob("Mysticeti", face = "bold", size = 17))
plotMP1

#Odontoceti
#Save p-values as object
modularity_z_phylo_odont <- list()

for (c in 1:length(categories_list)){  
  modularity_z_phylo_odont[[c]] <- CR_compare_odont_phylo[[c]][["pairwise.z"]]
}

modularity_pvals_phylo_odont <- list()

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_odont[[c]] <- CR_compare_odont_phylo[[c]][["pairwise.P"]]
}


#Remove null hypothesis
for (c in 1:length(categories_list)){  
  modularity_z_phylo_odont[[c]] <- modularity_z_phylo_odont[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_z_phylo_odont[[c]] <- modularity_z_phylo_odont[[c]][,-1]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_odont[[c]] <- modularity_pvals_phylo_odont[[c]][-1,]
}

for (c in 1:length(categories_list)){  
  modularity_pvals_phylo_odont[[c]] <- modularity_pvals_phylo_odont[[c]][,-1]
}

#Set correct row and col names for both
for (c in 1:length(categories_list)){  
  rownames(modularity_z_phylo_odont[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_z_phylo_odont[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  rownames(modularity_pvals_phylo_odont[[c]]) <- module_hyp_list
}

for (c in 1:length(categories_list)){  
  colnames(modularity_pvals_phylo_odont[[c]]) <- module_hyp_list
}

#Get upper triangles only - half matrix, eliminates redundant info
modularity_z_phylo_odont_upper_tri <- list()

for (c in 1:length(categories_list)){    
  modularity_z_phylo_odont_upper_tri[[c]] <- get_upper_tri(modularity_z_phylo_odont[[c]])
}

modularity_pvals_phylo_odont_upper_tri <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_odont_upper_tri[[c]] <- get_upper_tri(modularity_pvals_phylo_odont[[c]])
}

#Melt to make table in the format needed for heatmap
modularity_z_phylo_odont_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_z_phylo_odont_melt[[c]] <-  melt(modularity_z_phylo_odont_upper_tri[[c]], value.name = "z", na.rm = TRUE)
}

modularity_pvals_phylo_odont_melt <- list()

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_odont_melt[[c]] <-  melt(modularity_pvals_phylo_odont_upper_tri[[c]], value.name = "p", na.rm = TRUE)
}


#Add column to main table

for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_odont_melt[[c]]$z <- modularity_z_phylo_odont_melt[[c]]$z
}

#Create columns where only significant values are shown
for (c in 1:length(categories_list)){    
  modularity_pvals_phylo_odont_melt[[c]] <- modularity_pvals_phylo_odont_melt[[c]] %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                            p_if_sig = ifelse(sig_p, p, NA),
                                                                                            z_if_sig = ifelse(sig_p, z, NA)) %>%
    mutate_at(vars(starts_with("z")), list(~ round(., 2)))
}

modularity_pvals_phylo_odont_melt 

#Check if all NAs in some categories
#Check NAs first - if TRUE do not plot

for (c in 1:length(categories_list)){    
  print(all(is.na(modularity_pvals_phylo_odont_melt[[c]]$p_if_sig)))
}

#Nice heatmap plot for each category
modularity_pvals_phylo_odont_early_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_odont_melt[[1]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_odont_melt[[1]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[1])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
modularity_pvals_phylo_odont_early_heatmap_ggplot 

modularity_pvals_phylo_odont_late_new_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_odont_melt[[2]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_odont_melt[[2]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[2])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_odont_late_new_heatmap_ggplot 

modularity_pvals_phylo_odont_immature_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_odont_melt[[3]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_odont_melt[[3]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[3])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_odont_immature_heatmap_ggplot 

modularity_pvals_phylo_odont_adult_heatmap_ggplot <- ggplot(data = modularity_pvals_phylo_odont_melt[[4]], aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_pvals_phylo_odont_melt[[4]]$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle (levels(categories)[4])+ 
  theme(plot.title = element_text(face = 3, hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank())+
  guides(fill = "none")
modularity_pvals_phylo_odont_adult_heatmap_ggplot 

plotMP2 <- ggarrange(modularity_pvals_phylo_odont_early_heatmap_ggplot, modularity_pvals_phylo_odont_late_new_heatmap_ggplot, modularity_pvals_phylo_odont_immature_heatmap_ggplot, modularity_pvals_phylo_odont_adult_heatmap_ggplot,
                     nrow = 2, ncol = 2, common.legend = F)
plotMP2 <- annotate_figure(plotMP2, top = text_grob("Odontoceti", face = "bold", size = 17))
plotMP2

##Plot modules best hypotheses----

#All data+Odont
#Plot on surface
shade3d(refmesh_phylo_all, col = "white", alpha = 0.5)
spheres3d(shape_array[,,41], col = col_modules_2dev, type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/all_phylo_odont_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/all_phylo_odont_modules1.png")
rgl.snapshot(filename = "Output/4b-Modularity phylo means/all_phylo_odont_modules2.png")
play3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6, movie = "all_phylo_odont_modules" ,dir = "Output/4b-Modularity phylo means/")
clear3d()

#Myst
shade3d(myst_fetus, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Ff3",Ids)], col =  col_modules_2dev, type = "s",
          radius = 1, aspect = T, main = "mean",axes = F, main = F, fov = 0)
rgl.snapshot(filename = "Output/4b-Modularity phylo means/myst_modules.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/myst_modules1.png") 
rgl.snapshot(filename = "Output/4b-Modularity phylo means/myst_modules2.png") 
play3d(spin3d(axis = c(0, 1,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 1,1), rpm = 10), duration = 6, movie = "myst_modules" ,dir = "Output/4b-Modularity phylo means/")

##Modularity by growth stage and group ----
#Mysticeti
#Compare CR
CR_compare_myst_cat_phylo <- compare.CR(m2dev_myst_phylo[[1]], m2dev_myst_phylo[[2]],m2dev_myst_phylo[[3]],m2dev_myst_phylo[[4]], CR.null = TRUE)

summary(CR_compare_myst_cat_phylo)

#Odontoceti
#Compare CR
CR_compare_odont_cat_phylo <- compare.CR(m2dev_odont_phylo[[1]], m2dev_odont_phylo[[2]],m2dev_odont_phylo[[3]],m2dev_odont_phylo[[4]], CR.null = TRUE)

summary(CR_compare_odont_cat_phylo)

#Compare CR all
CR_compare_myst_odont_cat_phylo <- compare.CR(m2dev_myst_phylo[[1]], m2dev_myst_phylo[[2]],m2dev_myst_phylo[[3]],m2dev_myst_phylo[[4]],
                                              m2dev_odont_phylo[[1]], m2dev_odont_phylo[[2]],m2dev_odont_phylo[[3]],m2dev_odont_phylo[[4]],  CR.null = TRUE)


summary(CR_compare_myst_odont_cat_phylo)

#Save compare CR results to file - include pairwise modularity scores
sink("Output/4b-Modularity phylo means/Compare_CR_myst_odont_cat_phylo.txt")
print("Mysticeti")
summary(CR_compare_myst_cat_phylo)
print("Odontoceti")
summary(CR_compare_odont_cat_phylo)
print("All")
summary(CR_compare_myst_odont_cat_phylo)
sink() 

###Line plots for CR of groups at each stage ----
#Make dataset with all CRs
CR_group_cat_phylo1 <- list()
CR_group_cat_phylo2 <- list()

for (c in 1:length(categories_list)){    
  CR_group_cat_phylo1[[c]] <- m2dev_myst_phylo[[c]][["CR"]]
  CR_group_cat_phylo2[[c]] <- m2dev_odont_phylo[[c]][["CR"]]
}

CR_group_cat_phylo <- data.frame(CR = rbind(CR_group_cat_phylo1[[1]],CR_group_cat_phylo1[[2]],CR_group_cat_phylo1[[3]],CR_group_cat_phylo1[[4]],
                                          CR_group_cat_phylo2[[1]],CR_group_cat_phylo2[[2]],CR_group_cat_phylo2[[3]],CR_group_cat_phylo2[[4]]))
  
#Add labels and other attributes to tibble as columns
CR_group_cat_phylo <- CR_group_cat_phylo %>% 
  mutate(category= rep(categories_list, times = 2), group = rep(groups_list, each = 4))
CR_group_cat_phylo 

#Nice line plot by group
CR_group_cat_phylo_ggplot <- ggplot(CR_group_cat_phylo, aes(x = category, y = CR)) + 
  geom_line(aes(x = category, y = CR, linetype = group, colour = group, group = group), linewidth = 1, inherit.aes = F)+
  geom_point(aes(color = group,
                 shape = group),size = 4, fill = "white", stroke = 1.5)+
  scale_shape_manual(name = "Groups", labels = levels(groups), values = shapes)+
  scale_colour_manual(name = "Groups", labels =levels(groups),
                      values = mypalette_groups, aesthetics = c("colour","fill"))+   
  scale_linetype_manual(name = "Groups", labels = levels(groups),
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("CR")+
  xlab("Growth stage")+
  ggtitle ("Modularity index (CR) group and growth stage")+ 
  scale_x_discrete(labels = levels(as.factor(categories_list_short)))+
  ylim(0.7, 1.05)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),legend.key = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal", strip.text.x = element_text(size=12))+
  guides(colour = "none", linetype = "none", shape = "none")
CR_group_cat_phylo_ggplot

#Annotate plot
CR_group_cat_phylo_ggplot <-
  CR_group_cat_phylo_ggplot +
  annotate("rect",  xmin = -Inf, xmax = Inf, ymin = 1.0, ymax = 1.05,
           alpha = .5,fill = "gray70")+
  annotate("text", x = 4, y = c(1.02, 0.98), label = c("Fully integrated", "Modular"), size = 5)
CR_group_cat_phylo_ggplot

#Add phylopics
CR_group_cat_phylo_ggplot  <- 
  CR_group_cat_phylo_ggplot + # 1 line per taxon, alphabetical order
  add_phylopic(myst, alpha = 1, x = 2.3, y = 0.95, ysize = 0.02, fill = mypalette_groups[1])+
  add_phylopic(odont, alpha = 1, x = 3.2, y = 0.75, ysize = 0.018, fill = mypalette_groups[2])
CR_group_cat_phylo_ggplot

###Heatmap plot significant difference modularity group and category ----
#All data
#Save p-values as object
modularity_grp_cat_z_phylo_all <- CR_compare_myst_odont_cat_phylo[["pairwise.z"]]
modularity_grp_cat_pvals_phylo_all <- CR_compare_myst_odont_cat_phylo[["pairwise.P"]]

#Remove null hypothesis
modularity_grp_cat_z_phylo_all <- modularity_grp_cat_z_phylo_all[-1,]
modularity_grp_cat_z_phylo_all <- modularity_grp_cat_z_phylo_all[,-1]
modularity_grp_cat_pvals_phylo_all <- modularity_grp_cat_pvals_phylo_all[-1,]
modularity_grp_cat_pvals_phylo_all <- modularity_grp_cat_pvals_phylo_all[,-1]  

#Set correct row and col names for both

#Make lists for labels
groups_list <- str_to_lower(levels(groups))
groups_list_short <- c("Myst", "Odont")
categories_list <- levels(gdf$category)
categories_list_short <- c("1","2","3","4")

modularity_grp_cat_names <- paste0(rep(groups_list_short, each =4), "_", rep(categories_list_short, times  =2))

rownames(modularity_grp_cat_z_phylo_all) <- modularity_grp_cat_names
rownames(modularity_grp_cat_pvals_phylo_all) <- modularity_grp_cat_names
colnames(modularity_grp_cat_z_phylo_all) <- modularity_grp_cat_names
colnames(modularity_grp_cat_pvals_phylo_all) <- modularity_grp_cat_names

#Get upper triangles only - half matrix, eliminates redundant info
modularity_grp_cat_z_phylo_all_upper_tri <- get_upper_tri(modularity_grp_cat_z_phylo_all)
modularity_grp_cat_pvals_phylo_all_upper_tri <- get_upper_tri(modularity_grp_cat_pvals_phylo_all)

#Melt to make table in the format needed for heatmap
modularity_grp_cat_z_phylo_all_melt <- melt(modularity_grp_cat_z_phylo_all_upper_tri, value.name = "z", na.rm = TRUE)
modularity_grp_cat_pvals_phylo_all_melt <- melt(modularity_grp_cat_pvals_phylo_all_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
modularity_grp_cat_pvals_phylo_all_melt$z <- modularity_grp_cat_z_phylo_all_melt$z

#Create columns where only significant values are shown
modularity_grp_cat_pvals_phylo_all_melt <- modularity_grp_cat_pvals_phylo_all_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                                  p_if_sig = ifelse(sig_p, p, NA),
                                                                                  z_if_sig = ifelse(sig_p, z, NA)) %>%
  mutate_at(vars(starts_with("z")), list(~ round(., 2)))
modularity_grp_cat_pvals_phylo_all_melt 

#Nice heatmap plot
modularity_grp_cat_pvals_phylo_all_heatmap_ggplot <- ggplot(data = modularity_grp_cat_pvals_phylo_all_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = z_if_sig), color = "white", size = 4.5) +
  scale_fill_gradient2(low = mypalette_seq_modularity[9], high = mypalette_seq_modularity[1], mid = mypalette_seq_modularity[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(min(modularity_grp_cat_pvals_phylo_all_melt$p_if_sig), 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_modularity[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Modularity z-score differences by group and category")+ 
  theme(plot.title = element_text(face = 2, hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = "inside", legend.position.inside = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
modularity_grp_cat_pvals_phylo_all_heatmap_ggplot 

###### 

#Next - ch. 5b - phyloPCA means
