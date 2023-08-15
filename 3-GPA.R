
#===========================================================#
#                                                           #
#     SKULL MODULARITY - MYSTICETI & ODONTOCETI             #
#                                                           #
#===========================================================#


#CH.3 - Prepare final dataset for analysis, run GPA whole skull, rostrum  and braincase

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
library(car)
library(Rvcg)
library(scales)

#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#DATA PREP ----

###SET WD to root folder from console!! -->

#Import classifiers
classifiers <- read_csv("Data/specimens_all.csv")
glimpse(classifiers)

#Make sure the specimens are in the same order before proceeding!!!
identical(dimnames(final_dataset)[[3]], classifiers$specimen, attrib.as.set = T)

#Create new columns needed for analyses
#Calculate max TL and BZW for each genus
#Use genus2 for more even grouping
max_TL_BZW <- classifiers %>% group_by(genus2) %>% summarize(max_TL = max(TL_m), max_BZW = max(BZW_mm))

#Count number of specimens per genus
count_genus <- classifiers %>% count(genus2)

#Create new columns with max TL and BZW repeated for all specimens of each genus

max_TL <- c()

for (i in 1:length(count_genus$n)) {
  x<-rep(max_TL_BZW$max_TL[i],count_genus$n[i])
  max_TL <- c(max_TL,x)
}

max_BZW <- c()

for (i in 1:length(count_genus$n)) {
  x<-rep(max_TL_BZW$max_BZW[i],count_genus$n[i])
  max_BZW <- c(max_BZW,x)
}

classifiers$max_TL <- max_TL
classifiers$max_BZW <- max_BZW

#Calculate % max TL and max BZW for each specimen in each genus
classifiers <- classifiers %>% mutate(TL_100 = (TL_m*100/max_TL), BZW_100 = (BZW_mm*100/max_BZW))   
glimpse(classifiers)

#Count how many stages and category per genus2 - useful to decide which binning to use in analyses
count_category <-  classifiers %>% group_by(genus2) %>% count(category)
count_stage <-  classifiers %>% group_by(genus2) %>% count(stage)

View(count_category) #category more even
View(count_stage)

genera_list <- levels(as.factor(classifiers$genus2))
families_list <- levels(as.factor(classifiers$family))

count_category <- count_category %>% 
  mutate(group = if_else(genus2 %in% genera_list[c(1:4,7,13)], "mysticeti", "odontoceti"),
         family = if_else(genus2 %in% genera_list[c(1)],families_list[1], if_else(genus2 %in% genera_list[c(4)],families_list[8], 
                                                                                  if_else(genus2 %in% genera_list[c(5,15)],families_list[7],  if_else(genus2 %in% genera_list[c(2,3,7,13)],families_list[2],
                                                                                                                                                      if_else(genus2 %in% genera_list[c(16,20)],families_list[10], if_else(genus2 %in% genera_list[c(23)],families_list[6], 
                                                                                                                                                                                                                           if_else(genus2 %in% genera_list[c(22)],families_list[12], if_else(genus2 %in% genera_list[c(10,14)],families_list[13],
                                                                                                                                                                                                                                                                                             if_else(genus2 %in% genera_list[c(11,21)],families_list[11],families_list[3]))))))))))

View(count_category)

#Check which genera have all multiple specimens for each category - useful to choose for later analyses

count_category_n2 <- count_category %>% filter(n >= 2)
View(count_category_n2)

#Add prenatal/postantal column to check which genera are ok for simplified traj analyses (at least 2 per group)
count_category_n2 <- count_category_n2 %>% mutate(early = if_else(category == "1-early", 1, 0),
                                                  late_new = if_else(category == "2-late/new", 1, 0), immature = if_else(category == "3-immature", 1, 0),
                                                  adult = if_else(category == "4-adult", 1, 0), prenatal = if_else(category %in% c("1-early", "2-late/new"), 1, 0),
                                                  postnatal = if_else(category %in% c("3-immature", "4-adult"), 1, 0))

View(count_category_n2)

#Check if groups have all categories or stages
count_category <- count_category %>% mutate(early = if_else(category == "1-early", 1, 0),
                                            late_new = if_else(category == "2-late/new", 1, 0), immature = if_else(category == "3-immature", 1, 0),
                                            adult = if_else(category == "4-adult", 1, 0), prenatal = if_else(category %in% c("1-early", "2-late/new"), 1, 0),
                                            postnatal = if_else(category %in% c("3-immature", "4-adult"), 1, 0))
count_category <- count_category %>% mutate(stage = ifelse(category %in% c("3-immature", "4-adult"),"postantal", "prenatal"))

prenatal_groups <- count_category %>% group_by(group) %>% filter(prenatal == 1) %>% summarise(sum(n))
postnatal_groups <- count_category  %>% group_by(group) %>% filter(postnatal == 1) %>% summarise(sum(n))

early_groups <- count_category %>% group_by(group) %>% filter(early == 1) %>% summarise(sum(n))
late_groups <- count_category %>% group_by(group) %>% filter(late_new == 1) %>% summarise(sum(n))
immature_groups <- count_category %>% group_by(group) %>% filter(immature == 1) %>% summarise(sum(n))
adult_groups <- count_category %>% group_by(group) %>% filter(adult == 1) %>% summarise(sum(n))

prenatal_postnatal_groups <- Reduce(intersect, list(prenatal_groups$group,postnatal_groups$group))
all_categories_groups <- Reduce(intersect, list(early_groups$group,late_groups$group,
                                                immature_groups$group, adult_groups$group))

count_group_category <- count_category %>%  group_by(group) %>% 
  count(early, late_new, immature, adult) %>% mutate(category = c("adult", "immature", "late_new", "early"))
count_group_stages <- count_category %>%  group_by(group) %>% 
  count(prenatal, postnatal) %>% mutate(category = c("postnatal", "prenatal"))

prenatal_families <- count_category %>% group_by(family) %>% filter(prenatal == 1) %>% summarise(sum(n))
postnatal_families <- count_category  %>% group_by(family) %>% filter(postnatal == 1) %>% summarise(sum(n))

early_families <- count_category %>% group_by(family) %>% filter(early == 1) %>% summarise(sum(n))
late_families <- count_category %>% group_by(family) %>% filter(late_new == 1) %>% summarise(sum(n))
immature_families <- count_category %>% group_by(family) %>% filter(immature == 1) %>% summarise(sum(n))
adult_families <- count_category %>% group_by(family) %>% filter(adult == 1) %>% summarise(sum(n))

prenatal_postnatal_families <- Reduce(intersect, list(prenatal_families$family,postnatal_families$family))
all_categories_families <- Reduce(intersect, list(early_families$family,late_families$family,
                                                  immature_families$family, adult_families$family))

count_family_category <- count_category %>%  group_by(family, category) %>% summarise(sum(n)) %>% 
  filter(n() >= 4) 

count_family_stages <- count_category %>%  group_by(family, stage) %>% summarise(sum(n))%>% 
  filter(n() >= 2) 

##Order shape data by category, useful for plot legend
#Check levels/category names
as.factor(classifiers$category)

#Order shape data as category
order_dataframe <- geomorph.data.frame(raw_data = final_dataset, category = classifiers$category)
#Check specimens order
dimnames(order_dataframe$raw_data)[[3]]

#Order dataframe
order_dataframe$category <- factor(order_dataframe$category,
                                   levels = c("1-early", "2-late/new", "3-immature", "4-adult"))

#Create new shape data object ordered by category
shape_array <- order_dataframe$raw_data[,,order(order_dataframe$category)]
#Check it worked
dimnames(shape_array)[[3]]

##Order classifiers data by category, useful for plot legend
#Make factor for variable
classifiers$category <- factor(classifiers$category, 
                               levels = c("1-early", "2-late/new", "3-immature", "4-adult")) #copy from string printed with the code above
#Order
classifiers <- classifiers[order(classifiers$category),]
View(classifiers)

#Check specimens and classifiers are in the same order
identical(dimnames(shape_array)[[3]], classifiers$specimen, attrib.as.set = T) 

##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(shape_array)

#Import simplified ply
refmesh_all <- vcgImport("Data/refmesh_all.ply")

#Define fixed LMs and shape array only with LMs
fixed_LMs <- c(1:64)

shape_array_LM <- shape_array[fixed_LMs,,]

#Plot on surface
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[fixed_LMs,,41], col =  "firebrick", type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(shape_array[-fixed_LMs,,41], col =  "tomato", type = "s",
          radius = 0.4, aspect = T, main = "mean",axes = F, main = F, fov = 0)
text3d(shape_array_LM[,1,41], shape_array_LM[,2,41], shape_array_LM[,3,41], 
       texts = fixed_LMs, pos = 4, offset = 1, font = 2) #change pos

rgl.snapshot(filename = "Output/landmarks_dorsal.png") 
rgl.snapshot(filename = "Output/landmarks_lateral1.png") 
rgl.snapshot(filename = "Output/landmarks_lateral2.png") 
rgl.snapshot(filename = "Output/landmarks_ventral.png") 
rgl.snapshot(filename = "Output/landmarks_posterior.png") 
rgl.snapshot(filename = "Output/landmarks_anterior.png") 

play3d(spin3d(axis = c(0, 1,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6, movie = "landmarks" ,dir = "Output/")

#Check for outliers in raw data shape array, they would be displayed in red
#Might be due to absent bones, check misstable list
plotOutliers(shape_array)
#Plot outliers landmarks to search for possible problems - check file name to find number in classifiers
checkLM(shape_array, path="Data/ply/", pt.size = 5, suffix=".ply", render = "s", begin = 162, point = "s")

#Save txt file of output
sink("Output/outliers_raw.txt")
print(plotOutliers(shape_array))
sink() 
#Mark where outliers end based on plot

#Save specimens names as object
specimens <- dimnames(shape_array)[[3]]

#Save Id as object, useful for later analysis
Ids <- classifiers$code

#Save growth stages as factor, useful for later analysis
categories <- as.factor(classifiers$category)
#Check how many colors are needed
levels(categories) #4
#To capitalize or change spelling
levels(categories) <- c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult")

#Save families as factor, useful for later analysis
families <- as.factor(classifiers$family)
#Check how many colors are needed
length(levels(families)) #13
#To capitalize or change spelling
levels(families) <- c("Balaenidae" ,     "Balaenopteridae" , "Delphinidae"  ,   "Delphininae" ,    "Globicephalinae", "Inioidea",    
                      "Monodontidae" ,   "Neobalaenidae" , "Orcininae"     ,  "Phocoenidae",     "Physeteroidea" ,  "Platanistidae" , 
                      "Ziphiidae")

#Save groups as factor, useful for later analysis
groups <- as.factor(classifiers$group)
#Check how many colors are needed
groups #2
#To capitalize or change spelling
levels(groups) <- c("Mysticeti","Odontoceti")


##Create project palettes----
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_hueyellowblue <- brewer.pal(9,"YlGnBu")
image(1:9, 1, as.matrix(1:9), col = mypalette_hueyellowblue, xlab = "Y-B Hues",
      ylab = "", yaxt = "n")

mypalette_huepurple <- brewer.pal(9,"RdPu")
image(1:9, 1, as.matrix(1:9), col = mypalette_huepurple, xlab = "Purple Hues",
      ylab = "", yaxt = "n")

mypalette_huegreen <- brewer.pal(9,"YlGn")
image(1:9, 1, as.matrix(1:9), col = mypalette_huegreen, xlab = "Green Hues",
      ylab = "", yaxt = "n")

mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", yaxt = "n")

#Palette for 27 genera - based on genus2, different shades for myst and odont
mypalette_myst <- colorRampPalette(c(mypalette_huepurple[3:6], mypalette_huepurple[8]))
mypalette_myst(6)
plot(rep(1,6),col=mypalette_myst(6),pch=19,cex=3)
mypalette_myst2 <- mypalette_myst(6)

mypalette_odont <- colorRampPalette(c(mypalette_hueyellowblue[2:4],mypalette_huegreen[4:6],
                                      mypalette_hueyellowblue[5:6], mypalette_huegreen[7:8],
                                      mypalette_hueyellowblue[8])) #mix different hues to help
mypalette_odont(21)
plot(rep(1,21),col=mypalette_odont(21),pch=19,cex=3)
mypalette_odont2 <- mypalette_odont(21)

#Palette for families - toextarct gorup colors
levels(families)

mypalette_families <- c(mypalette_myst2[4], mypalette_myst2[3], mypalette_odont2[6], mypalette_odont2[14], #"balaenidae"      "balaenopteridae" "delphinidae"     "delphininae" 
                        mypalette_odont2[10], mypalette_odont2[17], mypalette_odont2[20], #"globicephalinae" "inioidea" "monodontidae"    
                        mypalette_myst2[2], mypalette_odont2[3], mypalette_odont2[13], mypalette_odont2[2], #"neobalaenidae" "orcininae"       "phocoenidae"     "physeteroidea"
                        mypalette_odont2[8], mypalette_odont2[16]) #"platanistidae" "ziphiidae"

plot(rep(1,13),col=mypalette_families ,pch=19,cex=3, main = "Families colors", ylab = "", xlab = "" ,cex.main = 2,
     yaxt = "n", xaxt = "n")
title(xlab = "1-Balaenidae 2-Balaenopteridae 3-Delphinidae 4-Delphininae     
5-Globicephalinae 6-Inioidea 7-Monodontidae 8-Neobalaenidae 
9-Orcininae 10-Phocoenidae 11-Physeteroidea 12-Platanistidae 13-Ziphiidae", cex.lab = 1.3, font.lab = 3, line = -3)
text(x = seq_along(1:13), y = 1.05, labels = seq_along(1:13))

#Palette for categories - early, late/new, immature, adult
mypalette_category <- c(mypalette_blue[3,], mypalette_blue[7,], mypalette_blue[13,], mypalette_blue[18,])
image(1:4, 1, as.matrix(1:4), col = mypalette_category, main = "Categories colors", 
      xlab =  "1-early 2-late/new 3-immature 4-adult", cex.lab = 1.3, cex.main =2,
      ylab = "", yaxt = "n")

#Palette for groups (Mysticeti, Odontoceti)
mypalette_groups <- c(mypalette_myst2[4], mypalette_odont2[6])
image(1:2, 1, as.matrix(1:2), col = mypalette_groups, main = "Groups colors", xlab = "1-Mysticeti 2-Odontoceti", 
      cex.lab = 1.2, cex.main =2,ylab = "", yaxt = "n", xaxt = "n")
axis(1, at = c(1, 2))


#Create shape palette 2 groups (Mysticeti, Odontoceti), 4 categories and 13 families
shapes <- c(21,22)
shapes_cat <- c(23,24,22,21)
shapes_fam <- c(17, 19, 0, 7, 13, 2, 4, 15, 9, 5, 3, 6, 11) #mysticeti solid, odontoceti empty

##Images for plots
myst <- readPNG("Data/megaptera.png")
odont <- readPNG("Data/stenella.png")

#GPA ALIGNMENT ----

##Whole skull ----
#Procrustes alignment, should also show mean config coordinates
gpa <- gpagen(shape_array) 

#Save Centroid size as object
Csize <- gpa$Csize 
#Log-transform Centroid size as object
logCsize <- log10(Csize) 

#Save mean shape to create links
mean_shape <- gpa$consensus 

#Coordinates of all specimens after GPA alignment
coords <- gpa$coords 

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords, mean = TRUE, label = F, plot.param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))
#Save screenshot of 3D viewer
rgl.snapshot(filename = "Output/plot_gpa_points.png") 
rgl.snapshot(filename = "Output/plot_gpa_points1.png") 

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords)
#Plot landmarks from outliers in 3D to check how they look
spheres3d(coords[,,28], r = 0.002)

#checkLM(shape_array, path="", pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = coords,
                           Id = classifiers$code, genus = classifiers$genus2, category = classifiers$category,
                           family = classifiers$family, group = classifiers$group,
                           TL_100 = classifiers$TL_100, BZW_100 = classifiers$BZW_100, 
                           feeding = classifiers$Feeding_BL20, size = logCsize)
glimpse(gdf)

##Rostrum and braincase ----
#Create separate gpa aligment for the 2 parts of the skull
#Avoid alignment problems due to align whole skull
gpa_rostrum <- gpagen(shape_array[rostrum,,])
gpa_braincase <- gpagen(shape_array[braincase,,])

#Create gdf with all data needed
#Use size of WHOLE SKULL for allometry (modules not fully independent)
gdf_rostrum <- geomorph.data.frame(coords = gpa_rostrum$coords, size = logCsize, Id = gdf$Id,
                                   genus = gdf$genus, category = gdf$category, group = gdf$group,
                                   family = gdf$family, feeding = gdf$feeding, grp_cat = group_category_grp)

#Create gdf with all data needed
gdf_braincase <- geomorph.data.frame(coords = gpa_braincase$coords, size = logCsize,  Id = gdf$Id,
                                     genus = gdf$genus, category = gdf$category, group = gdf$group,
                                     family = gdf$family, feeding = gdf$feeding, grp_cat = group_category_grp)


###Clean up environment before proceeding
Rdata_1 <- list(absent_curves, absent_LMs,absentcurve, absentLM, inherit.name=TRUE)
Rdata_2<-list(subsampled.lm,newpts, newpts2,slidedlms, final_dataset, slid.lms, inherit.name=TRUE)

save(Rdata_1, file = "absent_data.RData")
save(Rdata_2, file = "import_arrays.RData")

rm(Rdata_1, Rdata_2, absent_curves, absent_LMs,absentcurve, absentLM,
   subsampled.lm,newpts, newpts2,slidedlms, final_dataset, slid.lms )

#PREPARE WARP MESH  ----

#Find specimen closer to mean, useful to create warp mesh
findMeanSpec(coords) #number below specimen name is the number of the specimen in the array
#If specimen ply not good find closet one by species/age in list
#dimnames(coords)[3]

#Create object containing only that specimen coordinates
warp_specimen <- coords[,,149] #number displayed by findMeanSpec
warp_specimen 

#Import simplified mesh to create warp mesh on
ref_mesh <- vcgImport("Data/refmesh.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(ref_mesh$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(warp_specimen)

#Import simplified ply
#Less faces, no holes or isolated triangles
refmesh_all <- vcgImport("Data/refmesh_myst.ply")

#Define fixed LMs and shape array only with LMs
fixed_LMs <- c(1:64)

shape_array_LM <- shape_array[fixed_LMs,,]

#Plot on surface
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[fixed_LMs,,52], col =  "firebrick", type = "s",
          radius = 8, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(shape_array[-fixed_LMs,,52], col =  "tomato", type = "s",
          radius = 6, aspect = T, main = "mean",axes = F, main = F, fov = 0)
text3d(shape_array_LM[,1,52], shape_array_LM[,2,52], shape_array_LM[,3,52], 
       texts = fixed_LMs, pos = 4, offset = 1, font = 2) #change pos

rgl.snapshot(filename = "Output/landmarks_dorsal.png") 
rgl.snapshot(filename = "Output/landmarks_lateral1.png") 
rgl.snapshot(filename = "Output/landmarks_lateral2.png") 
rgl.snapshot(filename = "Output/landmarks_ventral.png") 
rgl.snapshot(filename = "Output/landmarks_posterior.png") 
rgl.snapshot(filename = "Output/landmarks_anterior.png") 

play3d(spin3d(axis = c(1, 0,0), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(1, 0,0), rpm = 10), duration = 6, movie = "landmarks" ,dir = "Output/")

#####

#Next - ch. 4 - PCAs