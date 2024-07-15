#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
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
max_TL_BZW <- classifiers %>% group_by(genus) %>% summarize(max_TL = max(TL_m), max_BZW = max(BZW_mm))

#Count number of specimens per genus
count_genus <- classifiers %>% count(genus)

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

#Count how many stages and category per genus - useful to check binning used in analyses
count_category <-  classifiers %>% group_by(genus) %>% count(category)

View(count_category) #category more even

#Lists genera and families
genera_list <- levels(as.factor(classifiers$genus))
families_list <- levels(as.factor(classifiers$family))

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


#Check for outliers in raw data shape array, they would be displayed in red
#Might be due to absent bones, check misstable list
plotOutliers(shape_array)
#Plot outliers landmarks to search for possible problems - check file name to find number in classifiers
checkLM(shape_array, path="Data/ply/", pt.size = 5, suffix=".ply", render = "s", begin = 162, point = "s")

#Save txt file of output
sink("Output/3-GPA/outliers_raw.txt")
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
levels(families) <- str_to_title(levels(families))

#Save groups as factor, useful for later analysis
groups <- as.factor(classifiers$group)
#Check how many colors are needed
groups #2
#To capitalize or change spelling
levels(groups) <- str_to_title(levels(groups))


##Create project palettes ----
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


#Palette for 27 genera - based on genus, different shades for myst and odont
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

# Initialize an empty vector to store the results
families_numbered <- character(length(levels(families)))

# Loop through each word and paste the sequential number
for (i in seq_along(levels(families))) {
  families_numbered[i] <- paste(i, levels(families)[i], sep = "-")
}

#Vector with numbered families
families_numbered

#Insert line breaks after every 3 words
families_labels <- sapply(seq(1, length(families_numbered), by = 3), function(i) {
  paste(families_numbered[i:min(i+2, length(families_numbered))], collapse = " ")
})

#Combine the broken labels with line breaks
families_labels  <- paste(families_labels, collapse = "\n")

#Palette for families - to extarct group colors
levels(families)

mypalette_families <- c(mypalette_myst2[4], mypalette_myst2[3], mypalette_odont2[6], mypalette_odont2[14], #"balaenidae"      "balaenopteridae" "delphinidae"     "delphininae" 
                        mypalette_odont2[10], mypalette_odont2[17], mypalette_odont2[20], #"globicephalinae" "inioidea" "monodontidae"    
                        mypalette_myst2[2], mypalette_odont2[3], mypalette_odont2[13], mypalette_odont2[2], #"neobalaenidae" "orcininae"       "phocoenidae"     "physeteroidea"
                        mypalette_odont2[8], mypalette_odont2[16]) #"platanistidae" "ziphiidae"

plot(rep(1,length(levels(families))),col=mypalette_families ,pch=19,cex=3, main = "Families colors", ylab = "", xlab = "" ,cex.main = 2,
     yaxt = "n", xaxt = "n")
title(xlab = families_labels, cex.lab = 1.3, font.lab = 3, line = -3)
text(x = seq_along(1:length(levels(families))), y = 1.05, labels = seq_along(1:length(levels(families))))

# Initialize an empty vector to store the results
categories_numbered <- character(length(levels(categories)))

# Loop through each word and paste the sequential number
for (i in seq_along(levels(categories))) {
  categories_numbered[i] <- paste(i, levels(categories)[i], sep = "-")
}

#Vector with numbered categories
categories_numbered

#Combine the broken labels with line breaks
categories_labels  <- paste(categories_numbered, collapse = " ")

#Palette for categories - early, late/new, immature, adult
mypalette_category <- c(mypalette_blue[3,], mypalette_blue[7,], mypalette_blue[13,], mypalette_blue[18,])
image(1:4, 1, as.matrix(1:4), col = mypalette_category, main = "Categories colors", 
      xlab =  categories_labels, cex.lab = 1.3, cex.main =2,
      ylab = "", yaxt = "n")

#Palette for groups (Mysticeti, Odontoceti)
mypalette_groups <- c(mypalette_families[2], mypalette_families[12])
image(1:2, 1, as.matrix(1:2), col = mypalette_groups, main = "Groups colors", xlab = "1-Mysticeti 2-Odontoceti", 
      cex.lab = 1.2, cex.main =2,ylab = "", yaxt = "n", xaxt = "n")
axis(1, at = c(1, 2))


#Create shape palette 2 groups (Mysticeti, Odontoceti), 4 categories and 13 families
#?pch to see shapes
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
rgl.snapshot(filename = "Output/3-GPA/plot_gpa_points.png") 


#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords)
#Plot landmarks from outliers in 3D to check how they look
#spheres3d(coords[,,28], r = 0.002)

#checkLM(shape_array, path="", pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = coords,
                           Id = classifiers$code, genus = classifiers$genus, category = classifiers$category,
                           family = classifiers$family, group = classifiers$group,
                           TL_100 = classifiers$TL_100, BZW_100 = classifiers$BZW_100, 
                           size = logCsize)
glimpse(gdf)

##Rostrum and braincase ----
#Create separate gpa aligment for the 2 parts of the skull
#Avoid alignment problems due to align whole skull

#First create list all bones
modules_all <- rep('other', dim(gdf$coords)[[1]]) 

#Put selected landmarks in each module
modules_all[premaxilla]<-'premaxilla' 
modules_all[maxilla]<-'maxilla' 
modules_all[nasals]<-'nasal' 
modules_all[orbit]<-'orbit' 
modules_all[squamosal]<-'squamosal' 
modules_all[palatine]<-'palatine' 
modules_all[interparietal]<-'interparietal'
modules_all[supraoccipital]<-'supraoccipital' 
modules_all[exoccipital]<-'exoccipital'
modules_all[condyles]<-'condyles'
modules_all[basioccipital]<-'basioccipital' 
modules_all

#Create rostrum and briancase partitions based on developmental hypothesis in Goswami et al. 2022
#Rostrum -> Neural crest
#Neural crest - maxilla, premaxilla, palatine, vomer, nasal, squamosal, orbit
#Briancase  -> Mesoderm
#Mesoderm - supraoccipital, exoccipital, interparietal, condyle, basioccipital
landmarks <- 1:dim(gdf$coords)[[1]]

rostrum <- landmarks[which(modules_all %in% c("maxilla", "premaxilla", "palatine", "nasal", "squamosal", "orbit"))]  

braincase <- landmarks[-rostrum]

#Plot on surface to check assigment
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[rostrum,,41], col =  mypalette_paired[5], type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(shape_array[braincase,,41], col =  mypalette_paired[1], type = "s",
          radius = 0.7, aspect = T, main = "mean",axes = F, main = F, fov = 0)

#Extract modules from general alignment
coords_rostrum <- gdf$coords[rostrum,,]
coords_braincase <- gdf$coords[braincase,,]

#Create gdf with all data needed
#Use size of WHOLE SKULL for allometry (modules not fully independent)
gdf_rostrum <- geomorph.data.frame(coords = coords_rostrum, size = logCsize, Id = gdf$Id,
                                   genus = gdf$genus, category = gdf$category, group = gdf$group,
                                   family = gdf$family)

#Create gdf with all data needed
gdf_braincase <- geomorph.data.frame(coords = coords_braincase, size = logCsize,  Id = gdf$Id,
                                     genus = gdf$genus, category = gdf$category, group = gdf$group,
                                     family = gdf$family)


#PLOT LANDMARKS ON SURFACES ----
##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(shape_array)
#Simplify ply and save in Data folder

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
clear3d()

rgl.snapshot(filename = "Output/3-GPA/landmarks_dorsal.png") 
rgl.snapshot(filename = "Output/3-GPA/landmarks_lateral1.png") 
rgl.snapshot(filename = "Output/3-GPA/landmarks_lateral2.png") 
rgl.snapshot(filename = "Output/3-GPA/landmarks_ventral.png") 
rgl.snapshot(filename = "Output/3-GPA/landmarks_posterior.png") 
rgl.snapshot(filename = "Output/3-GPA/landmarks_anterior.png") 

play3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6, movie = "landmarks" ,dir = "Output/3-GPA/")


#Landmarks plotted on early fetus and adult groups
#Select specimens from PCA plot and prepare mesh file
specimens[match("Ff3", Ids)]
specimens[match("Sa1", Ids)]

specimens[match("Ttrf3", Ids)]
specimens[match("Tada1", Ids)]

#Import simplified ply
myst_fetus <- vcgImport("Data/myst_fetus.ply")
myst_adult <- vcgImport("Data/myst_adult.ply")
odont_fetus <- vcgImport("Data/odont_fetus.ply")
odont_adult <- vcgImport("Data/odont_adult.ply")

#Plot on surface
shade3d(myst_fetus, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Ff3", Ids)], col =  col_modules, type = "s",
          radius = 0.8, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/3-GPA/myst_fetus.png") 
clear3d()

shade3d(myst_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Sa1", Ids)], col =  col_modules, type = "s",
          radius = 10, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/3-GPA/myst_adult.png") 
clear3d()

shade3d(odont_fetus, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Ttrf3", Ids)], col =  col_modules, type = "s",
          radius = 0.5, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/3-GPA/odont_fetus.png") 
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/3-GPA/odont_adult.png") 
clear3d()


#####

#Next - ch. 4 - Modularity analyses