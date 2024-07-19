#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.3 - Prepare final dataset for analysis, run GPA, calculate mean shapes

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
library(abind)

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

#Count how many stages and category per genus - useful to check binning used in analyses
count_category <-  classifiers %>% group_by(genus) %>% count(category)

View(count_category)

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
checkLM(shape_array, path="Data/ply/", pt.size = 5, suffix=".ply", render = "s", begin = 169, point = "s")

#Save txt file of output
sink("Output/1-2-3/outliers_raw.txt")
print(plotOutliers(shape_array))
sink() 
#Mark where outliers end based on plot
#Save plot
jpeg("Output/1-2-3/outliers raw data.jpg", width = 875, height = 860, units = "px", pointsize = 16, quality = 100)
plotOutliers(shape_array)
dev.off()

#Save specimens names as object
specimens <- dimnames(shape_array)[[3]]

#Save specimens Ids as object
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

jpeg("Output/1-2-3/colors families.jpg", width = 875, height = 860, units = "px", pointsize = 16, quality = 100)
plot(rep(1,length(levels(families))),col=mypalette_families ,pch=19,cex=3, main = "Families colors", ylab = "", xlab = "" ,cex.main = 2,
     yaxt = "n", xaxt = "n")
title(xlab = families_labels, cex.lab = 1.3, font.lab = 3, line = -3)
text(x = seq_along(1:length(levels(families))), y = 1.05, labels = seq_along(1:length(levels(families))))
dev.off()

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
jpeg("Output/1-2-3/colors categories.jpg", width = 875, height = 860, units = "px", pointsize = 16, quality = 100)
image(1:4, 1, as.matrix(1:4), col = mypalette_category, main = "Categories colors", 
      xlab =  categories_labels, cex.lab = 1.3, cex.main =2,
      ylab = "", yaxt = "n")
dev.off()

#Palette for groups (Mysticeti, Odontoceti)
mypalette_groups <- c(mypalette_families[2], mypalette_families[12])
jpeg("Output/1-2-3/colors groups.jpg",width = 875, height = 860, units = "px", pointsize = 16, quality = 100)
image(1:2, 1, as.matrix(1:2), col = mypalette_groups, main = "Groups colors", xlab = "1-Mysticeti 2-Odontoceti", 
      cex.lab = 1.2, cex.main =2,ylab = "", yaxt = "n", xaxt = "n")
axis(1, at = c(1, 2))
dev.off()

#Create shape palette 2 groups (Mysticeti, Odontoceti), 4 categories and 13 families
#?pch to see shapes
shapes <- c(21,22)
shapes_cat <- c(23,24,22,21)
shapes_fam <- c(17, 19, 0, 7, 13, 2, 4, 15, 9, 5, 3, 6, 11) #mysticeti solid, odontoceti empty

##Images for plots
myst <- readPNG("Data/megaptera.png")
odont <- readPNG("Data/stenella.png")

#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa <- gpagen(shape_array) 

#Log-transform Centroid size as object
logCsize <- log10(gpa$Csize) 

#Coordinates of all specimens after GPA alignment
coords <- gpa$coords 

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords, mean = TRUE, label = F, plot_param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))
#Save screenshot of 3D viewer
rgl.snapshot(filename = "Output/1-2-3/plot_gpa_points.png") 

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords)
#Plot landmarks from outliers in 3D to check how they look
#spheres3d(coords[,,148], r = 0.002)

#Save plot
jpeg("Output/1-2-3/outliers gpa aligned data.jpg", width = 875, height = 860, units = "px", pointsize = 16, quality = 100)
plotOutliers(coords)
dev.off()

#checkLM(shape_array, path="", pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = coords,
                           genus = classifiers$genus, category = classifiers$category,
                           family = classifiers$family, group = classifiers$group,
                           size = logCsize)
glimpse(gdf)

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


rgl.snapshot(filename = "Output/1-2-3/landmarks_dorsal.png") 
rgl.snapshot(filename = "Output/1-2-3/landmarks_lateral1.png") 
rgl.snapshot(filename = "Output/1-2-3/landmarks_lateral2.png") 
rgl.snapshot(filename = "Output/1-2-3/landmarks_ventral.png") 
rgl.snapshot(filename = "Output/1-2-3/landmarks_posterior.png") 
rgl.snapshot(filename = "Output/1-2-3/landmarks_anterior.png") 

play3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(0, 0,1), rpm = 10), duration = 6, movie = "landmarks" ,dir = "Output/1-2-3/")

#clear3d()

#Landmarks plotted on early fetus and adult groups
#Representative specimens chosen from dataset
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

rgl.snapshot(filename = "Output/1-2-3/myst_fetus.png") 
clear3d()

shade3d(myst_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Sa1", Ids)], col =  col_modules, type = "s",
          radius = 10, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/1-2-3/myst_adult.png") 
clear3d()

shade3d(odont_fetus, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Ttrf3", Ids)], col =  col_modules, type = "s",
          radius = 0.5, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/1-2-3/odont_fetus.png") 
clear3d()

shade3d(odont_adult, col = "white", alpha = 0.5)
spheres3d(shape_array[,,match("Tada1", Ids)], col =  col_modules, type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/1-2-3/odont_adult.png") 
clear3d()

#================================================================#
#             MSHAPE LOOP FOR EACH GENUS IN EACH CATEGORY        #
#================================================================#

#MEAN SHAPES GENUS BY CATEGORY----

#Create mean shapes for each category ----
#1 - Select coords for each category
#Create rows list first
categories_list <- levels(gdf$category)

#Empty object
rows_categories <- list()

#Loop
for (c in 1:length(categories_list)){
  rows_categories[[c]] <- which(gdf$category == categories_list[c])
}


#Empty object
category_coords <- list()

#Loop
for (s in 1:length(rows_categories)){
  category_coords[[s]] <- gdf$coords[,,rows_categories[[s]]]
}

#2 - Select genera for each category
#Empty object
category_genera <- list()

#Loop
for (s in 1:s){
  category_genera[[s]] <- gdf$genus[rows_categories[[s]]]
}

#3 - Combine coords and genera in new gdf per category
#Empty object
gdf_category <- list()

#Loop
for (s in 1:s){
  gdf_category[[s]] <- geomorph.data.frame(coords = category_coords[[s]], genus = category_genera[[s]])
}

#4 - Select rows for each genus in each category
#Create list for each category first

#Empty object
genera_categories_list <- list()

#Loop
for (s in 1:s){
  genera_categories_list[[s]] <- levels(as.factor(gdf_category[[s]]$genus))
}

#Empty object
category_genus_rows_1 <- list()
category_genus_rows_2 <- list()
category_genus_rows_3 <- list()
category_genus_rows_4 <- list()

#Loop
for (k in 1:length(genera_categories_list[[1]])){
  category_genus_rows_1[[k]] <- which(gdf_category[[1]]$genus == genera_categories_list[[1]][k])
}

#Loop
for (k in 1:length(genera_categories_list[[2]])){
  category_genus_rows_2[[k]] <- which(gdf_category[[2]]$genus == genera_categories_list[[2]][k])
}

#Loop
for (k in 1:length(genera_categories_list[[3]])){
  category_genus_rows_3[[k]] <- which(gdf_category[[3]]$genus == genera_categories_list[[3]][k])
}

#Loop
for (k in 1:length(genera_categories_list[[4]])){
  category_genus_rows_4[[k]] <- which(gdf_category[[4]]$genus == genera_categories_list[[4]][k])
}

#5 - Make arrays with mshape for each genus in each category
#Empty object
#First number = number of landmark points, second number = 3 dimensions, third number = number of category
coords_early <- array(0, dim = c(dim(gdf$coords)[[1]],3,length(genera_categories_list[[1]]))) 
coords_late_new <- array(0, dim = c(dim(gdf$coords)[[1]],3,length(genera_categories_list[[2]])))
coords_immature <- array(0, dim = c(dim(gdf$coords)[[1]],3,length(genera_categories_list[[3]])))
coords_adult <- array(0, dim = c(dim(gdf$coords)[[1]],3,length(genera_categories_list[[4]])))

#Loop
for (k in 1:length(genera_categories_list[[1]])){
  coords_early[,,k] <- mshape(gdf_category[[1]]$coords[,,category_genus_rows_1[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[2]])){
  coords_late_new[,,k] <- mshape(gdf_category[[2]]$coords[,,category_genus_rows_2[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[3]])){
  coords_immature[,,k] <- mshape(gdf_category[[3]]$coords[,,category_genus_rows_3[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[4]])){
  coords_adult[,,k] <- mshape(gdf_category[[4]]$coords[,,category_genus_rows_4[[k]]])
}

#Assign dimnames as genera names
dimnames(coords_early) [[3]] <- as.list(genera_categories_list[[1]])
dimnames(coords_late_new) [[3]] <- as.list(genera_categories_list[[2]])
dimnames(coords_immature) [[3]] <- as.list(genera_categories_list[[3]])
dimnames(coords_adult) [[3]] <- as.list(genera_categories_list[[4]])

#Calculate mean logCS for each genus in each category for analysis
#Empty object
logCsize_early <- list()
logCsize_late_new <- list()
logCsize_immature <- list()
logCsize_adult <- list()

#Loop
for (k in 1:length(genera_categories_list[[1]])){
  logCsize_early[[k]] <- mean(logCsize[category_genus_rows_1[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[2]])){
  logCsize_late_new[[k]] <- mean(logCsize[category_genus_rows_2[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[3]])){
  logCsize_immature[[k]] <- mean(logCsize[category_genus_rows_3[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[4]])){
  logCsize_adult[[k]] <- mean(logCsize[category_genus_rows_4[[k]]])
}

#Unlist
logCsize_early <- unlist(logCsize_early)
logCsize_late_new <- unlist(logCsize_late_new)
logCsize_immature <- unlist(logCsize_immature)
logCsize_adult <- unlist(logCsize_adult)

#Create grouping variable
genera_list <- levels(as.factor(gdf$genus))

groups_early <- if_else(dimnames(coords_early)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_late_new <- if_else(dimnames(coords_late_new)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_immature <- if_else(dimnames(coords_immature)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_adult <- if_else(dimnames(coords_adult)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")

genus_family_df <- data.frame(genus = classifiers$genus, family = classifiers$family, category = classifiers$category)

genus_family_df <- genus_family_df %>%  group_by(category) %>% distinct

families_early <- genus_family_df$family[genus_family_df$category == categories_list[1]]
families_late_new <- genus_family_df$family[genus_family_df$category == categories_list[2]]
families_immature <- genus_family_df$family[genus_family_df$category == categories_list[3]]
families_adult <-   genus_family_df$family[genus_family_df$category == categories_list[4]]
  
#GDF mean shapes ----
#Create gdf mean shapes 
gdf_mean_early <- geomorph.data.frame(coords = coords_early, genus = dimnames(coords_early)[[3]], size = logCsize_early, group = groups_early, family = families_early, category = rep(categories_list[1], times = length(logCsize_early)))
gdf_mean_late_new <- geomorph.data.frame(coords = coords_late_new, genus = dimnames(coords_late_new)[[3]], size = logCsize_late_new, group = groups_late_new, family = families_late_new, category = rep(categories_list[2], times = length(logCsize_late_new)))
gdf_mean_immature <- geomorph.data.frame(coords = coords_immature, genus = dimnames(coords_immature)[[3]], size = logCsize_immature, group = groups_immature, family = families_immature, category = rep(categories_list[3], times = length(logCsize_immature)))
gdf_mean_adult <- geomorph.data.frame(coords = coords_adult, genus = dimnames(coords_adult)[[3]], size = logCsize_adult, group = groups_adult, family = families_adult, category = rep(categories_list[4], times = length(logCsize_adult)))

#List gdfs
gdf_mean_shapes <- list(gdf_mean_early, gdf_mean_late_new, gdf_mean_immature, gdf_mean_adult)

glimpse(gdf_mean_shapes)

#Align all mean shapes from list with another GPA - ensures consistent alignment
#Change coords dimnames so that they are unique
gdf_mean_shapes1 <- gdf_mean_shapes

for (c in 1:length(categories_list)){
  dimnames(gdf_mean_shapes1[[c]]$coords)[[3]] <- paste0(dimnames(gdf_mean_shapes1[[c]]$coords)[[3]], "_", categories_list_short[c])
}

#Create common coords array
coords_means <- abind(gdf_mean_shapes1[[1]]$coords, gdf_mean_shapes1[[2]]$coords, gdf_mean_shapes1[[3]]$coords, gdf_mean_shapes1[[4]]$coords)
glimpse(coords_means)

#Perform alignment
gpa_means <- gpagen(coords_means)

##Make data frame for analyses in geomorph
gdf_mean_all <- geomorph.data.frame(coords = gpa_means$coords,
                                    genus = c(gdf_mean_shapes[[1]]$genus, gdf_mean_shapes[[2]]$genus, gdf_mean_shapes[[3]]$genus, gdf_mean_shapes[[4]]$genus), 
                                    category = c(gdf_mean_shapes[[1]]$category, gdf_mean_shapes[[2]]$category, gdf_mean_shapes[[3]]$category, gdf_mean_shapes[[4]]$category),
                                    family = c(gdf_mean_shapes[[1]]$family, gdf_mean_shapes[[2]]$family, gdf_mean_shapes[[3]]$family, gdf_mean_shapes[[4]]$family), 
                                    group = c(gdf_mean_shapes[[1]]$group, gdf_mean_shapes[[2]]$group, gdf_mean_shapes[[3]]$group, gdf_mean_shapes[[4]]$group),
                                    size = c(gdf_mean_shapes[[1]]$size, gdf_mean_shapes[[2]]$size, gdf_mean_shapes[[3]]$size, gdf_mean_shapes[[4]]$size)) #keep original logCS calculation from raw data
glimpse(gdf_mean_all)

##Create gdfs for groups ----
#Mysticeti
rows_mysticeti_means <- which(gdf_mean_all$group == "mysticeti")

gdf_mean_mysticeti <- geomorph.data.frame(coords = gpa_means$coords[,,rows_mysticeti_means],
                                        genus = gdf_mean_all$genus[rows_mysticeti_means], category = gdf_mean_all$category[rows_mysticeti_means],
                                        family =  gdf_mean_all$family[rows_mysticeti_means], group = gdf_mean_all$group[rows_mysticeti_means],
                                        size = gdf_mean_all$size[rows_mysticeti_means]) #overall skull size
glimpse(gdf_mean_mysticeti)

#Odontoceti
gdf_mean_odontoceti <- geomorph.data.frame(coords = gpa_means$coords[,,-rows_mysticeti_means],
                                          genus = gdf_mean_all$genus[-rows_mysticeti_means], category = gdf_mean_all$category[-rows_mysticeti_means],
                                          family =  gdf_mean_all$family[-rows_mysticeti_means], group = gdf_mean_all$group[-rows_mysticeti_means],
                                          size = gdf_mean_all$size[-rows_mysticeti_means]) #overall skull size
glimpse(gdf_mean_odontoceti)



#####

#Next - ch. 4 - Modularity analyses