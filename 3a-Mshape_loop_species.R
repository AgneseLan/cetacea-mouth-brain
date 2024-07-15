#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.3a - Create mean shapes for analyses 

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

jpeg("Output/1-2-3/colors families.jpg", width = 1000, height = 1000)
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
jpeg("Output/1-2-3/colors categories.jpg", width = 1000, height = 1000)
image(1:4, 1, as.matrix(1:4), col = mypalette_category, main = "Categories colors", 
      xlab =  categories_labels, cex.lab = 1.3, cex.main =2,
      ylab = "", yaxt = "n")
dev.off()

#Palette for groups (Mysticeti, Odontoceti)
mypalette_groups <- c(mypalette_families[2], mypalette_families[12])
jpeg("Output/1-2-3/colors groups.jpg", width = 1000, height = 1000)
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

#================================================================#
#             MSHAPE LOOP FOR EACH GENUS IN EACH CATEGORY        #
#================================================================#

#MEAN SHAPES GENUS BY CATEGORY----

#Create common data frame with data
means_dataframe <- geomorph.data.frame(raw_coords = shape_array, genus = classifiers$genus, category = classifiers$category) 

#1 - Select coords for each category
#Create rows list first
categories_list <- levels(means_dataframe$category)

#Empty object
rows_categories <- list()

#Loop
for (c in 1:length(categories_list)){
  rows_categories[[c]] <- which(means_dataframe$category == categories_list[c])
}

#Empty object
category_coords <- list()

#Loop
for (s in 1:length(rows_categories)){
  category_coords[[s]] <- means_dataframe$raw_coords[,,rows_categories[[s]]]
}

#2 - Select genera for each category
#Empty object
category_genera <- list()

#Loop
for (s in 1:s){
  category_genera[[s]] <- means_dataframe$genus[rows_categories[[s]]]
}

#3 - Combine coords and genera in new dataframe per category
#Empty object
means_dataframe_category <- list()

#Loop
for (s in 1:s){
  means_dataframe_category[[s]] <- geomorph.data.frame(raw_coords = category_coords[[s]], genus = category_genera[[s]])
}

#4 - Select rows for each genus in each category
#Create list for each category first

#Empty object
genera_categories_list <- list()

#Loop
for (s in 1:s){
genera_categories_list[[s]] <- levels(as.factor(means_dataframe_category[[s]]$genus))
}

#Empty object
category_genus_rows_1 <- list()
category_genus_rows_2 <- list()
category_genus_rows_3 <- list()
category_genus_rows_4 <- list()

#Loop
for (k in 1:length(genera_categories_list[[1]])){
  category_genus_rows_1[[k]] <- which(means_dataframe_category[[1]]$genus == genera_categories_list[[1]][k])
  }

#Loop
for (k in 1:length(genera_categories_list[[2]])){
  category_genus_rows_2[[k]] <- which(means_dataframe_category[[2]]$genus == genera_categories_list[[2]][k])
}

#Loop
for (k in 1:length(genera_categories_list[[3]])){
  category_genus_rows_3[[k]] <- which(means_dataframe_category[[3]]$genus == genera_categories_list[[3]][k])
}

#Loop
for (k in 1:length(genera_categories_list[[4]])){
  category_genus_rows_4[[k]] <- which(means_dataframe_category[[4]]$genus == genera_categories_list[[4]][k])
}

#5 - Make arrays with mshape for each genus in each category
#Empty object
#First number = number of landmark points, second number = 3 dimensions, third number = number of category
coords_early <- array(0, dim = c(dim(means_dataframe$raw_coords)[[1]],3,length(genera_categories_list[[1]]))) 
coords_late_new <- array(0, dim = c(dim(means_dataframe$raw_coords)[[1]],3,length(genera_categories_list[[2]])))
coords_immature <- array(0, dim = c(dim(means_dataframe$raw_coords)[[1]],3,length(genera_categories_list[[3]])))
coords_adult <- array(0, dim = c(dim(means_dataframe$raw_coords)[[1]],3,length(genera_categories_list[[4]])))

#Loop
for (k in 1:length(genera_categories_list[[1]])){
  coords_early[,,k] <- mshape(means_dataframe_category[[1]]$raw_coords[,,category_genus_rows_1[[k]]])
  }

#Loop
for (k in 1:length(genera_categories_list[[2]])){
  coords_late_new[,,k] <- mshape(means_dataframe_category[[2]]$raw_coords[,,category_genus_rows_2[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[3]])){
  coords_immature[,,k] <- mshape(means_dataframe_category[[3]]$raw_coords[,,category_genus_rows_3[[k]]])
}

#Loop
for (k in 1:length(genera_categories_list[[4]])){
  coords_adult[,,k] <- mshape(means_dataframe_category[[4]]$raw_coords[,,category_genus_rows_4[[k]]])
}

#Assign dimnames as genera names
dimnames(coords_early) [[3]] <- as.list(genera_categories_list[[1]])
dimnames(coords_late_new) [[3]] <- as.list(genera_categories_list[[2]])
dimnames(coords_immature) [[3]] <- as.list(genera_categories_list[[3]])
dimnames(coords_adult) [[3]] <- as.list(genera_categories_list[[4]])


##Dataframes mean shapes by category ----

#Create grouping variable
genera_list <- levels(as.factor(means_dataframe$genus))

groups_early <- if_else(dimnames(coords_early)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_late_new <- if_else(dimnames(coords_late_new)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_immature <- if_else(dimnames(coords_immature)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")
groups_adult <- if_else(dimnames(coords_adult)[[3]] %in%  genera_list[c(1:4,7,13)], "mysticeti", "odontoceti")

#Create means_dataframe mean shapes 
means_dataframe_early <- geomorph.data.frame(coords = coords_early, genus = dimnames(coords_early)[[3]], group = groups_early, category = rep(categories_list[1], times = length(groups_early)))
means_dataframe_late_new <- geomorph.data.frame(coords = coords_late_new, genus = dimnames(coords_late_new)[[3]], group = groups_late_new, category = rep(categories_list[2], times = length(groups_late_new)))
means_dataframe_immature <- geomorph.data.frame(coords = coords_immature, genus = dimnames(coords_immature)[[3]], group = groups_immature, category = rep(categories_list[3], times = length(groups_immature)))
means_dataframe_adult <- geomorph.data.frame(coords = coords_adult, genus = dimnames(coords_adult)[[3]], group = groups_adult, category = rep(categories_list[4], times = length(groups_adult)))

#List means_dataframes
means_dataframe_shapes <- list(means_dataframe_early, means_dataframe_late_new, means_dataframe_immature, means_dataframe_adult)



#Create means_dataframes for groups
#Find rows for Mysticeti
rows_mysticeti_means <- list()

for (c in 1:length(categories_list)){  
  rows_mysticeti_means[[c]] <-which(means_dataframe_mean_shapes[[c]][["group"]] == "mysticeti")
}  

#Create means_dataframe mean shapes
means_dataframe_mean_early_myst <- geomorph.data.frame(coords = coords_early[,,rows_mysticeti_means[[1]]], genus = dimnames(coords_early[,,rows_mysticeti_means[[1]]])[[3]], 
                                           size = logCsize_early[rows_mysticeti_means[[1]]], group = groups_early[rows_mysticeti_means[[1]]], category = rep(categories_list[1], times = length(logCsize_early[rows_mysticeti_means[[1]]])))
means_dataframe_mean_late_new_myst <- geomorph.data.frame(coords = coords_late_new[,,rows_mysticeti_means[[2]]], genus = dimnames(coords_late_new[,,rows_mysticeti_means[[2]]])[[3]], 
                                              size = logCsize_late_new[rows_mysticeti_means[[2]]], group = groups_late_new[rows_mysticeti_means[[2]]], category = rep(categories_list[2], times = length(logCsize_late_new[rows_mysticeti_means[[2]]])))
means_dataframe_mean_immature_myst <- geomorph.data.frame(coords = coords_immature[,,rows_mysticeti_means[[3]]], genus = dimnames(coords_immature[,,rows_mysticeti_means[[3]]])[[3]], 
                                              size = logCsize_immature[rows_mysticeti_means[[3]]], group = groups_immature[rows_mysticeti_means[[3]]], category = rep(categories_list[3], times = length(logCsize_immature[rows_mysticeti_means[[3]]])))
means_dataframe_mean_adult_myst <- geomorph.data.frame(coords = coords_adult[,,rows_mysticeti_means[[4]]], genus = dimnames(coords_adult[,,rows_mysticeti_means[[4]]])[[3]], 
                                           size = logCsize_adult[rows_mysticeti_means[[4]]], group = groups_adult[rows_mysticeti_means[[4]]], category = rep(categories_list[4], times = length(logCsize_adult[rows_mysticeti_means[[4]]])))

#List means_dataframes
means_dataframe_mean_shapes_myst <- list(means_dataframe_mean_early_myst, means_dataframe_mean_late_new_myst, means_dataframe_mean_immature_myst, means_dataframe_mean_adult_myst)

#Find rows for Odontoceti
rows_odontoceti_means <- list()

for (c in 1:length(categories_list)){  
  rows_odontoceti_means[[c]] <-which(means_dataframe_mean_shapes[[c]][["group"]] == "odontoceti")
}  

#Create means_dataframe mean shapes
means_dataframe_mean_early_odont <- geomorph.data.frame(coords = coords_early[,,rows_odontoceti_means[[1]]], genus = dimnames(coords_early[,,rows_odontoceti_means[[1]]])[[3]], 
                                           size = logCsize_early[rows_odontoceti_means[[1]]], group = groups_early[rows_odontoceti_means[[1]]], category = rep(categories_list[1], times = length(logCsize_early[rows_odontoceti_means[[1]]])))
means_dataframe_mean_late_new_odont <- geomorph.data.frame(coords = coords_late_new[,,rows_odontoceti_means[[2]]], genus = dimnames(coords_late_new[,,rows_odontoceti_means[[2]]])[[3]], 
                                              size = logCsize_late_new[rows_odontoceti_means[[2]]], group = groups_late_new[rows_odontoceti_means[[2]]], category = rep(categories_list[2], times = length(logCsize_late_new[rows_odontoceti_means[[2]]])))
means_dataframe_mean_immature_odont <- geomorph.data.frame(coords = coords_immature[,,rows_odontoceti_means[[3]]], genus = dimnames(coords_immature[,,rows_odontoceti_means[[3]]])[[3]], 
                                              size = logCsize_immature[rows_odontoceti_means[[3]]], group = groups_immature[rows_odontoceti_means[[3]]], category = rep(categories_list[3], times = length(logCsize_immature[rows_odontoceti_means[[3]]])))
means_dataframe_mean_adult_odont <- geomorph.data.frame(coords = coords_adult[,,rows_odontoceti_means[[4]]], genus = dimnames(coords_adult[,,rows_odontoceti_means[[4]]])[[3]], 
                                           size = logCsize_adult[rows_odontoceti_means[[4]]], group = groups_adult[rows_odontoceti_means[[4]]], category = rep(categories_list[4], times = length(logCsize_adult[rows_odontoceti_means[[4]]])))

#List means_dataframes
means_dataframe_mean_shapes_odont <- list(means_dataframe_mean_early_odont, means_dataframe_mean_late_new_odont, means_dataframe_mean_immature_odont, means_dataframe_mean_adult_odont)

