#===========================================================#
#                                                           #
#          SKULL MODULARITY - MYSTICETI & ODONTOCETI        #
#                                                           #
#===========================================================#

#CH.2 - Assigning coordinates to landmarks of absent bones
#Code adapted from Ellen Coombs

#LOAD LIBRARIES ----
#always do this first!!
library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library("devtools")
library(SURGE)
library(magick)
library(RColorBrewer) 

#ABSENT BONES ----

#Add the data for absent bones for specific species 
#Both curves and fixed LMs

###SET WD to root from console!! -->

#Import LMs list - curves listed in curve_table
LM_table <-  read_csv("Data/LMs_all.csv")

#Import sets of absent curves and LMs and open file to check what bones have absent data
absent_curves <- read.csv("Data/absent_curves_all.csv")
absent_LMs <- read.csv("Data/absent_LMs_all.csv")
View(absent_curves)
View(absent_LMs)

#Make list of specimens using dimnames
specimens_list <- as.character(unlist(dimnames(slidedlms)[3]))

#Check it is still same as pts list
all.equal(specimens_list, ptslist2)

#Save dimnames to file - check specimens names associated with classifiers
write.table(specimens_list, file = "Output/specimens_list1.csv", sep = ",", col.names = NA)

#Check both absent curves and LMs have specimens matching dataset list
all.equal(absent_curves$specimen, specimens_list)
all.equal(absent_LMs$specimen, specimens_list)

#If not, check strings in Excel using exported csv list and "if" function, then reimport tables

##Absent curves####
#Look for bones with absent curves in curve list - column names
# colnames(absent_curves)
# check curve_table$bone for names of bones
curve_nasal_l <- my_curves$Curves[which(curve_table$bone%in%c("nasal_l"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_nasal_r <- my_curves$Curves[which(curve_table$bone%in%c("nasal_r"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_interparietal <- my_curves$Curves[which(curve_table$bone%in%c("interparietal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_condyle_l <- my_curves$Curves[which(curve_table$bone%in%c("condyle_l"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_condyle_r <- my_curves$Curves[which(curve_table$bone%in%c("condyle_r"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_palatine_l <- my_curves$Curves[which(curve_table$bone%in%c("palatine_l"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_palatine_r <- my_curves$Curves[which(curve_table$bone%in%c("palatine_r"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_basioccipital_ll <- my_curves$Curves[which(curve_table$bone%in%c("basioccipital_ll"))]%>%unlist(.)%>%unique(.)%>%sort(.)
curve_basioccipital_lr <- my_curves$Curves[which(curve_table$bone%in%c("basioccipital_lr"))]%>%unlist(.)%>%unique(.)%>%sort(.)

#Create new object for absent curves
absentcurve <- slidedlms

absentcurve[curve_nasal_l,,41] #specimen number on the end, test if it worked - check specimen number form absent_curves file
absentcurve[curve_basioccipital_ll,,89]
absentcurve[curve_condyle_l,,142]

#Loop to substitute coordinates for semilandmarks in absent bone curves
#Put first landmark of curve in matrix for each bone

#Left nasal 
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$nasal_l[i]))
    absentcurve[curve_nasal_l,c(1:3),i] <- matrix(absentcurve[6,c(1:3),i], nrow = length(curve_nasal_l), ncol=3, byrow=TRUE)
}

#Right nasal
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$nasal_r[i]))
    absentcurve[curve_nasal_r,c(1:3),i] <- matrix(absentcurve[15,c(1:3),i], nrow = length(curve_nasal_r), ncol=3, byrow=TRUE)
}

#Interparietal
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$interparietal[i]))
    absentcurve[curve_interparietal,c(1:3),i] <- matrix(absentcurve[57,c(1:3),i], nrow = length(curve_interparietal), ncol=3, byrow=TRUE)
}

#Condyle left
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$condyle_l[i]))
    absentcurve[curve_condyle_l,c(1:3),i] <- matrix(absentcurve[36,c(1:3),i], nrow = length(curve_condyle_l), ncol=3, byrow=TRUE)
}

#Condyle right
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$condyle_r[i]))
    absentcurve[curve_condyle_r,c(1:3),i] <- matrix(absentcurve[40,c(1:3),i], nrow = length(curve_condyle_r), ncol=3, byrow=TRUE)
}

#Palatine left
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$palatine_l[i]))
    absentcurve[curve_palatine_l,c(1:3),i] <- matrix(absentcurve[51,c(1:3),i], nrow = length(curve_palatine_l), ncol=3, byrow=TRUE)
}

#Palatine right
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$palatine_r[i]))
    absentcurve[curve_palatine_r,c(1:3),i] <- matrix(absentcurve[46,c(1:3),i], nrow = length(curve_palatine_r), ncol=3, byrow=TRUE)
}

#Basioccipital lateral left
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$basioccipital_ll[i]))
    absentcurve[curve_basioccipital_ll,c(1:3),i] <- matrix(absentcurve[56,c(1:3),i], nrow = length(curve_basioccipital_ll), ncol=3, byrow=TRUE)
}

#Basioccipital lateral right
for (i in 1:nrow(absent_curves)){
  if( !is.na(absent_curves$basioccipital_lr[i]))
    absentcurve[curve_basioccipital_lr,c(1:3),i] <- matrix(absentcurve[54,c(1:3),i], nrow = length(curve_basioccipital_lr), ncol=3, byrow=TRUE)
}


absentcurve[curve_nasal_l,,41] #check if it worked with specimen number with missing curve
absentcurve[curve_basioccipital_ll,,89]
absentcurve[curve_condyle_l,,142]

##Absent LMs####
#Look for absent bones first
# colnames(absent_LMs)
# check LM_table$bone for names of bones
LMs_nasal_l <- LM_table$lm[which(LM_table$bone%in%c("nasal_l"))]
LMs_nasal_r <- LM_table$lm[which(LM_table$bone%in%c("nasal_r"))]
LMs_interparietal <- LM_table$lm[which(LM_table$bone%in%c("interparietal"))]
LMs_condyle_l <- LM_table$lm[which(LM_table$bone%in%c("condyle_l"))]
LMs_condyle_r <- LM_table$lm[which(LM_table$bone%in%c("condyle_r"))]
LMs_palatine_l <- LM_table$lm[which(LM_table$bone%in%c("palatine_l"))]
LMs_palatine_r <- LM_table$lm[which(LM_table$bone%in%c("palatine_r"))]
LMs_vomer <- LM_table$lm[which(LM_table$bone%in%c("vomer"))]

#Create new object for absent LMs
absentLM <- absentcurve

absentLM[LMs_nasal_l,,41] #specimen number on the end, test if it worked - check specimen number form absent_LMs file
absentLM[LMs_palatine_l,,146]

#Loop to substitute coordinates for landmarks in absent bones
#Put first landmark in matrix for each bone

#Left nasal
  for (i in 1:nrow(absent_LMs)){
    if( !is.na(absent_LMs$nasal_l[i]))
      absentLM[LMs_nasal_l,c(1:3),i] <- matrix(absentLM[6,c(1:3),i], nrow = length(LMs_nasal_l), ncol=3, byrow=TRUE) #number (40) here is the LM that is missing 
  }

#Right nasal
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$nasal_r[i]))
    absentLM[LMs_nasal_r,c(1:3),i] <- matrix(absentLM[15,c(1:3),i], nrow = length(LMs_nasal_r), ncol=3, byrow=TRUE) 
}

#Interparietal
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$interparietal[i]))
    absentLM[LMs_interparietal,c(1:3),i] <- matrix(absentLM[57,c(1:3),i], nrow = length(LMs_interparietal), ncol=3, byrow=TRUE) 
}

#Condyle left
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$condyle_l[i]))
    absentLM[LMs_condyle_l,c(1:3),i] <- matrix(absentLM[36,c(1:3),i], nrow = length(LMs_condyle_l), ncol=3, byrow=TRUE) 
}

#Condyle right
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$condyle_r[i]))
    absentLM[LMs_condyle_r,c(1:3),i] <- matrix(absentLM[40,c(1:3),i], nrow = length(LMs_condyle_r), ncol=3, byrow=TRUE) 
}

#Palatine left
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$palatine_l[i]))
    absentLM[LMs_palatine_l,c(1:3),i] <- matrix(absentLM[51,c(1:3),i], nrow = length(LMs_palatine_l), ncol=3, byrow=TRUE) 
}

#Palatine right
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$palatine_r[i]))
    absentLM[LMs_palatine_r,c(1:3),i] <- matrix(absentLM[46,c(1:3),i], nrow = length(LMs_palatine_r), ncol=3, byrow=TRUE) 
}

#Vomer
for (i in 1:nrow(absent_LMs)){
  if( !is.na(absent_LMs$vomer[i]))
    absentLM[LMs_vomer,c(1:3),i] <- matrix(absentLM[64,c(1:3),i], nrow = length(LMs_vomer), ncol=3, byrow=TRUE) 
}

absentLM[LMs_nasal_l,,41] #specimen number on the end, test if it worked - check specimen number form absent_LMs file
absentLM[LMs_palatine_l,,146]

absentLM[curve_nasal_l,,41] #check curves still ok

#Create new object for analyses with all missing data, include only shape data
final_dataset <- absentLM

#Check plotting of absent bones
#Look up number for specimens with absent bones in absent_curves and absent_LMs

###SET WD to ply from console!! -->

#Plot
checkLM(final_dataset, path="", pt.size = 7, suffix=".ply", render = "s", begin = 89, point = "p")

###Define modules - List of points and curves for different bones####
maxilla <- c(LM_table$lm[which(LM_table$bone%in%c("maxilla"))], 
             my_curves$Curves[which(curve_table$bone%in%c("maxilla"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
premaxilla <- c(LM_table$lm[which(LM_table$bone%in%c("premaxilla"))], 
             my_curves$Curves[which(curve_table$bone%in%c("premaxilla"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
nasals <- c(LM_table$lm[which(LM_table$bone%in%c("nasal_l", "nasal_r"))], 
            my_curves$Curves[which(curve_table$bone%in%c("nasal_l", "nasal_r"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
orbit <- c(LM_table$lm[which(LM_table$bone%in%c("frontal", "jugal"))], 
           my_curves$Curves[which(curve_table$bone%in%c("frontal"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
squamosal <- c(LM_table$lm[which(LM_table$bone%in%c("squamosal"))], 
               my_curves$Curves[which(curve_table$bone%in%c("squamosal"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
interparietal <-  c(LM_table$lm[which(LM_table$bone%in%c("interparietal"))], 
                    my_curves$Curves[which(curve_table$bone%in%c("interparietal"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
supraoccipital <- c(LM_table$lm[which(LM_table$bone%in%c("supraoccipital"))], 
                    my_curves$Curves[which(curve_table$bone%in%c("supraoccipital"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
exoccipital <- c(LM_table$lm[which(LM_table$bone%in%c("exoccipital"))], 
                 my_curves$Curves[which(curve_table$bone%in%c("exoccipital"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
condyles <- c(LM_table$lm[which(LM_table$bone%in%c("condyle_l", "condyle_r"))], 
              my_curves$Curves[which(curve_table$bone%in%c("condyle_l", "condyle_r"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
basioccipital <- c(LM_table$lm[which(LM_table$bone%in%c("basioccipital", "basioccipital_lr", "basioccipital_ll"))],
                   my_curves$Curves[which(curve_table$bone%in%c("basioccipital", "basioccipital_ll", "basioccipital_lr"))]) %>% 
                    unlist(.)%>%unique(.)%>%sort(.)
palatine <-  c(LM_table$lm[which(LM_table$bone%in%c("palatine_l", "palatine_r", "vomer"))], 
               my_curves$Curves[which(curve_table$bone%in%c("palatine_l", "palatine_r"))]) %>% unlist(.)%>%unique(.)%>%sort(.)

#Check
spheres3d(final_dataset[nasals,,30], radius=1, color = "red") 
spheres3d(final_dataset[-nasals,,30], radius=1, color = "grey")

spheres3d(final_dataset[basioccipital,,30], radius=1, color = "red") 
spheres3d(final_dataset[-basioccipital,,30], radius=1, color = "grey")

####Create color scale for plots#####
#Easier then repeating the code every time
#Define which curves belong to which modules for a pretty plot

###SET WD to root from console!! -->

#Get first palette
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

#Create curve module and LM module objects
#Labels must follow order of landmarks
curvemodules <- c()
for (i in 1:nrow(curve_table)) {
  x<-rep(curve_table$bone[i],curve_table$ptswanted[i])
  curvemodules <- c(curvemodules,x)
}

LMmodules <- LM_table$bone

#Color scale
col_modules <- c(LMmodules, curvemodules)

col_modules <- as.factor(col_modules)

modules_list <- levels(col_modules) #save list of modules for later
modules_list <- str_to_title(modules_list)

#Change to colors
levels(col_modules) <- c(mypalette_paired[5], mypalette_paired[5], mypalette_paired[5], #basioccipital
                         mypalette_paired[10], mypalette_paired[10], #condyles
                         mypalette_paired[9], #exoccipital
                        mypalette_paired[3], mypalette_paired[1], mypalette_paired[3], #frontal, interparietal, jugal
                         mypalette_paired[8], mypalette_paired[4],  mypalette_paired[4], #maxilla, nasals
                        mypalette_paired[7], mypalette_paired[7], mypalette_paired[6], #palatines, premax
                        mypalette_paired[12], mypalette_paired[2], mypalette_paired[7]) #squamosal, socc., vomer

#Plot the landmarks and curves colored by module
spheres3d(final_dataset[,,1], radius = 5, col = col_modules)

#Initialize an empty vector to store the results
modules_numbered <- character(length(modules_list))

#Loop through each word and paste the sequential number
for (i in seq_along(modules_list)) {
  modules_numbered[i] <- paste(i, modules_list[i], sep = "-")
}

#Vector with numbered modules
modules_numbered

#Insert line breaks after every 3 words
modules_labels <- sapply(seq(1, length(modules_numbered), by = 3), function(i) {
  paste(modules_numbered[i:min(i+2, length(modules_numbered))], collapse = " ")
})

#Combine the broken labels with line breaks
modules_labels  <- paste(modules_labels, collapse = "\n")

#Modules color vector
modules_col_vector <- c(mypalette_paired[5], mypalette_paired[5], mypalette_paired[5], #basioccipital
                        mypalette_paired[10], mypalette_paired[10], #condyles
                        mypalette_paired[9], #exoccipital
                        mypalette_paired[3], mypalette_paired[1], mypalette_paired[3], #frontal, interparietal, jugal
                        mypalette_paired[8], mypalette_paired[4],  mypalette_paired[4], #maxilla, nasals
                        mypalette_paired[7], mypalette_paired[7], mypalette_paired[6], #palatines, premax
                        mypalette_paired[12], mypalette_paired[2], mypalette_paired[7]) #squamosal, socc., vomer

jpeg("Output/1-2-3/colors modules.jpg", width = 950, height = 860, units = "px", pointsize = 16, quality = 100)
plot(rep(1,length(modules_list)),col=modules_col_vector,pch=19,cex=3, main = "Modules colors", ylab = "", xlab = "" ,cex.main = 2,
     yaxt = "n", xaxt = "n")
title(xlab = modules_labels, cex.lab = 1.3, font.lab = 1, line = -3)
text(x = seq_along(1:length(modules_list)), y = 1.05, labels = seq_along(1:length(modules_list)))
dev.off()

#Save coordinates to file
save(final_dataset, file = "Output/final_dataset.RData")

###### 
#Next - ch. 3 - GPA and mean shapes
