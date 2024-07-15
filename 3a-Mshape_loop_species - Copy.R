#================================================================#
#             MSHAPE LOOP FOR EACH GENUS IN EACH CATEGORY        #
#================================================================#

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

#Calculate mean CS for each genus in each category for analysis
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


#GDF mean shapes ----
#Create gdf mean shapes 
gdf_mean_early <- geomorph.data.frame(coords = coords_early, genus = dimnames(coords_early)[[3]], size = logCsize_early, group = groups_early, family = families_categories_list[[1]], category = rep(categories_list[1], times = length(logCsize_early)))
gdf_mean_late_new <- geomorph.data.frame(coords = coords_late_new, genus = dimnames(coords_late_new)[[3]], size = logCsize_late_new, group = groups_late_new, family = families_categories_list[[2]], category = rep(categories_list[2], times = length(logCsize_late_new)))
gdf_mean_immature <- geomorph.data.frame(coords = coords_immature, genus = dimnames(coords_immature)[[3]], size = logCsize_immature, group = groups_immature, family = families_categories_list[[3]], category = rep(categories_list[3], times = length(logCsize_immature)))
gdf_mean_adult <- geomorph.data.frame(coords = coords_adult, genus = dimnames(coords_adult)[[3]], size = logCsize_adult, group = groups_adult, family = families_categories_list[[4]], category = rep(categories_list[4], times = length(logCsize_adult)))

#List gdfs
gdf_mean_shapes <- list(gdf_mean_early, gdf_mean_late_new, gdf_mean_immature, gdf_mean_adult)

#Create gdfs for groups
#Find rows for Mysticeti
rows_mysticeti_means <- list()

for (c in 1:length(categories_list)){  
  rows_mysticeti_means[[c]] <-which(gdf_mean_shapes[[c]][["group"]] == "mysticeti")
}  

#Create gdf mean shapes
gdf_mean_early_myst <- geomorph.data.frame(coords = coords_early[,,rows_mysticeti_means[[1]]], genus = dimnames(coords_early[,,rows_mysticeti_means[[1]]])[[3]], 
                                           size = logCsize_early[rows_mysticeti_means[[1]]], group = groups_early[rows_mysticeti_means[[1]]], category = rep(categories_list[1], times = length(logCsize_early[rows_mysticeti_means[[1]]])))
gdf_mean_late_new_myst <- geomorph.data.frame(coords = coords_late_new[,,rows_mysticeti_means[[2]]], genus = dimnames(coords_late_new[,,rows_mysticeti_means[[2]]])[[3]], 
                                              size = logCsize_late_new[rows_mysticeti_means[[2]]], group = groups_late_new[rows_mysticeti_means[[2]]], category = rep(categories_list[2], times = length(logCsize_late_new[rows_mysticeti_means[[2]]])))
gdf_mean_immature_myst <- geomorph.data.frame(coords = coords_immature[,,rows_mysticeti_means[[3]]], genus = dimnames(coords_immature[,,rows_mysticeti_means[[3]]])[[3]], 
                                              size = logCsize_immature[rows_mysticeti_means[[3]]], group = groups_immature[rows_mysticeti_means[[3]]], category = rep(categories_list[3], times = length(logCsize_immature[rows_mysticeti_means[[3]]])))
gdf_mean_adult_myst <- geomorph.data.frame(coords = coords_adult[,,rows_mysticeti_means[[4]]], genus = dimnames(coords_adult[,,rows_mysticeti_means[[4]]])[[3]], 
                                           size = logCsize_adult[rows_mysticeti_means[[4]]], group = groups_adult[rows_mysticeti_means[[4]]], category = rep(categories_list[4], times = length(logCsize_adult[rows_mysticeti_means[[4]]])))

#List gdfs
gdf_mean_shapes_myst <- list(gdf_mean_early_myst, gdf_mean_late_new_myst, gdf_mean_immature_myst, gdf_mean_adult_myst)

#Find rows for Odontoceti
rows_odontoceti_means <- list()

for (c in 1:length(categories_list)){  
  rows_odontoceti_means[[c]] <-which(gdf_mean_shapes[[c]][["group"]] == "odontoceti")
}  

#Create gdf mean shapes
gdf_mean_early_odont <- geomorph.data.frame(coords = coords_early[,,rows_odontoceti_means[[1]]], genus = dimnames(coords_early[,,rows_odontoceti_means[[1]]])[[3]], 
                                           size = logCsize_early[rows_odontoceti_means[[1]]], group = groups_early[rows_odontoceti_means[[1]]], category = rep(categories_list[1], times = length(logCsize_early[rows_odontoceti_means[[1]]])))
gdf_mean_late_new_odont <- geomorph.data.frame(coords = coords_late_new[,,rows_odontoceti_means[[2]]], genus = dimnames(coords_late_new[,,rows_odontoceti_means[[2]]])[[3]], 
                                              size = logCsize_late_new[rows_odontoceti_means[[2]]], group = groups_late_new[rows_odontoceti_means[[2]]], category = rep(categories_list[2], times = length(logCsize_late_new[rows_odontoceti_means[[2]]])))
gdf_mean_immature_odont <- geomorph.data.frame(coords = coords_immature[,,rows_odontoceti_means[[3]]], genus = dimnames(coords_immature[,,rows_odontoceti_means[[3]]])[[3]], 
                                              size = logCsize_immature[rows_odontoceti_means[[3]]], group = groups_immature[rows_odontoceti_means[[3]]], category = rep(categories_list[3], times = length(logCsize_immature[rows_odontoceti_means[[3]]])))
gdf_mean_adult_odont <- geomorph.data.frame(coords = coords_adult[,,rows_odontoceti_means[[4]]], genus = dimnames(coords_adult[,,rows_odontoceti_means[[4]]])[[3]], 
                                           size = logCsize_adult[rows_odontoceti_means[[4]]], group = groups_adult[rows_odontoceti_means[[4]]], category = rep(categories_list[4], times = length(logCsize_adult[rows_odontoceti_means[[4]]])))

#List gdfs
gdf_mean_shapes_odont <- list(gdf_mean_early_odont, gdf_mean_late_new_odont, gdf_mean_immature_odont, gdf_mean_adult_odont)

