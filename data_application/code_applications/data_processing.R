###############################################################################################################
######## save Python output into R rdata file to speed up reading OTU files ######## ######## ######## ######## 
###############################################################################################################


# ################################
# ## 88 soils ####################
# ################################

# soil_otu = as.matrix(read.csv("../data_application/88soils/88soils_genus_table.txt", sep='\t', row.names=1)) 
# save(soil_otu, file='../data_application/88soils/88soils_genus_table.RData')

# ################################
# ## oral age ####################
# ################################

# oral_otu = as.matrix(read.csv("../data_application/oral/oral_genus_table.txt", sep='\t', row.names=1))
# save(oral_otu, file='../data_application/oral/oral_genus_table.RData')

# ################################
# ## central park ################
# ################################

# park_otu = as.matrix(read.csv("../data_application/centralPark/soi_centralPark_genus_table.txt", sep='\t', row.names=1)) 
# save(park_otu, file='../data_application/centralPark/soi_centralPark_genus_table.RData')


# ################################
# ## Vitamin D ####################
# ################################

# vit_otu = as.matrix(read.csv("../data_application/vitaminD/vitaminD_genus_table.txt", sep='\t', row.names=1)) 
# save(vit_otu, file='../data_application/vitaminD/vitaminD_genus_table.RData')

################################
## skin age ####################
################################

# skin_otu = as.matrix(read.csv("../data_application/skin/skin_genus_table.txt", sep='\t', row.names=1))
# save(skin_otu, file='../data_application/skin/skin_genus_table.RData')


