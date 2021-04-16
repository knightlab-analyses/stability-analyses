#########################################################################################################
### This is to estimate stablity & MSE using bootstrap on BMI_Lin_2014 dataset #################
#########################################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../../../code_method/cv_method.R')
source('../../../code_method/getStability.R')
source('../../../code_method/stab_data_applications.R')
source('../../../code_method/bootstrap_test_compLasso_rf.R')

#####################################
##### data preparation ##############
#####################################
soil_otu = as.matrix(read.csv("../../88soils/88soils_genus_table.txt", sep='\t', row.names=1)) 
save(soil_otu, file='../../88soils/88soils_genus_table.RData')

load('../../88soils/88soils_genus_table.RData')
count = soil_otu

## filter 1% + add pesudo count 
x <- count[, colMeans(count > 0) >= 0.01] 
x[x == 0] <- 0.5
x <- x/rowSums(x) # relative abundance
taxa <- log(x)
print(paste('number of features:', dim(taxa)[2], sep=':'))

# # metadata
mf <- read.csv("../../88soils/88soils_modified_metadata.txt", sep='\t', row.names=1)
y <- mf$ph[match(rownames(count), rownames(mf))]
print(paste('number of samples:', length(y), sep=':'))

# save processed data
save(y, taxa, file='../../88soils/soil_ph.RData')



####################################
## compositional lasso #############
# ####################################
print('compLasso')
out_compLasso = boot_stab(num_boot = 100, method = 'compLasso',
                dat_file = '../data_application/88soils/soils_ph.RData')

save(out_compLasso, file=paste0(dir, '/soils_ph_compLasso.RData'))

#############################################################
# Lasso with glmnet default lambda sequence #############
#############################################################
print('Lasso')
out_lasso = boot_stab(num_boot = 100, method = 'lasso', lambda.grid = NULL, 
                dat_file = '../data_application/88soils/soils_ph.RData')

save(out_lasso, file=paste0(dir, '/soils_ph_lasso.RData', sep=''))


# ##############################################################
# ## Elastic Net with self defined lambda (NULL not allowed) #############
# ##############################################################
print('Elnet')
out_elnet = boot_stab(num_boot = 100, method = 'elnet', 
                 	  dat_file = '../../88soils/soils_ph.RData')

save(out_elnet, file=paste0(dir, '/soils_ph_elnet.RData', sep=''))


##############################################################
## Random Forest with Altman feature selection #############
##############################################################
print('RF')
out_rf = boot_stab(num_boot = 100, method = 'RF', method.perm='altmann', mtry.grid=seq(30, 60, 5), 
                   dat_file = '../../88soils/soils_ph.RData')

save(out_rf, file=paste0(dir, '/soils_ph_rf.RData', sep=''))

########################################################################
############### double bootstrap: random forest ##################
# ######################################################################
print('double boot RF')
boot_rf = boot_stab_data(num_boot=100, method = 'RF',
	 					 data_file='../../88soils/soils_ph.RData')
save(boot_rf, file=paste0(dir, '/soils_ph_boot_rf.RData'))

########################################################################
############### double bootstrap: compositional lasso ##################
# ######################################################################
print('double boot compLasso')
boot_compLasso = boot_stab_data(num_boot=100, method = 'compLasso',
	 							data_file='../../88soils/soils_ph.RData')
save(boot_compLasso, file=paste0(dir, '/soils_ph_boot_compLasso.RData'))










