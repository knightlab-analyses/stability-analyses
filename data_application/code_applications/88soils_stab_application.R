#########################################################################################################
### This is to estimate stablity & MSE using bootstrap on BMI_Lin_2014 dataset #################
#########################################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('cv_method.R')
source('getStability.R')
source('stab_data_applications.R')
source('bootstrap_test_compLasso_rf.R')

#####################################
##### data preparation ##############
#####################################
load('../data_application/88soils/88soils_genus_table.RData')
count = soil_otu

## filter 1% + add pesudo count 
x <- count[, colMeans(count > 0) >= 0.01] 
x[x == 0] <- 0.5
x <- x/rowSums(x) # relative abundance
taxa <- log(x)
print(paste('number of features:', dim(taxa)[2], sep=':'))

# # metadata
mf <- read.csv("../data_application/88soils/88soils_modified_metadata.txt", sep='\t', row.names=1)
y <- mf$ph[match(rownames(count), rownames(mf))]
print(paste('number of samples:', length(y), sep=':'))

# save processed data
save(y, taxa, file='../data_application/88soils/soil_ph.RData')


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
                 	  dat_file = '../data_application/88soils/soils_ph.RData')

save(out_elnet, file=paste0(dir, '/soils_ph_elnet.RData', sep=''))


##############################################################
## Random Forest with Altman feature selection #############
##############################################################
print('RF')
out_rf = boot_stab(num_boot = 100, method = 'RF', method.perm='altmann', mtry.grid=seq(30, 60, 5), 
                   dat_file = '../data_application/88soils/soils_ph.RData')

save(out_rf, file=paste0(dir, '/soils_ph_rf.RData', sep=''))

########################################################################
############### double bootstrap: random forest ##################
# ######################################################################
print('double boot RF')
boot_rf = boot_stab_data(num_boot=100, method = 'RF',
	 					 data_file='../data_application/88soils/soils_ph.RData')
save(boot_rf, file=paste0(dir, '/soils_ph_boot_rf.RData'))

########################################################################
############### double bootstrap: compositional lasso ##################
# ######################################################################
print('double boot compLasso')
boot_compLasso = boot_stab_data(num_boot=100, method = 'compLasso',
	 							data_file='../data_application/88soils/soils_ph.RData')
save(boot_compLasso, file=paste0(dir, '/soils_ph_boot_compLasso.RData'))










