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
count <- as.matrix(read.table("../../code_Lin/cvs/data/combo_count_tab.txt")) 

# filter 1% + add pesudo count 
depth <- sapply(strsplit(colnames(count), "\\."), length)
x <- count[, depth == 6 & colMeans(count > 0) >= 0.01]
x[x == 0] <- 0.5
x <- x/rowSums(x) # relative abundance
taxa <- log(x)
print(paste('number of features:', dim(taxa)[2], sep=':'))

# metadata
demo <- read.delim("../../code_Lin/cvs/data/demographic.txt")
y <- demo$bmi[match(rownames(count), demo$pid)]
print(paste('number of samples:', length(y), sep=':'))

# save processed data
save(y, taxa, file='../../BMI/BMI_Lin_2014.RData')

# ####################################
# ## compositional lasso #############
# ####################################
print('compLasso')
out_compLasso = boot_stab(num_boot = 100, method = 'compLasso',
                dat_file = '../../BMI/BMI_Lin_2014.RData')

save(out_compLasso, file=paste0(dir, '/BMI_compLasso.RData'))

##############################################################
## Lasso with glmnet default lambda sequence #############
##############################################################
print('Lasso')
out_lasso = boot_stab(num_boot = 100, method = 'lasso', lambda.grid = NULL, 
                dat_file = '../../BMI/BMI_Lin_2014.RData')

save(out_lasso, file=paste0(dir, '/BMI_lasso.RData', sep=''))


# ##############################################################
# ## Elastic Net with self defined lambda (NULL not allowed) #############
# ##############################################################
print('Elnet')
out_elnet = boot_stab(num_boot = 100, method = 'elnet', 
                 	  dat_file = '../../BMI/BMI_Lin_2014.RData')

save(out_elnet, file=paste0(dir, '/BMI_elnet.RData', sep=''))


# ##############################################################
# ## Random Forest with Altman feature selection #############
# ##############################################################
print('RF')
out_rf = boot_stab(num_boot = 100, method = 'RF', method.perm='altmann', 
                       dat_file = '../../BMI/BMI_Lin_2014.RData')

save(out_rf, file=paste0(dir, '/BMI_rf.RData', sep=''))


#############################################################
# Random Forest with Janitza feature selection #############
#############################################################
print("RF_JNT")
dir = '../data_application'
out_rf_jnt= boot_stab(num_boot = 100, method = 'RF', method.perm='janitza', 
                      dat_file = '../data_application/BMI_Lin_2014.RData')

save(out_rf_jnt, file=paste0(dir, '/BMI_rf_jnt.RData', sep=''))


########################################################################
############### double bootstrap: random forest ##################
# ######################################################################
print('double boot RF')
boot_rf = boot_stab_data(num_boot=100, method = 'RF',
	 					 data_file='../data_application/BMI_Lin_2014.RData')
save(boot_rf, file=paste0(dir, '/BMI_boot_rf.RData'))

########################################################################
############### double bootstrap: compositional lasso ##################
# ######################################################################
print('double boot compLasso')
boot_compLasso = boot_stab_data(num_boot=100, method = 'compLasso',
	 							data_file='../data_application/BMI_Lin_2014.RData')
save(boot_compLasso, file=paste0(dir, '/BMI_boot_compLasso.RData'))






