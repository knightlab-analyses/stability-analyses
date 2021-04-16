#########################################################################################################
### This is to estimate stablity & MSE using bootstrap on BMI_Lin_2014 dataset #################
#########################################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../../../code_method/cv_method_binary_update.R')
source('../../../code_method/getStability.R')
source('../../../code_method/stab_data_applications_binary.R') # boot_stab()
source('../../../code_method/bootstrap_test_compLasso_rf_binary.R') # boot_stab_data()

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
y <- as.factor(ifelse(y >= median(y), 1, 0)) # transform to be binary
print(paste('number of samples:', length(y), sep=':'))

# save processed data
save(y, taxa, file='../../BMI/BMI_Lin_2014_binary.RData')

#------------------------
# compare methods
#------------------------
i <- as.numeric(Sys.getenv("PBS_ARRAYID"))

if (i == 1){
    # ####################################
    # ## Generalized compositional lasso #############
    # ####################################
    print('GencompLasso')
    out_GenCompLasso = boot_stab(num_boot = 100, method = 'GenCompLasso',
                    dat_file = '../../BMI/BMI_Lin_2014_binary.RData')

    save(out_GenCompLasso, file=paste0(dir, '/BMI_binary_GenCompLasso.RData'))
}else if (i == 2){

    ##############################################################
    ## Lasso with glmnet default lambda sequence #############
    ##############################################################
    print('Lasso')
    out_lasso = boot_stab(num_boot = 100, method = 'lasso', lambda.grid = NULL, 
                    dat_file = '../../BMI/BMI_Lin_2014_binary.RData')

    save(out_lasso, file=paste0(dir, '/BMI_binary_lasso.RData', sep=''))
}else if (i == 3){
    # ##############################################################
    # ## Elastic Net with self defined lambda (NULL not allowed) #############
    # ##############################################################
    print('Elnet')
    out_elnet = boot_stab(num_boot = 100, method = 'elnet', 
                          dat_file = '../../BMI/BMI_Lin_2014_binary.RData')

    save(out_elnet, file=paste0(dir, '/BMI_binary_elnet.RData', sep=''))
}else if (i == 4){
    # ##############################################################
    # ## Random Forest with Altman feature selection #############
    # ##############################################################
    print('RF')
    out_rf = boot_stab(num_boot = 100, method = 'RF', method.perm='altmann', 
                           dat_file = '../../BMI/BMI_Lin_2014_binary.RData')

    save(out_rf, file=paste0(dir, '/BMI_binary_rf.RData', sep=''))    
}else if (i == 5){
    ########################################################################
    ############### double bootstrap: random forest ##################
    # ######################################################################
    print('double boot RF')
    boot_rf = boot_stab_data(num_boot=100, method = 'RF',
                             data_file='../../BMI/BMI_Lin_2014_binary.RData')
    save(boot_rf, file=paste0(dir, '/BMI_binary_boot_rf.RData'))
}else if (i == 6){
    ########################################################################
    ############### double bootstrap: generalized compositional lasso ######
    # ######################################################################
    print('double boot compLasso')
    boot_compLasso = boot_stab_data(num_boot=100, method = 'GenCompLasso',
    	 							data_file='../../BMI/BMI_Lin_2014_binary.RData')
    save(boot_compLasso, file=paste0(dir, '/BMI_binary_boot_compLasso.RData'))
}    





