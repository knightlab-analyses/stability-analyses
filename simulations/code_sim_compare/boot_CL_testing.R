#########################################################################################################
### This is to use bootstrap for testing on compositional lasso #################
#########################################################################################################

#args = commandArgs(trailingOnly=TRUE)
#print(args)
#dir = args[1]

source('cv_method.R')
source('getStability.R')
source('bootstrap_test_compLasso_rf.R')

dir = '../sim_data'

toe_lin = boot_stab(num_boot=100, sim_file= paste0(dir, '/sim_toeplitz_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'compLasso', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)

save(toe_lin, file=paste0(dir, '/boot_toe_compLasso.RData'))

block_lin = boot_stab(num_boot=100, sim_file= paste0(dir, '/sim_block_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'compLasso', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)

save(block_lin, file=paste0(dir, '/boot_block_compLasso.RData'))






