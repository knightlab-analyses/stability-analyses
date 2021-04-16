#########################################################################################################
### This is to use bootstrap for testing on random forests #################
#########################################################################################################

args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../../code_method/cv_method.R')
source('../../code_method/getStability.R')
source('../../code_method/bootstrap_test_compLasso_rf.R')

# the function boot_stab is probably upddated to be boot_stab_sim() now

toe_rf = boot_stab(num_boot=100, sim_file= paste0(dir, '/sim_toeplitz_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'RF', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)

save(toe_rf, file=paste0(dir, '/boot_toe_RF.RData'))

block_rf = boot_stab(num_boot=100, sim_file= paste0(dir, '/sim_block_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'RF', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)

save(block_rf, file=paste0(dir, '/boot_block_RF.RData'))







