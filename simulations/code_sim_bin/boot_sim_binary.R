#-----------------------------------------------------------------------------------------------
# hypothesis testing based on boostrapped confidence interval in selected simulation scenarios
#-----------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

library(FSA)
source('../../code_method/cv_method_binary_update.R')
source('../../code_method/getStability.R')
source('../../code_method/bootstrap_test_compLasso_rf_binary.R')

M <- as.numeric(Sys.getenv("PBS_ARRAYID"))

if (M == 1){
    print('toeplitz: random forests')
    toe_rf = boot_stab_sim(num_boot=100, sim_file= paste0(dir, '/sim_toeplitz_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'RF', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)
    save(toe_rf, file=paste0(dir, '/binary_update/boot_toe_RF_binary.RData'))
}else if (M == 2){
    print('toeplitz: generalized compositional lasso')
    toe_genCompLasso = boot_stab_sim(num_boot=100, sim_file= paste0(dir, '/sim_toeplitz_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'GenCompLasso', seednum=31, ratio.training=0.8, fold.cv=10, 
				   lambda.coda=seq(0.1, 0.2, 0.01))
    save(toe_genCompLasso, file=paste0(dir, '/binary_update/boot_toe_genCompLasso_binary.RData'))
}else if (M == 3){
    print('block: random forests')
    block_rf = boot_stab_sim(num_boot=100, sim_file= paste0(dir, '/sim_block_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'RF', seednum=31, ratio.training=0.8, fold.cv=10, 
				   mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05)

    save(block_rf, file=paste0(dir, '/binary_update/boot_block_RF_binary.RData'))
}else if (M == 4){
    print('block: generalized compositional lasso')
    block_genCompLasso = boot_stab_sim(num_boot=100, sim_file= paste0(dir, '/sim_block_corr0.5P_1000_N_100.RData', sep=''), 
				   method= 'GenCompLasso', seednum=31, ratio.training=0.8, fold.cv=10, 
				   lambda.coda=seq(0.1, 0.2, 0.01))
    save(block_genCompLasso, file=paste0(dir, '/binary_update/boot_block_genCompLasso_binary.RData'))
}    






