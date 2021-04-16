#####################################################################################
### run all methods on simulated data with independent correlation ##################
#####################################################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]
#dir = '../sim_data'

library(FSA)
source('../../code_method/cv_method_binary_update.R')
source('../../code_method/getStability.R')
source('cv_sim_apply_binary_update.R')


dim.list = list()
size = c(50, 100, 500, 1000)
idx = 0
for (P in size){
    for (N in size){
        idx = idx + 1
        dim.list[[idx]] = c(P=P, N=N)
    }
}

files = NULL
for (dim in dim.list){
    p = dim[1]
    n = dim[2]
    files = cbind(files, paste0(dir, '/sim_independent_', paste('P', p, 'N', n, sep='_'), '.RData'))
}

#------------------------
# compare methods
#------------------------
i <- as.numeric(Sys.getenv("PBS_ARRAYID"))

# ###################
# #### Lasso ########
# ###################
print('Lasso')
file_name = files[[i]]
results_ind_lasso = sim_evaluate_cv(sim_file=file_name, method='lasso', family='binomial')

save(file=paste0(dir, '/binary_update/ind_Lasso_binary_', i, '.RData'), results_ind_lasso, file_name)

# # #########################
# # #### Elastic Net ########
# # #########################
print('Elnet')
file_name = files[[i]]
results_ind_elnet = sim_evaluate_cv(sim_file=file_name, method='elnet', family='binomial')

save(file=paste0(dir, '/binary_update/ind_Elnet_binary_', i, '.RData'), results_ind_elnet, file_name)

# ############################
# #### Random Forests ########
# ############################
print('Random Forests')
file_name = files[[i]]
results_ind_rf = sim_evaluate_cv(sim_file=file_name, method='RF')
save(file=paste0(dir, '/binary_update/ind_RF_binary_', i, '.RData'), results_ind_rf, file_name)

# #########################################################
# #### Generalized compositional Lasso ########
# #########################################################
print('Generalized Compositional Lasso')
file_name = files[[i]]
results_ind_GenCompLasso = sim_evaluate_cv (sim_file=file_name, method='GenCompLasso')

save(file=paste0(dir, '/binary_update/ind_GenCompLasso_binary_', i, '.RData'),
     results_ind_GenCompLasso, file_name)

