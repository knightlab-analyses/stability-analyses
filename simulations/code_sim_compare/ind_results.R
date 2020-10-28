#####################################################################################
### run all methods on simulated data with independent correlation ##################
#####################################################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
#dir = args[1]
dir = '../sim_data'

source('cv_method.R')
source('getStability.R')
source('cv_sim_apply.R')


library(FSA)
library(foreach)
library(doParallel)
numCores <- detectCores() - 2 
registerDoParallel(numCores)  

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

##################
### Lasso ########
##################
print('Lasso')
results_ind_lasso = foreach(i = iter(files)) %dopar%{
   sim_evaluate_cv(sim_file=i, method='lasso')
}

save(file=paste0(dir, '/independent_Lasso.RData'), results_ind_lasso)

#########################################################
#### Lasso with default lambda sequence ########
#########################################################
print('Lasso')
results_ind_lasso = foreach(i = iter(files)) %dopar%{
   sim_evaluate_cv(sim_file=i, method='lasso', lambda.grid=NULL)
}

save(file=paste0(dir, '/independent_Lasso_lambdaNULL.RData'), results_ind_lasso)


#########################
#### Elastic Net ########
#########################
print('Elnet')
results_ind_elnet= foreach(i = iter(files)) %dopar%{
   print(i)
   sim_evaluate_cv(sim_file=i, method='elnet')
}

save(file=paste0(dir, '/independent_Elnet.RData'), results_ind_elnet)

########################## #########################
#### Elastic Net with default lambda sequence ########
########################## #########################
this cannot be done, as Elastic Net code requires lambda
print('Elnet')
results_ind_elnet= foreach(i = iter(files)) %dopar%{
   print(i)
   sim_evaluate_cv(sim_file=i, method='elnet', lambda.grid=NULL)
}

save(file=paste0(dir, '/independent_Elnet_lambdaNULL.RData'), results_ind_elnet)

############################
#### Random Forests ########
############################
print('RF')
results_ind_rf = foreach(i = iter(files)) %dopar%{
   print(i)
   sim_evaluate_cv(sim_file=i, method='RF')
}

save(file=paste0(dir, '/independent_RF.RData'), results_ind_rf)

############################################################
#### Random Forests with Janitza permutation method ########
############################################################
print('RF')
results_ind_rf = foreach(i = iter(files)) %dopar%{
   print(i)
   sim_evaluate_cv(sim_file=i, method='RF', method.perm='janitza')
}

save(file=paste0(dir, '/independent_RF_janitza.RData'), results_ind_rf)
