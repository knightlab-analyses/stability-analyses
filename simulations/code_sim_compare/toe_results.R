#####################################################################################
### run all methods on simulated data with Toeplitz correlation ##################
#####################################################################################
# args = commandArgs(trailingOnly=TRUE)
# print(args)
# dir = args[1]
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

## correlation strength
rou.list = seq(0.1, 0.9, 0.2)

files = NULL
for (rou in rou.list){
  for (dim in dim.list){
    p = dim[1]
    n = dim[2]
    files = cbind(files, paste0(dir, '/sim_toeplitz_corr', rou, paste('P', p, 'N', n, sep='_'), '.RData', sep=''))
  }
}


##################
### Lasso ########
##################
print('Lasso')
results_toe_lasso = foreach(i = iter(files)) %dopar%{
   print(i)
   sim_evaluate_cv(sim_file=i, method='lasso')
}

save(file=paste0(dir, '/toe_Lasso.RData'), results_toe_lasso)


#########################
#### Elastic Net ########
#########################
print('Elnet')
results_toe_elnet= foreach(i = iter(files)) %dopar%{
   sim_evaluate_cv(sim_file=i, method='elnet')
}

save(file=paste0(dir, '/toe_Elnet.RData'), results_toe_elnet)

############################
#### Random Forests ########
############################
print('Random Forests')
results_toe_rf = foreach(i = iter(files)) %dopar%{
   sim_evaluate_cv(sim_file=i, method='RF')
}

save(file=paste0(dir, '/toe_RF.RData'), results_toe_rf)

