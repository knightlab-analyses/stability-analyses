#####################################################################
#### Compositional Lasso on simulation results ######################
#####################################################################
library(FSA) # for se()
source('cv_method.R') 
source('getStability.R')
source('cv_sim_apply.R')

dir = '../sim_data'
dim.list = list()
size = c(50, 100, 500, 1000)
idx = 0
for (P in size){
    for (N in size){
        idx = idx + 1
        dim.list[[idx]] = c(P=P, N=N)
    }
}

###########################################
#### Independent simulations ##############
###########################################
files = NULL
for (dim in dim.list){
    p = dim[1]
    n = dim[2]
    files = cbind(files, paste0(dir, '/sim_independent_', paste('P', p, 'N', n, sep='_'), '.RData'))
}

results_ind_compLasso = NULL
for (i in 1:length(files)){ # parallel computing not working
	print(i)
	results_ind_compLasso[[i]] = sim_evaluate_cv(sim_file=files[i], method='compLasso')
}
save(file=paste0(dir, '/independent_compLasso.RData'), results_ind_compLasso)

###########################################
#### Toeplitz simulations ##############
###########################################
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

results_toe_compLasso = NULL
for (i in 1:length(files)){ # parallel computing not working
	print(i)
	results_toe_compLasso[[i]] = sim_evaluate_cv(sim_file=files[i], method='compLasso')
}
save(file=paste0(dir, '/toe_compLasso.RData'), results_toe_compLasso)

###########################################
#### Block simulations ####################
###########################################
rou.list = seq(0.1, 0.9, 0.2)

files = NULL
for (rou in rou.list){
  for (dim in dim.list){
    p = dim[1]
    n = dim[2]
    files = cbind(files, paste0(dir, '/sim_block_corr', rou, paste('P', p, 'N', n, sep='_'), '.RData', sep=''))
  }
}

results_block_compLasso = NULL
for (i in 1:length(files)){ # parallel computing not working
	print(i)
	results_block_compLasso[[i]] = sim_evaluate_cv(sim_file=files[i], method='compLasso')
}
save(file=paste0(dir, '/block_compLasso.RData'), results_block_compLasso)









