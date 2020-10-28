##########################################################
### estimate correlation between stability index #########
##########################################################

source('cv_method.R') 
source('getStability.R')
source('cv_sim_apply.R')

## set up parallel computing
library(foreach)
library(doParallel)
numCores <- detectCores() - 2 # 6 cores
registerDoParallel(numCores)  # use multicore, set to the number of our cores

boot_stab_sim = function(num_boot=100, sim_file, method, seednum=31, ratio.training=0.8, fold.cv=10, 
					  family='gaussian', lambda.grid=exp(seq(-4, -2, 0.2)), alpha.grid=seq(0.1, 0.9, 0.1), 
	                  mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05){
	
	# load simulated data
	load(sim_file, dat <- new.env())

	idx.start = 1; idx.stop = 100 # 100 repetitions for each simulated scenario
	rou = dat$sim_array[[1]]$rou # rou, n, p are same across all repetitions
	n = dat$sim_array[[1]]$n
	p = dat$sim_array[[1]]$p

	## get a vector of stability index from bootstrapped data
	stab_index = rep(0, num_boot)
	b = 0
	for (i in idx.start:idx.stop){
		b = b + 1
		print(paste('index', i, sep=':'))
	    sub = dat$sim_array[[i]]

	    # bootsrap with parallelization
		selections = foreach (i=1:num_boot) %dopar% {
			N = length(sub$Y) # number of samples
			boot_ids = sample(N, size=N, replace=TRUE)
			boot_Z = sub$Z[boot_ids, ]
			boot_Y = sub$Y[boot_ids, ]

		    ## select features from lasso/elnet
		    if (method == 'compLasso'){
		        result.lin = cons_lasso_cv(y=boot_Y, datx=boot_Z, seednum=i, ratio.training=ratio.training)
		        select.lin = result.lin$coef.chosen  # since 1 represents intercept

		    } else if (method == 'RF'){
		        result.rf = randomForest_cv(y=boot_Y, datx=boot_Z, seednum=i, fold.cv=fold.cv, 
	    							        num_trees=num_trees, mtry.grid = mtry.grid, pval_thr=pval_thr)
	            select.rf = result.rf$coef.chosen 
		    }
		}

		# calculate stability index from bootstrapped data
		stability_table = matrix(rep(0, num_boot * p), ncol=p)
		for (j in 1:num_boot){
		stability_table[j, selections[[j]]] = 1 
		}

		stab_index[b] = round(getStability(stability_table)$stability, 2)
	}	

	results=list(rou=rou, n=n, p=p, num_boot=num_boot, method=method, stab_index=stab_index)

}


########################################################
## double bootstrap applied to real data application ###
########################################################
boot_stab_data = function(num_boot=100, data_file, method, seednum=31, ratio.training=0.8, fold.cv=10, 
					      family='gaussian', lambda.grid=exp(seq(-4, -2, 0.2)), alpha.grid=seq(0.1, 0.9, 0.1), 
	                      mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05){
	
	# load simulated data
	set.seed(seednum)
	load(data_file)
	p = dim(taxa)[2]

	stab_index = rep(0, num_boot)
	MSE_list = list()
	# first loop of bootstrap to generate num_boot bootstrapped datasets
	for (k in 1:num_boot){
        print(paste('num_boot', k, sep=':'))
		N = length(y) # number of samples
		sample_ids = seq(1, N, 1)

		# bootstrapped samples
	    boot_ids = sample(sample_ids, size=N, replace=TRUE)
	    boot_taxa = taxa[boot_ids, ]
	    boot_mf = y[boot_ids]

	    ## second loop of bootstrap to perform variable selection
		results = foreach (i=1:num_boot) %dopar% {
			boot_ids_second = sample(N, size=N, replace=TRUE)
			boot_Z = boot_taxa[boot_ids_second, ]
			boot_Y = boot_mf[boot_ids_second]

		    ## select features from lasso/elnet
		    if (method == 'compLasso'){
		        result.lin = cons_lasso_cv(y=boot_Y, datx=boot_Z, seednum=i, ratio.training=ratio.training)
		        output.lin = c(result.lin$MSE, result.lin$coef.chosen) 

		    } else if (method == 'RF'){
		        result.rf = randomForest_cv(y=boot_Y, datx=boot_Z, seednum=i, fold.cv=fold.cv, 
	    							        num_trees=num_trees, mtry.grid = mtry.grid, pval_thr=pval_thr)
	            output.rf = c(result.rf$MSE, result.rf$coef.chosen) 
		    }
		}

		# reformat results (stability & MSE)
		stability_table = matrix(rep(0, num_boot * p), ncol=p)
		results_mse = results_chosen = list()
		for (b in 1:num_boot){
			results_mse[b] = results[[b]][1]
			results_chosen[[b]] = results[[b]][-1]
			stability_table[b, results_chosen[[b]]] = 1 
		}

		stab_index[k] = round(getStability(stability_table)$stability, 2)
		MSE_list[[k]] = results_mse
	}


	results=list(num_boot=num_boot, method=method, stab_index=stab_index, MSE_list=MSE_list)

}
