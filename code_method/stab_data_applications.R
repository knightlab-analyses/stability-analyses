##########################################################
### estimate correlation between stability index #########
##########################################################

#source('cv_method.R') 
#source('getStability.R')

## set up parallel computing
library(foreach)
library(doParallel)
#numCores <- detectCores() - 2 # 6 cores
numCores <- detectCores() - 4 # for old computer
registerDoParallel(numCores)  # use multicore, set to the number of our cores

boot_stab = function(num_boot=100, dat_file, method, ratio.training=0.8, fold.cv=10, 
	                 family='gaussian', lambda.grid=exp(seq(-4, -2, 0.2)), alpha.grid=seq(0.1, 0.9, 0.1), 
	                 mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05, method.perm='altmann'){
	
	# load data
	load(dat_file)

	# bootsrap with parallelization
	results = foreach (i=1:num_boot) %dopar% {
		print(paste('bootnum', i, sep=":"))
		set.seed(i) # ensure each bootstrapped data is the same across methods
		N = length(y) # number of samples
		sample_ids = seq(1, N, 1)

	    # bootstrapped samples
	    boot_ids = sample(sample_ids, size=N, replace=TRUE)
	    boot_taxa = taxa[boot_ids, ]
	    boot_mf = y[boot_ids]

	    ## select features from lasso/elnet
	    if (method == 'compLasso'){
	        result.lin = cons_lasso_cv(y=boot_mf, datx=boot_taxa, seednum=i, ratio.training=ratio.training)
	        output.lin = c(result.lin$MSE, result.lin$coef.chosen) 

	    } else if (method == 'RF'){
	        result.rf = randomForest_cv(y=boot_mf, datx=boot_taxa, seednum=i, fold.cv=fold.cv, 
    							        num_trees=num_trees, mtry.grid = mtry.grid, pval_thr=pval_thr, method.perm=method.perm)
            output.rf = c(result.rf$MSE, result.rf$coef.chosen) 

	    } else if (method == 'lasso'){
			result.lasso = lasso_cv(y=boot_mf, datx=boot_taxa, seednum=i,family=family, lambda.choice='lambda.1se',
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    output.lasso = c(result.lasso$MSE, result.lasso$coef.chosen) 

	    } else if (method == 'elnet'){
		    result.elnet = elnet_cv(y=boot_mf, datx=boot_taxa, seednum=i,family=family, alpha.grid=alpha.grid,
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    output.elnet = c(result.elnet$MSE, result.elnet$coef.chosen) 
	    }
	}

	# reformat results
	p = dim(taxa)[2]
	stability_table = matrix(rep(0, num_boot * p), ncol=p)
	results_mse = results_chosen = list()
	for (b in 1:num_boot){
		results_mse[b] = results[[b]][1]
		results_chosen[[b]] = results[[b]][-1]
		stability_table[b, results_chosen[[b]]] = 1 
	}

	stab_index = round(getStability(stability_table)$stability, 2)
	MSE_mean = round(mean(unlist(results_mse), na.rm=T),2)
	MSE_se = round(FSA::se(unlist(results_mse), na.rm=T),2)

	results_list=list(num_boot=num_boot, method=method, stab_index=stab_index, stab_table=stability_table,
				 results_chosen=results_chosen, lists_mse=results_mse, MSE_mean=MSE_mean, MSE_se=MSE_se)

}
