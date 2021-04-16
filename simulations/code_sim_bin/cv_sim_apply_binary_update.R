#######################################################################################
### apply different feature selection methods to simulated data #######################
#######################################################################################

# library(FSA) # for se()
# source('cv_method_binary_update.R') 
# source('getStability.R')

sim_evaluate_cv = function(sim_file, method, seednum=31, ratio.training=0.8,
                           fold.cv=10, family='binomial',
                           lambda.grid=exp(seq(-4, -2, 0.2)), 
                           alpha.grid=seq(0.1, 0.9, 0.1), 
	                       mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05,
                           method.perm='altmann', 
                           lambda.coda=seq(0.1, 0.2, 0.01), data.split=FALSE){
	# load simulated data
	load(sim_file, dat <- new.env())

	idx.start = 1; idx.stop = 100 # 100 repetitions for each simulated scenario
	rou = dat$sim_array[[1]]$rou # rou, n, p are same across all repetitions
	n = dat$sim_array[[1]]$n
	p = dat$sim_array[[1]]$p

	# evaluating different methods	
	fp = fn = roc = NULL
	stability.table = matrix(rep(0, (idx.stop - idx.start + 1) * p), ncol=p)
	colnames(stability.table) = paste('V', seq(1:p), sep='')

	for (i in idx.start:idx.stop){ 
        print(paste0('simulation idx: ', i))
	    sub = dat$sim_array[[i]]
	    coef = sub$beta
	    coef.true = which(coef != 0)
	    y_binary = as.factor(ifelse(sub$Y >= median(sub$Y), 1, 0))
	    
		if (method == 'lasso'){
			result = lasso_cv(y=y_binary, datx=sub$Z, seednum=seednum,family=family, lambda.choice='lambda.1se',
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    select.features = result$coef.chosen 
		} else if (method == 'elnet'){
		    result = elnet_cv(y=y_binary, datx=sub$Z, seednum=seednum,family=family, alpha.grid=alpha.grid,
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    select.features = result$coef.chosen 
		} else if (method == 'RF'){
		    result = randomForest_cv(y=y_binary, datx=sub$Z, seednum=seednum, fold.cv=fold.cv, 
		                             num_trees=num_trees, mtry.grid = mtry.grid, pval_thr=pval_thr, method.perm=method.perm)
		    select.features = result$coef.chosen 
		} else if (method == 'GenCompLasso'){
            # note that generalized compositonal lasso do X to Z transformation within
			result = gen_cons_lasso_cv (y=y_binary, datx=sub$X, seednum=seednum,
                                        data.split=data.split,
                                        ratio.training=ratio.training,
                                        lambda.coda=lambda.coda)
			select.features = result$coef.chosen
		}

	    # false positives: shouldn't be chosen but chosen
	    fp = c(fp, length(setdiff(select.features, coef.true)))

	    # false negatives: # should be chosen yet not
	    fn = c(fn, length(setdiff(coef.true, select.features)))

	    # classification error
	    roc = c(roc, result$ROC)

	    # stability table
	    stability.table[i, select.features] = 1
	}

	# store results
	FP = paste(mean(fp), '(', round(se(fp),2), ')')
	FN = paste(mean(fn), '(', round(se(fn),2), ')')
	ROC = paste(round(mean(roc, na.rm=T),2), '(', round(se(roc, na.rm=T),2), ')')
	Stab = round(getStability(stability.table)$stability, 2)

	results=list(rou=rou, n=n, p=p, FP=FP, FN=FN, ROC=ROC, Stab=Stab, 
		         Stab.table=stability.table, FP.list=fp, FN.list=fn, ROC.list=roc)

}



