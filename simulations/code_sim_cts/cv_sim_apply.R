#######################################################################################
### apply different feature selection methods to simulated data #######################
#######################################################################################

#' library(FSA) # for se()
#' source('cv_method.R') 
#' source('getStability.R')
#' @title sim_evaluate_cv
#' @description 
#' @param sim_file The path of a .RData file includes the simulation data pre-generated.
#' @param method The machine learning algorithms to be used. This should be one of "lasso", "elnet", "RF", and "compLosso". 
#' @param seednum The seed number (default=31).
#' @param ratio.training The ratio of the whole data assigned for model training (default=0.8).
#' @param fold.cv The number of folds for cross-validation.(default=10)
#' @param family The family of linear regression models (only related to "lasso" or "elnet").
#' @param lambda.grid The tuning range for the regularization parameter "lambda" of "lasso" or "elnet".
#' @param alpha.grid The tuning range for the elastic-net mixing parameter "alpha".
#' @param mtry.grid The tuning range for the "RF" hyperparameter "mtry", that is, number of variables to possibly split at in each node. 
#' @param num_trees The number of decision trees built in the RF.
#' @param pval_thr The threshold for the estimated p value of RF importance scores.
#' @param method.perm The permutation method for estimating the p value of RF importance scores.
#' @details 

sim_evaluate_cv = function(sim_file, method, seednum=31, ratio.training=0.8, fold.cv=10, family='gaussian',
	                       lambda.grid=exp(seq(-4, -2, 0.2)), alpha.grid=seq(0.1, 0.9, 0.1), 
	                       mtry.grid=seq(5, 25, 5), num_trees = 500, pval_thr = 0.05, method.perm='altmann'){
	# load simulated data
	load(sim_file, dat <- new.env())

	idx.start = 1; idx.stop = 100 # 100 repetitions for each simulated scenario
	rou = dat$sim_array[[1]]$rou # rou, n, p are same across all repetitions
	n = dat$sim_array[[1]]$n
	p = dat$sim_array[[1]]$p

	# evaluating different methods	
	fp = fn = mse = OOB_rf = NULL
	stability.table = matrix(rep(0, (idx.stop - idx.start + 1) * p), ncol=p)
	colnames(stability.table) = paste('V', seq(1:p), sep='')

	for (i in idx.start:idx.stop){ 
	    sub = dat$sim_array[[i]]
	    coef = sub$beta
	    coef.true = which(coef != 0)
	    
		if (method == 'lasso'){
			result = lasso_cv(y=sub$Y, datx=sub$Z, seednum=seednum,family=family, lambda.choice='lambda.1se',
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    select.features = result$coef.chosen 
		} else if (method == 'elnet'){
		    result = elnet_cv(y=sub$Y, datx=sub$Z, seednum=seednum,family=family, alpha.grid=alpha.grid,
		                            ratio.training=ratio.training, fold.cv=fold.cv, lambda.grid=lambda.grid)
		    select.features = result$coef.chosen 
		} else if (method == 'RF'){
		    result = randomForest_cv(y=sub$Y, datx=sub$Z, seednum=seednum, fold.cv=fold.cv, 
		                             num_trees=num_trees, mtry.grid = mtry.grid, pval_thr=pval_thr, method.perm=method.perm)
		    select.features = result$coef.chosen 
		} else if (method == 'compLasso'){
			result = cons_lasso_cv(datx=sub$Z, y=sub$Y, seednum=seednum, ratio.training=ratio.training)
			select.features = result$coef.chosen
		}

	    # false positives: shouldn't be chosen but chosen
	    fp = c(fp, length(setdiff(select.features, coef.true)))

	    # false negatives: # should be chosen yet not
	    fn = c(fn, length(setdiff(coef.true, select.features)))

	    # prediction error
	    mse = c(mse, result$MSE)

	    if (method == 'RF'){
	    	OOB_rf = c(OOB_rf, result$OOB)
	    }

	    # stability table
	    stability.table[i, select.features] = 1
	}

	# store results
	FP = paste(mean(fp), '(', round(se(fp),2), ')') # FPR?
	FN = paste(mean(fn), '(', round(se(fn),2), ')') # FNR?
	MSE = paste(round(mean(mse, na.rm=T),2), '(', round(se(mse, na.rm=T),2), ')')
	Stab = round(getStability(stability.table)$stability, 2)

	results=list(rou=rou, n=n, p=p, FP=FP, FN=FN, MSE=MSE, Stab=Stab, 
		         Stab.table=stability.table, FP.list=fp, FN.list=fn, MSE.list=mse, OOB.list=OOB_rf)

}



