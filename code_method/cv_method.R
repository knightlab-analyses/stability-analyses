###########################################
#### methods for comparisions #############
###########################################
library(glmnet)
library(caret)
library(ranger) # faster random forest

# general reference with caret on glmnet & random forest
## http://rstudio-pubs-static.s3.amazonaws.com/251240_12a8ecea8e144fada41120ddcf52b116.html

# important note
## when training model on training set, as tuning parameters are set in final model, coefficient are chosen already
## use test set to predict the MSE
## as unclear how carret handle fitting and prediction, use carret only for parameter tuning in elastic net

###########################################
#### Lasso (tune lambda) ##################
###########################################
#' @title lasso_cv
#' @description Does a K-fold cross-validation for Lasso.
#' @param datx The input data matrix.
#' @param y The response variable.
#' @param seednum The seed number (default=31).
#' @param family The family of linear regression models (default="gaussian").
#' @param ratio.training The ratio of the whole data assigned for model training (default=0.8).
#' @param fold.cv The number of folds for cross-validation.(default=10)
#' @param lambda.grid The tuning range for the regularization parameter "lambda" of "lasso".
#' @return A list of Lasso model output including MSE, and selected features.

lasso_cv = function(datx, y, seednum=31, family=family, ratio.training=0.8, fold.cv=10, 
					lambda.grid, lambda.choice='lambda.1se'){
	# seednum: the seed number
	# ratio.training: the ratio of training set (parento principle training:test=8:2)
	# lambda.grid: possible candidate values for tuning parameter Lambda 
	# fold.cv: n-fold cross validation
	# lambda.choice: 'lambda.min' or 'lambda.1se' 

	set.seed(seednum)

	# # split data into training and tests sets
	nn <- length(y)
	trn <- sample(1:nn, ratio.training*nn)
	x.train <- datx[trn, ]
	x.test <- datx[-trn, ]
	y.train <- y[trn]
	y.test <- y[-trn]

	# use cross-validation on training data & fit on test data
	cv.fit <- cv.glmnet(x.train, y.train, family=family, alpha=1, lambda=lambda.grid, nfolds=fold.cv)
	pred.fit <- predict(cv.fit, s=lambda.choice, newx=x.test, type='response')

	# covariates chosen 
	coef.fit <- predict(cv.fit, s=lambda.choice, newx=x.test, type='coefficients')
	coef.chosen = which(coef.fit != 0) 
    coef.chosen = coef.chosen - 1 # index for 1 representing intercept

	# evaluate results
	MSE <- mean((y.test - pred.fit)^2)


	result = list(MSE=MSE, coef.chosen=coef.chosen)
	return(result)
}

###########################################
#### Elastic Net (tune lambda and alpha) ##
###########################################
# reference 1 (tune both parameters with caret)
##  https://stats.stackexchange.com/questions/268885/tune-alpha-and-lambda-parameters-of-elastic-nets-in-an-optimal-way
# reference 2 (extract final model and prediction with caret)
## https://topepo.github.io/caret/model-training-and-tuning.html

#' @title elnet_cv
#' @description Does a K-fold cross-validation for elastic-net.
#' @param datx The input data matrix.
#' @param y The response variable.
#' @param seednum The seed number (default=31).
#' @param family The family of linear regression models (default="gaussian").
#' @param ratio.training The ratio of the whole data assigned for model training (default=0.8).
#' @param fold.cv The number of folds for cross-validation.(default=10)
#' @param lambda.grid The tuning range for the regularization parameter "lambda" of "lasso".
#' @param alpha.grid The tuning range for the elastic-net mixing parameter "alpha".
#' @return A list of elnet model output including MSE, and selected features.

elnet_cv = function(datx, y, seednum=31, alpha.grid, lambda.grid, family=family, 
	                ratio.training=0.8, fold.cv=10){
	# seednum: the seed number
	# ratio.training: the ratio of training set (parento principle training:test=8:2)
	# alpha.grid: possible candidate values for tuning parameter alpha
	# lambda.grid: possible candidate values for tuning parameter Lambda 
	# fold.cv: n-fold cross validation

	set.seed(seednum)

	# # split data into training and tests sets
	nn <- length(y)
	trn <- sample(1:nn, ratio.training*nn)
	x.train <- datx[trn, ]
	x.test <- datx[-trn, ]
	y.train <- y[trn]
	y.test <- y[-trn]

	data.train = as.data.frame(cbind(y.train, x.train)) 
	colnames(data.train) = c('y', paste('V', seq(1, dim(datx)[2]), sep=''))

	# tune parameters with CV on training data with caret
	trnCtrl <- trainControl(method = "cv", number = fold.cv)
	srchGrid <- expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)

	my_train <- train(y ~., data.train,
              		  method = "glmnet",
                      tuneGrid = srchGrid,
                      trControl = trnCtrl)

	# fit the model with best tuning parameters on test data with glmnet
	fit = glmnet(x.train, y.train, alpha=my_train$bestTune$alpha, lambda=my_train$bestTune$lambda)
	pred.fit <- predict(fit, x.test, type='response')
	
	# covariates chosen 
	coef.fit <- predict(fit, x.test, type='coefficients')
	coef.chosen = which(coef.fit != 0) 
    coef.chosen = coef.chosen - 1 # index for 1 representing intercept

	# evaluate results
	MSE <- mean((y.test - pred.fit)^2)


	result = list(MSE=MSE, coef.chosen=coef.chosen)
	return(result)
}


###############################################################
#### Random Forests (tune # variables at each random split ####
###############################################################
# reference on regression random forst
## https://uc-r.github.io/random_forests
# OOB error is different from test error (see above website)

#' @title randomForest_cv
#' @description Does a K-fold cross-validation for random forests.
#' @param datx The input data matrix.
#' @param y The response variable.
#' @param sim_file The path of a .RData file includes the simulation data pre-generated.
#' @param seednum The seed number (default=31).
#' @param fold.cv The number of folds for cross-validation.(default=10)
#' @param ratio.training The ratio of the whole data assigned for model training (default=0.8).
#' @param mtry.grid The tuning range for the hyperparameter "mtry", that is, number of variables to possibly split at in each node. 
#' @param num_trees The number of trees grow in random forests.
#' @param pval_thr The threshold for the estimated p value of RF importance scores.
#' @param method.perm The permutation method for estimating the p value of RF importance scores.
#' @return A list of RF model output including the "mtry" grid, selected features, MSE, OOB, and p-values of selected features.

randomForest_cv = function(datx, y, seednum=31, fold.cv=5, ratio.training=0.8, mtry.grid=10, num_trees=500, 
						   pval_thr=0.05, method.perm='altmann'){
	# mtry: number of variables to randomly sample at each split
	# num_trees: number of trees to grow in random forests
	# pval_thr: threshold for permutation test
	# note that permutation p-value can use "altmann method" for all types of data; 'Janita' for high-dimensitional data only
	# ref on permutation methods: http://finzi.psych.upenn.edu/library/ranger/html/importance_pvalues.html

	set.seed(seednum)  

	# split data into training and tests sets
	data = as.data.frame(cbind(y, datx)) 
	colnames(data) = c('y', paste('V', seq(1, dim(datx)[2]), sep=''))
	inTraining = createDataPartition(data$y, p = ratio.training, list=FALSE)
	train <- data[inTraining, ]
	test <- data[-inTraining, ]  

	# tune parameter with cross validation
	hyper.grid <- expand.grid(mtry = mtry.grid, OOB_RMSE = 0)
	for (i in 1:nrow(hyper.grid)){
	  model = ranger::ranger(y ~., data = train,
	                 num.trees=500, mtry=hyper.grid$mtry[i],
	                 seed=seednum, importance = 'permutation')
	  hyper.grid$OOB_RMSE[i] = sqrt(model$prediction.error)
	}
	OOB = min(hyper.grid$OOB_RMSE) # out of bag error
	position = which.min(hyper.grid$OOB_RMSE)

	# permutation test on tuned random forst model to obtain chosen features
	if (method.perm == 'altmann'){ # for all data types
		rf.model <- ranger::ranger(y ~., data=test, num.trees = num_trees,
					       mtry = hyper.grid$mtry[position], importance = 'permutation')
	    table = as.data.frame(importance_pvalues(rf.model, method = "altmann", 
                              formula = y ~ ., data = test)) 
	} else if (method.perm == 'janitza'){ # for high dimensional data only
		rf.model <- ranger::ranger(y ~., data=test, num.trees = num_trees,
					       mtry = hyper.grid$mtry[position], importance = 'impurity_corrected')
	    table = as.data.frame(importance_pvalues(rf.model, method = "janitza", 
                              formula = y ~ ., data = test)) 		
	}

    coef.chosen = which(table$pvalue < pval_thr) 

    # if nothing been selected
    if (identical(coef.chosen, integer(0))){
    	coef.chosen = 0
    }

    # obtain additional prediction error to make comparable to other methods
    pred_rf = predict(rf.model, test)
    MSE <- mean((test$y - pred_rf$predictions)^2)

	result = list(mtry=mtry.grid[position], coef.chosen=coef.chosen, MSE=MSE, OOB=OOB, p.value=table$pvalue)
	return(result)

}


###############################################################
#### Compositional Lasso by Lin et al 2014 ####################
###############################################################
dyn.load("../../code_Lin/cvs/cdmm.dll")
source("../../code_Lin/cvs/cdmm.R")

#' @title cons_lasso_cv
#' @description Does a K-fold cross-validation for elastic-net.
#' @param datx The input data matrix.
#' @param y The response variable.
#' @param seednum The seed number (default=31).
#' @param ratio.training The ratio of the whole data assigned for model training (default=0.8).
#' @return A list of compLasso model output including MSE, and selected features.

cons_lasso_cv = function(datx, y, seednum, ratio.training=0.8){
	set.seed(seednum)
	z = datx
	n = length(y)

	itrn = sample(n, ratio.training*n)
	itst = setdiff(1:n, itrn)
	ans <- cv.cdmm(y[itrn], z[itrn, ], refit=TRUE) # proposed method (default: 10 fold CV)
	bet <- ans$bet; int <- ans$int
	pe <- mean((y[itst] - int - z[itst, ] %*% bet)^2)

	coef.chosen = which(bet != 0) # the first beta refer to 1st feature
	MSE = pe

	result = list(MSE=MSE, coef.chosen=coef.chosen)
	return(result)
}

###########################################################
#### Lasso (double cv -- no validation set) ###############
###########################################################
lasso_double_cv = function(datx, y, seednum=31, family=family, fold.cv=10, lambda.grid){
	set.seed(seednum)

	## double loop for cross-vlidation 
	flds <- caret::createFolds(y, k = fold.cv, list = TRUE, returnTrain = FALSE)

	## outer loop for estimating MSE
	MSE = STAB = matrix(rep(0, fold.cv * length(lambda.grid)), nrow=fold.cv)
	rownames(MSE) = rownames(STAB) = paste('fold', seq(1:fold.cv), sep='')
	colnames(MSE) = colnames(STAB) = paste('lambda', seq(1:length(lambda.grid)), sep='')
	for (b in 1:fold.cv){
		idx.test = flds[[b]]
		x.train = datx[-idx.test, ]
		y.train = y[-idx.test]
		x.test = datx[idx.test, ]
		y.test = y[idx.test]

		fit = glmnet(x.train, y.train, family=family, alpha=1, lambda=lambda.grid)
		pred.fit <- predict(fit, newx=x.test, s=lambda.grid, type='response')
		MSE[b, ] <- colMeans((y.test - pred.fit)^2)

		## inner loup for estimating stability
		flds.inn <- createFolds(y.train, k = fold.cv, list = TRUE, returnTrain = FALSE)
		stab.table.inn = list()
		for (bb in 1:fold.cv){
		idx.test.inn = flds.inn[[bb]]
		stab.x.train = x.train[-idx.test.inn, ]
		stab.y.train = y.train[-idx.test.inn]
		stab.x.test = x.train[idx.test.inn, ]
		stab.y.test = y.train[idx.test.inn]

		fit.inn = glmnet(stab.x.train, stab.y.train, family=family, alpha=1, lambda=lambda.grid)
		coef.fit.inn <- predict(fit.inn, newx=stab.x.test, s=lambda.grid, type='coefficients')
		coef.chosen.table = as.matrix(coef.fit.inn)
		coef.chosen.table = coef.chosen.table[-1, ] # no need of intercept
		coef.chosen.table[coef.chosen.table != 0] = 1 # convert to binary table
		coef.chosen.table[coef.chosen.table == 0] = 0
		colnames(coef.chosen.table) = paste('lambda', seq(1:length(lambda.grid)), sep='')
		stab.table.inn[[bb]] = coef.chosen.table
		}

		# calculate stability
		table.tmp = matrix(rep(0, dim(datx)[2] * fold.cv), nrow=fold.cv)
		for (i in 1:length(lambda.grid)){
			for (bb in 1:fold.cv){
				table.tmp[bb, ] = stab.table.inn[[bb]][, i]
			}
			STAB[b, i] = round(getStability(table.tmp)$stability, 2)
		}
	}

	result = list(lambda.grid=lambda.grid, MSE.list=MSE, STAB.list=STAB, MSE.value=colMeans(MSE), STAB.value=colMeans(STAB))
	return(result)
}


###########################################################
#### RF (double cv -- no validation set) ###############
###########################################################
randomForest_double_cv = function(datx, y, seednum=31, fold.cv=5, mtry.grid=10, num_trees=500, pval_thr=0.05){
	set.seed(seednum)  

	## double loop for cross-vlidation 
	data = as.data.frame(cbind(y, datx)) 
	colnames(data) = c('y', paste('V', seq(1, dim(datx)[2]), sep='')) 
	flds <- createFolds(data$y, k = fold.cv, list = TRUE, returnTrain = FALSE)

	## outer loop for estimating MSE
	MSE = STAB = matrix(rep(0, fold.cv * length(mtry.grid)), nrow=fold.cv)
	rownames(MSE) = rownames(STAB) = paste('fold', seq(1:fold.cv), sep='')
	colnames(MSE) = colnames(STAB) = paste('mtry', seq(1:length(mtry.grid)), sep='')
	for (b in 1:fold.cv){
		idx.test = flds[[b]]
		train <- data[-idx.test, ]
		test <- data[idx.test, ] 		

		for (mtry in mtry.grid){
			fit = ranger(y ~., data = train, num.trees=500, mtry=mtry, seed=seednum, importance = 'permutation')
			pred.fit <- predict(fit, test)
			MSE[b, ] <- mean((test$y - pred.fit$predictions)^2)
		}
		
		## inner loup for estimating stability
		flds.inn <- createFolds(train$y, k = fold.cv, list = TRUE, returnTrain = FALSE)
		stab.table.inn = list()
		for (bb in 1:fold.cv){
		idx.test.inn = flds.inn[[bb]]
		stab.train <- train[-idx.test.inn, ]
		stab.test <- train[idx.test.inn, ]

			coef.chosen.table = matrix(rep(0, dim(datx)[2] * length(mtry.grid)), ncol=length(mtry.grid))
			colnames(coef.chosen.table) = paste('mtry', seq(1:length(mtry.grid)), sep='')
			idx = 0
			for (mtry in mtry.grid){
				idx = idx + 1
				fit.inn = ranger(y ~., data = stab.train, num.trees=500, mtry=mtry, seed=seednum, importance = 'permutation')
				table = as.data.frame(importance_pvalues(fit.inn, method = "altmann", formula = y ~ ., data = stab.train)) 
				coef.chosen = which(table$pvalue < pval_thr)
				coef.chosen.table[coef.chosen, idx] = 1
				stab.table.inn[[bb]] = coef.chosen.table
			}

		}

		# calculate stability
		table.tmp = matrix(rep(0, dim(datx)[2] * fold.cv), nrow=fold.cv)
		for (i in 1:length(mtry.grid)){
			for (bb in 1:fold.cv){
				table.tmp[bb, ] = stab.table.inn[[bb]][, i]
			}
			STAB[b, i] = round(getStability(table.tmp)$stability, 2)
		}
	}

	result = list(mtry.grid=mtry.grid, MSE.list=MSE, STAB.list=STAB, MSE.value=colMeans(MSE), STAB.value=colMeans(STAB))
	return(result)

}
