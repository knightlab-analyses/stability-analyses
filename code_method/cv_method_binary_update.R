###########################################
#### methods for binary outcome #############
###########################################
library(glmnet)
library(caret)
library(ranger) # faster random forest
library(pROC) # calculate ROC (need to load this library to avoid error)

# general reference with caret on glmnet & random forest
## http://rstudio-pubs-static.s3.amazonaws.com/251240_12a8ecea8e144fada41120ddcf52b116.html

# important note
## when training model on training set, as tuning parameters are set in final model, coefficient are chosen already
## use test set to predict the MSE
## as unclear how carret handle fitting and prediction, use carret only for parameter tuning in elastic net

###########################################
#### Lasso (tune lambda) ##################
###########################################
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
	pred.fit <- predict(cv.fit, s=lambda.choice, newx=x.test, type='class') # predict class for binomial

	# covariates chosen 
	coef.fit <- predict(cv.fit, s=lambda.choice, newx=x.test, type='coefficients')
	coef.chosen = which(coef.fit != 0) 
    coef.chosen = coef.chosen - 1 # index for 1 representing intercept

    pred.class <- as.numeric(pred.fit)
    ROC <- pROC::roc(y.test, pred.class)$auc # in fact should be called "AUC" value instead of ROC


	result = list(ROC=ROC, coef.chosen=coef.chosen)
	return(result)
}

###########################################
#### Elastic Net (tune lambda and alpha) ##
###########################################
# reference 1 (tune both parameters with caret)
##  https://stats.stackexchange.com/questions/268885/tune-alpha-and-lambda-parameters-of-elastic-nets-in-an-optimal-way
# reference 2 (extract final model and prediction with caret)
## https://topepo.github.io/caret/model-training-and-tuning.html


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

	#data.train = as.data.frame(cbind(y.train, x.train)) 
    data.train = data.frame(y.train, x.train) # avoid conversion of data types with binary response
	colnames(data.train) = c('y', paste('V', seq(1, dim(datx)[2]), sep=''))

	# tune parameters with CV on training data with caret
	trnCtrl <- trainControl(method = "cv", number = fold.cv)
	srchGrid <- expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)

    data.train$y = as.factor(data.train$y)
	my_train <- train(y ~., data.train,
              		  method = "glmnet", 
                      tuneGrid = srchGrid,
                      trControl = trnCtrl)

	# fit the model with best tuning parameters on test data with glmnet
	fit = glmnet(x.train, y.train, alpha=my_train$bestTune$alpha, 
                 lambda=my_train$bestTune$lambda, family=family)
	pred.fit <- predict(fit, x.test, type='class') # predict class for binomial
	
	# covariates chosen 
	coef.fit <- predict(fit, x.test, type='coefficients')
	coef.chosen = which(coef.fit != 0) 
    coef.chosen = coef.chosen - 1 # index for 1 representing intercept

	# evaluate results
    pred.class <- as.numeric(unlist(pred.fit))
    ROC <- pROC::roc(y.test, pred.class)$auc

	result = list(ROC=ROC, coef.chosen=coef.chosen)
    #result = list(pred.fit=pred.fit, y.test=y.test)
	return(result)
}


###############################################################
#### Random Forests (tune # variables at each random split ####
###############################################################
# reference on regression random forst
## https://uc-r.github.io/random_forests
# OOB error is different from test error (see above website)


randomForest_cv = function(datx, y, seednum=31, fold.cv=5, ratio.training=0.8, mtry.grid=10, num_trees=500, 
						   pval_thr=0.05, method.perm='altmann'){
	# mtry: number of variables to randomly sample at each split
	# num_trees: number of trees to grow in random forests
	# pval_thr: threshold for permutation test
	# note that permutation p-value can use "altmann method" for all types of data; 'Janita' for high-dimensitional data only
	# ref on permutation methods: http://finzi.psych.upenn.edu/library/ranger/html/importance_pvalues.html

	set.seed(seednum)  

	# split data into training and tests sets
	#data = as.data.frame(cbind(y, datx)) 
    data = data.frame(y, datx) # avoid conversion of data types with binary response
	colnames(data) = c('y', paste('V', seq(1, dim(datx)[2]), sep=''))
	inTraining = createDataPartition(data$y, p = ratio.training, list=FALSE)
	train <- data[inTraining, ]
	test <- data[-inTraining, ]  

	# tune parameter with cross validation
	hyper.grid <- expand.grid(mtry = mtry.grid, OOB_RMSE = 0)
	for (i in 1:nrow(hyper.grid)){
	  model = ranger(y ~., data = train,
	                 num.trees=500, mtry=hyper.grid$mtry[i],
	                 seed=seednum, importance = 'permutation')
	  hyper.grid$OOB_RMSE[i] = sqrt(model$prediction.error)
	}
	OOB = min(hyper.grid$OOB_RMSE) # out of bag error
	position = which.min(hyper.grid$OOB_RMSE)

	# permutation test on tuned random forst model to obtain chosen features
	if (method.perm == 'altmann'){ # for all data types
		rf.model <- ranger(y ~., data=test, num.trees = num_trees,
					       mtry = hyper.grid$mtry[position], importance = 'permutation')
	    table = as.data.frame(importance_pvalues(rf.model, method = "altmann", 
                              formula = y ~ ., data = test)) 
	} else if (method.perm == 'janitza'){ # for high dimensional data only
		rf.model <- ranger(y ~., data=test, num.trees = num_trees,
					       mtry = hyper.grid$mtry[position], importance = 'impurity_corrected')
	    table = as.data.frame(importance_pvalues(rf.model, method = "janitza", 
                              formula = y ~ ., data = test)) 		
	}

    coef.chosen = which(table$pvalue < pval_thr) 

    # if nothing been selected
    if (identical(coef.chosen, integer(0))){
    	coef.chosen = 0
    }

    # predicted class: https://www.rdocumentation.org/packages/ranger/versions/0.12.1/topics/predict.ranger
    pred_rf = predict(rf.model, test)
    pred.class <- as.numeric(pred_rf$predictions)
    ROC <- pROC::roc(test$y, pred.class)$auc


	result = list(mtry=mtry.grid[position], coef.chosen=coef.chosen, ROC=ROC, OOB=OOB, p.value=table$pvalue)
	return(result)

}


# # ###############################################################
# # #### Generalized Compositional Lasso by Lu et al 2019 #########
# # ###############################################################
dir_coda  = '../code_coda/Microbiome-Variable-Selection-master/Microbiome_variable_selection_tutorial/'

source(file = paste0(dir_coda, 'CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R'))

source(file = paste0(dir_coda, 'functions.R'))

gen_cons_lasso_cv = function(datx, y, seednum=31, data.split=FALSE,
                             ratio.training=0.8, lambda.coda=seq(0.1, 0.2, 0.01)){
    # note that generalized compositional lasso do the X transformation within
    # thus use sub$X instead of sub$Z for generalized compositonal lasso method
    
    set.seed(seednum)
    colnames(datx) = paste('V', seq(1, dim(datx)[2]), sep='')
    
    if (data.split == TRUE){
        # # split data into training and tests sets
        nn <- length(y)
        trn <- sample(1:nn, ratio.training*nn)
        x.train <- datx[trn, ]
        x.test <- datx[-trn, ]
        y.train <- y[trn]
        y.test <- y[-trn]
        
        lambda_table <- lambdaRange_codalasso(X = x.train, y = y.train, 
                                              lambdaSeq = lambda.coda)
        lambda_optimal <- lambda_table[which.max(lambda_table$prop.explained.dev),
                                       'lambda']
        codalasso_sim <- coda_logistic_lasso(X = x.test, y = y.test, 
                                             lambda = lambda_optimal)
        ROC <- pROC::roc(y.test, codalasso_sim$`predicted class`)$auc
        sim.results_codalasso <- coda_lasso_wrapper(result = codalasso_sim, 
                                                    X = x.test)
        coef.chosen = tidyr::extract_numeric(sim.results_codalasso$varSelect)
     }else{
        lambda_table <- lambdaRange_codalasso(X = datx, y = y, 
                                              lambdaSeq = lambda.coda)

        # choose lambda explained most deviance (if several, choose smaller lambda)
        lambda_optimal <- lambda_table[which.max(lambda_table$prop.explained.dev),
                                       'lambda']
        codalasso_sim <- coda_logistic_lasso(X = datx, y = y, 
                                             lambda = lambda_optimal)
        ROC <- pROC::roc(y, codalasso_sim$`predicted class`)$auc
        sim.results_codalasso <- coda_lasso_wrapper(result = codalasso_sim, 
                                                    X = datx)
        coef.chosen = tidyr::extract_numeric(sim.results_codalasso$varSelect)
     }   

        result = list(ROC=ROC, coef.chosen=coef.chosen)
        return(result)   
}    
    








