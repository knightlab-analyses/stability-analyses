#' ## source code: https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R
#' @title getStability
#' @description Nogueira's stability measure from a binary matrix representing 
#' feature selection results in M bootstrapped datasets using a machine-learning method.
#' @param X A binary matrix X of size M * p, where a row represents a feature set
#' (for one of M data sets) and a column represents the selection of a given 
#' feature over all the M data sets.
#' @param alpha The level of significance (e.g. if alpha=0.05, we will get 95% 
#' confidence intervals). It's an optional argument and is set to 5% by default.
#' @details 
#' @return A list of stability measures is returned. The list inlcude stability, 
#' the variance of statbility, the upper/lower bound of the (1-alpha) 
#' confidence interval.
#' @examples 
#' ## extreme cases example 
#' d = 2 # number of features
#' M = 10 # number of bootstrap replicates
#'
#' ## case 1: when stability index undefined -- Z all zeros or all ones (Nogueria2018 p.13)
#' Z_all_missed = matrix(rep(0, M*d), nrow=M) # since K_bar = 0, thus SI undefined
#' getStability(Z_all_missed)$stability
#'
#' Z_all_selected = matrix(rep(1, M*d), nrow=M) # since K_bar = d = 3, thus SI underfined
#' getStability(Z_all_selected)$stability
#'
#' ## case 2: when stability index reaches maximum 1 -- each column of Z either all ones or all zeros (but not only zeros or only ones)(Nogueria2018 p.13)
#' # this was the case when we got the wrong SI almost 1 for soil datasets: since sampled dataset the same across all
#' d = 10
#' for (i in 1: (d-1)){
#' 	d_ones = sample(seq(1, d, 1), i)
#' 	Z_tmp = matrix(rep(0, M*d), nrow=M)
#' 	Z_tmp[, d_ones] = 1 # since Sf = 0 for all features, thus SI = 1
#' 	SI = getStability(Z_tmp)$stability
#' 	print(i)
#' 	print(paste('SI:', SI, sep=''))
#' }
#'
#' ## case 3: when stability index near minimum 0 (appedix D: - 1/(M-1), but as M goes to infinity, minimum asymptotically 0)
#' ## when for each column of feature, it receives same numbers of 0 and 1
#' d_alt = rep(c(0, 1), M/2)
#' Z_alt = matrix(rep(d_alt, d), ncol=d)
#' getStability(Z_alt)$stability # -0.1111111, which is - 1/(M-1)
#
#' d_alt_2 = c(rep(0, M/2), rep(1, M/2))
#' Z_alt_2 = matrix(rep(d_alt_2, d), ncol=d)
#' getStability(Z_alt_2)$stability # -0.1111111, which is - 1/(M-1)
#
#' # when M goes to infinity
#' M = 10000
#' d_alt = rep(c(0, 1), M/2)
#' Z_alt = matrix(rep(d_alt, d), ncol=d)
#' getStability(Z_alt)$stability # -0.00010001, very close to 0
#' @export
getStability <- function(X, alpha=0.05) {
## the input X is a binary matrix of size M*d where:
## M is the number of bootstrap replicates
## d is the total number of features
## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
## it's an optional argument and is set to 5% by default
### first we compute the stability

M<-nrow(X)
d<-ncol(X)
hatPF<-colMeans(X) # selection probability of each feature 
kbar<-sum(hatPF)
v_rand=(kbar/d)*(1-kbar/d) # kbar is the sum of selection probability on all features; v_rand is like the variance of bernoulli dist
stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate

## then we compute the variance of the estimate
ki<-rowSums(X)
phi_i<-rep(0,M)
for(i in 1:M){ 
	phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
}
phi_bar=mean(phi_i)
var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate

## then we calculate lower and upper limits of the confidence intervals
z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval

return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))

}










