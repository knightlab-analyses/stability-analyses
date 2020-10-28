######################################################################################
### simulate data with Independent & Toeplitz correlation structure ##################
######################################################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]


library(MASS)

## different size of P (number of features) & N (number of samples)
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
rou.list = c(0, seq(0.1, 0.9, 0.2)) # corr = 0 indicates independent; the others for Toeplitz


set.seed(31)
rep = 100 # number of replicates for simulated data

for (rou in rou.list){
	print(paste('rou', rou))
    for (dim in dim.list){
    	print(dim)

        sim_array = NULL
        for (b in 1:rep){
            #print(paste('b', b))

            p = dim[1]
            n = dim[2]

            # construct covariate matrix X
            theta = vector(mode='numeric', length=p)
            theta[1:5] = rep(log(0.5*p), 5)
            theta[6:p] = rep(0, p-5)

            Sigma = matrix(rep(0, p*p), nrow=p)
            for (i in 1:p){
                for (j in 1:p){
                    Sigma[i, j] = rou^abs(i-j)
                }
            }

            W = mvrnorm(n = n, mu=theta, Sigma=Sigma)
            X = matrix(rep(0, n*p), nrow=n)
            for (i in 1:n){
                for (j in 1:p){
                    X[i, j] = exp(W[i, j]) / sum(exp(W[i, ]))
                }
            }

            # construct response Y
            beta = c(c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2), rep(0, p-8))
            Z = log(X)
            sigma = 0.5
            epsilon = rnorm(n=n, mean=0, sd=sigma)
            Y = Z %*% beta + epsilon

            sim_array[[b]] = list(rep_idx=b, rou=rou, p=p, n=n, W=W, X=X, Z=Z, beta=beta, sigma=sigma, Y=Y)

            if (rou == 0){
                save(file=paste0(dir, '/sim_independent_', paste('P', p, 'N', n, sep='_'), '.RData'), sim_array)

            } else {
                save(file=paste0(dir, '/sim_toeplitz_corr', rou, paste('P', p, 'N', n, sep='_'), '.RData', sep=''), sim_array)
            }
            
        }
    }
}

