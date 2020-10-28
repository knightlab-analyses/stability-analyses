########################################################################
### simulate data with Toeplitz correlation structure ##################
########################################################################
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

B = 5 # number of blocks
corr = 0.09 # across-block correlation (decrease from 0.1 to 0.09 to ensure positive definteness

### three different levels of correlation strength ###
rou.list = seq(0.1, 0.9, 0.2)


set.seed(31)
rep = 100


for (rou in rou.list){
    print(paste('rou', rou))
    for (dim in dim.list){
        print(dim)

        sim_array = NULL
        for (b in 1:rep){
            n = dim[2]
            p = dim[1]

            # construct covariate matrix X
            theta = vector(mode='numeric', length=p)
            theta[1:5] = rep(log(0.5*p), 5)
            theta[6:p] = rep(0, p-5)

            # block design covariance
            Sigma = matrix(rep(0, p*p), nrow=p)
            for (i in 1:p){
              for (j in 1:p){
                if (i == j){Sigma[i, j] = rou}
                if (i !=j & (i-j) %% (p/B) == 0) {Sigma[i, j] = corr}
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
            # beta: if 5 blocks, each block has one true signal, while the last block has two true signals
            beta = rep(0, p)
            beta[p/B] = 1; beta[p/B * 2] = -0.8; beta[p/B*3] = 0.6; 
            beta[p/B*4] = - 1.5; beta[c(p/B*5-1, p/B*5)] = c(-0.5, 1.2)
            
            Z = log(X)
            sigma = 0.5
            epsilon = rnorm(n=n, mean=0, sd=sigma)
            Y = Z %*% beta + epsilon

            sim_array[[b]] = list(rep_idx=b, rou=rou, p=p, n=n, W=W, X=X, Z=Z, beta=beta, sigma=sigma, Y=Y)

            save(file=paste0(dir, '/sim_block_corr', rou, paste('P', p, 'N', n, sep='_'), '.RData', sep=''), sim_array)
        }
    }
}


