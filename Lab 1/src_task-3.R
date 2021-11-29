library(mixtools)
library(car)

# 3.1
n = 1000
mu = c(2,3)
cov_mat = cbind(c(1, 1.4),c(1.4,4))
X = rmvnorm(n,mu,cov_mat)
scatterplot(X[,1],X[,2])

#3 .2
library(Rlab)
X = matrix(0,1000,2)
mu1 = c(2,3)
mu2 = c(3,2)

cov_mat1 = cbind(c(0.04, 0.06),c(0.06,0.36))
cov_mat2 = cbind(c(0.16, 0.06),c(0.06,0.09))

p = 0.6
n = 1000
B = rbern(n,p)
for (i in 1:1000) {
    if (B[i]==1) {
        x = rmvnorm(1,mu1,cov_mat1)
        X[i,1] = x[1]
        X[i,2] = x[2]
    }
    else if (B[i]==0) {
        x = rmvnorm(1,mu2,cov_mat2)
        X[i,1] = x[1]
        X[i,2] = x[2]
    }
}

# 3.3
mdl = mvnormalmixEM(X)