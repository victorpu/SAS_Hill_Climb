#This file compares SAS_gs, SAS_gss, sparse Kmeans and IFPCA under the regime in section 4.2.
#Also need to call hill_climb and hill_climb_GSS before running this code.
require(parcor)
require(clusterGeneration)
require(sparcl)
require(bootSVD)
set.seed(123)

Rand_hill = rep(0,50)
sdiff_hill = rep(0,50)
Rand_gss = rep(0,50)
sdiff_gss = rep(0,50)

Rand_sparse = rep(0,50)
sdiff_sparse = rep(0,50)
true_clust = c(rep(1,30),rep(2,30),rep(3,30))
S = seq(1,50)
Sigma1 = diag(seq(1,5,length.out = 500))
Sigma2 = Sigma1
Sigma3 = Sigma1


for(b in 1:50){

    U1 = genQ(500, lim_attempts = 200)
    
    Sigma1t = t(U1)%*%Sigma1%*%U1
    U2 = genQ(500, lim_attempts = 200)
    
    Sigma2t = t(U2)%*%Sigma2%*%U2
    U3 = genQ(500, lim_attempts = 200)
    
    Sigma3t = t(U3)%*%Sigma3%*%U3
    
    
    mu = rep(0,500)
    Y1= mvrnorm(n =30, mu, Sigma1t, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    Y2 = mvrnorm(n =30, mu, Sigma2t, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    Y3 = mvrnorm(n =30, mu, Sigma3t, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    Y = rbind(Y1,Y2,Y3)
    mu1 = c(seq(1.02,2,length.out = 50),rep(0,450))
    mu2 = c(seq(2.02,3,length.out = 50),rep(0,450))
    mu3 = c(seq(3.02,4,length.out = 50),rep(0,450))
	
    X = matrix(rep(0,90*500),nrow=90)
    for(i in 1:30){
      X[i,] = mu1
    }
    for(i in 31:60){
      X[i,] = mu2
    }
    for(i in 61:90){
      X[i,] = mu3
    }
    X0 = X + Y
    
    fit = hill_climb(X0,k=3,nperms=25,nbins=100,itermax = 20,threshold=0)            
    Rand_hill[b] = RRand(true_clust, fit$best_result)$Rand
    set_hill = fit$feature_set
    sdiff_hill[b] = length(setdiff(S,set_hill)) + length(setdiff(set_hill,S))
    
    fit1 = hill_climb_GSS(X0,k=3,nperms=25,itermax=20,threshold=0,tolerance = 0)
    Rand_gss[b] = RRand(true_clust, fit1$result)$Rand
    set_gss = fit1$final_set
    sdiff_gss[b] = length(setdiff(S,set_gss)) + length(setdiff(set_gss,S))
    
    km.perm <- KMeansSparseCluster.permute(X0,K=3,wbounds=seq(2,12,len=20),nperms=5,silent=TRUE)
    km.out <- KMeansSparseCluster(X0,K=3,wbounds=km.perm$bestw,silent=TRUE)
    result1 = km.out[[1]]$Cs
    Rand_sparse[b] = RRand(true_clust, result1)$Rand
    set_sparse = which(km.out[[1]]$ws>0)
    sdiff_sparse[b] = length(setdiff(S,set_sparse)) + length(setdiff(set_sparse,S))
}

cat("The average rand index resulted from SAS_gs is",paste(mean(Rand_hill)),";")
cat("The s.d. of rand index resulted from SAS_gs is",paste(sd(Rand_hill)),";")
cat("The average rand index resulted from SAS_gss is",paste(mean(Rand_gss)),";")
cat("The s.d. of rand index resulted from SAS_gss is",paste(sd(Rand_gss)),";")
cat("The average rand index resulted from sparse Kmeans is",paste(mean(Rand_sparse)),";")
cat("The s.d. of rand index resulted from sparse Kmeans is",paste(sd(Rand_sparse)),";")

cat("The average symmetric difference resulted from SAS_gs is",paste(mean(sdiff_hill)),";")
cat("The s.d. of symmetric difference resulted from SAS_gs is",paste(sd(sdiff_hill)),";")
cat("The average symmetric difference resulted from SAS_gss is",paste(mean(sdiff_gss)),";")
cat("The s.d. of symmetric difference resulted from SAS_gss is",paste(sd(sdiff_gss)),";")
cat("The average symmetric difference resulted from sparse Kmeans is",paste(mean(sdiff_sparse)),";")
cat("The s.d. of symmetric difference resulted from sparse Kmeans is",paste(sd(sdiff_sparse)),";")