#Simulation 1 new
sink('simulation1.txt')
set.seed(234)
true_clust = c(rep(1,20),rep(2,20),rep(3,20))
for(mu in c(0.6,0.7,0.8,0.9,1.0)){
	for(p in c(100,200,500,1000)){
    cat("Current loop: p=",paste(p),"mu=",paste(mu),'\n')
    Rand_hill = rep(0,50)
	Rand_gss = rep(0,50)
	Rand_gss_avg = rep(0,50)
    Rand_sparse = rep(0,50)
    for(b in seq(1,50)){
      X1 = matrix(rnorm(n=50*20,mean=mu,sd=1),nrow=20)
      X2 = matrix(rnorm(n=50*20,mean=0,sd=1),nrow=20)
      X3 = matrix(rnorm(n=50*20,mean=-mu,sd=1),nrow=20)
      X = rbind(X1,X2,X3)
      if(p>50){
        X0 = matrix(rnorm((p-50)*60),nrow = 60)
        X = cbind(X,X0)
      }
      X = scale(X)
	  #nbins here varies, when p=100, we use nbins = 100; when p = 200,500, we use nbins = p/2; when p = 1000, we let nbins = 200.
      fit = hill_climb(X,k=3,nperms=30,nbins=min(p/2,100),itermax = 20,threshold=0)            
      Rand_hill[b] = RRand(true_clust, fit$best_result)$Rand
      
	  fit1 = hill_climb_GSS(X,k=3,nperms=30,itermax=20,threshold=0,tolerance = 0)
	  Rand_gss[b] = RRand(true_clust, fit1$result)$Rand
      
      km.perm <- KMeansSparseCluster.permute(X,K=3,nperms=30,silent=TRUE)
      km.out <- KMeansSparseCluster(X,K=3,wbounds=km.perm$bestw,silent=TRUE)
      result1 = km.out[[1]]$Cs
      Rand_sparse[b] = RRand(true_clust, result1)$Rand

    }
    cat('Rand index for hill_climb:', paste(Rand_hill, collapse=','), ';\n')
	cat('Rand index for hill_climb_GSS:', paste(Rand_gss, collapse=','), ';\n')
    cat('Rand index for SparseKmeans:', paste(Rand_sparse, collapse=','), ';\n') 
  }
}
sink()
