#This file generates data under the regime in 4.4, and computes the rand indexes resulted by SAS_gs and sparse-kmedoids.
#All the rand indexes are stored in "simulation_pam.txt", need to manually extract the data and compute mean and sd for each cell.
sink('simulation_pam.txt')
require(Rlab)
set.seed(234)
true_clust = c(rep(1,30),rep(2,30),rep(3,30))
for(mu in c(0.6,0.7,0.8,0.9,0.99)){
  for(p in c(30,60,100,150,200)){
    print(p)
    Rand_hill = rep(0,50)
    Rand_spam = rep(0,50)
    cat("Current loop: p=",paste(p),"mu=",paste(mu),'\n')
    for(b in seq(1,50)){
      X <- matrix(rbern(90*p, 0.1),nrow = 90)
      X[1:30,1:5] = matrix(rbern(5*30,mu),nrow = 30)
      X[31:60,6:10] = matrix(rbern(5*30,mu),nrow = 30)
      X[61:90,11:15] = matrix(rbern(5*30,mu),nrow = 30)
      
      clust1 = hill_climb_pam(X,k=3,nbins=p/2, nperms = 10,itermax = 10,threshold = 1)
      Rand_hill[b] = RRand(true_clust, clust1$best_result)$Rand
      
      pam.perm <- PamSparseCluster.permute(X,K=3,wbounds = seq(1.1,12,len=20),nperms = 10,silent = T)
      pam.out <- PamSparseCluster(X,K=3,wbounds=pam.perm$bestw,silent=TRUE)
      Rand_spam[b] = RRand(true_clust, pam.out[[1]]$Cs)$Rand
    }
    cat('Rand index for hill_climb:', paste(Rand_hill, collapse=','), ';\n')
    cat('Rand index for hill_spam:', paste(Rand_spam, collapse=','), ';\n')
  }
}
sink()