####On the dataset of Brain####
X = read.table("brain.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.table("brain.y.txt",header = FALSE)
y = t(y)
X = scale(X)
hill = hill_climb(X,k=5,nperms=25,nbins=550,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=5,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=5,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=5,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Breast ###
X = read.table("breast.x.txt",header = FALSE)
y = read.table("breast.y.txt",header = FALSE)
X = t(as.matrix(X))
hill = hill_climb(X,k=2,nperms=25,nbins=2200,itermax = 20,threshold=0) 
hii$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=25,itermax = 10,threshold=0,tolerance = 0) 
hill_GSS$best_result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Colon ###
X = read.table("colon.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.table("colon.y.txt",header = FALSE)
hill = hill_climb(X,k=2,nperms=25,nbins=2000,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=50,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Lung_gordon ###
X1 = read.table("lung_gordon1.x.txt",header = FALSE)
X2 = read.table("lung_gordon2.x.txt",header = FALSE)
X3 = read.table("lung_gordon3.x.txt",header = FALSE)
X4 = read.table("lung_gordon4.x.txt",header = FALSE)
X = rbind(t(as.matrix(X1)),t(as.matrix(X2)),t(as.matrix(X3)),t(as.matrix(X4)))
y = c(rep(1,150),rep(2,31))
hill = hill_climb(X,k=2,nperms=25,nbins=1200,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=50,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
km.perm <- KMeansSparseCluster.permute(X,K=2,nperms=5)
km.out <- KMeansSparseCluster(X,K=2,wbounds=km.perm$bestw)
km.out[[1]]$Cs

#### On the dataset of Lung(2) ###
X = read.table("lung2.x.txt",header = FALSE)
X = t(as.matrix(X))
y = c(rep(0,139),rep(1,64))
hill = hill_climb(X,k=2,nperms=25,nbins=1260,itermax = 20,threshold=0) 
hii$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=50,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Leukemia ###
X = read.table("leukemia.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.csv("leukemia.y.txt",header = FALSE)
y = t(y)
hill = hill_climb(X,k=2,nperms=25,nbins=350,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Lymphoma ###
X = read.table("lym.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.table("lym.y.txt",header = FALSE)

hill = hill_climb(X,k=3,nperms=25,nbins=400,itermax = 20,threshold=0) 
hill$best_result

hill_GSS = hill_climb_GSS(X,k=3,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=3,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=3,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of Prostate ###
X = read.table("prostate.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.table("prostate.y.txt",header = FALSE)

hill = hill_climb(X,k=2,nperms=25,nbins=600,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=2,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=pam.perm$bestw,silent=TRUE)
out[[1]]$Cs

#### On the dataset of SRBCT ###
X = read.table("srbct.x.txt",header = FALSE)
X = t(as.matrix(X))
y = read.table("srbct.y.txt",header = FALSE)
perm <- KMeansSparseCluster.permute(X,K=4,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=4,wbounds=perm$bestw,silent=TRUE)
out[[1]]$Cs
hill = hill_climb(X,k=4,nperms=25,nbins=230,itermax = 20,threshold=0) 
hill$best_result
hill_GSS = hill_climb_GSS(X,k=4,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result

#### On the dataset of Sucancer ###
X = read.table("SuCancer.txt",header = FALSE)
y = read.table("SuCancer.y.txt", header = FALSE)
X = t(as.matrix(X))

hill = hill_climb(X,k=2,nperms=25,nbins=790,itermax = 20,threshold=0) 
hill$best_result

hill_GSS = hill_climb_GSS(X,k=2,nperms=25,itermax = 20,threshold=0,tolerance = 0) 
hill_GSS$result
perm <- KMeansSparseCluster.permute(X,K=2,nperms = 25,silent = T)
out <- KMeansSparseCluster(X,K=2,wbounds=perm$bestw,silent=TRUE)
out[[1]]$Cs


