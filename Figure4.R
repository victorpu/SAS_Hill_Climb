#Figure 4 uses a typical dataset from section 4.3 and we project the data into the first two principal components of the submatrix and the whole matrix.
#We also illustrate how well the clustering results by SAS_gs and SAS_gss compared with those by sparse Kmeans and IF-PCA.
require(parcor)
require(clusterGeneration)
require(sparcl)
require(bootSVD)
require(ggplot2)
##### Generate a dataset under the setup in Section 4.3 ####
set.seed(123)
Sigma1 = diag(seq(1,2,length.out = 500))
Sigma2 = diag(seq(2,3,length.out = 500))
Sigma3 = diag(seq(3,4,length.out = 500))
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
write.table(X0, file = "X.txt",row.names = F, col.names = F)
X0 = as.matrix(read.table("X.txt",header = FALSE))
#### Apply SAS_gs, SAS_gss and IF-PCA on X0 and obtain their clustering results and important feature_set.####
  
fit = hill_climb(X0,k=3,nperms=30,nbins=100,itermax = 20,threshold=0)            
fit1 = hill_climb_GSS(X0,k=3,nperms=30,itermax=20,threshold=0,tolerance = 0)
km.perm <- KMeansSparseCluster.permute(X0,K=3,nperms=25,silent=TRUE)
km.out <- KMeansSparseCluster(X0,K=3,wbounds=km.perm$bestw,silent=TRUE)
result3 = km.out[[1]]$Cs

#Figure4(a) sim2_pc1
pc2 = prcomp(X0[,1:50], center = TRUE, scale. = TRUE) 
pcs2 = pc2$rotation[,1:2]
proj2 = X0[,1:50]%*%pcs2
proj2 = scale(proj2)
plot(proj2[61:90,1],proj2[61:90,2],xlim=c(min(proj2[,1])-0.01,max(proj2[,1])+0.01 ), 
     ylim=c(min(proj2[,2])-0.01,max(proj2[,2])+0.01),
     xlab = '1st Principal Compoment',
     ylab = '2nd Principal Component',
     main = 'First Two Principal Components of X[, 1:50]',
     col = 'darkgreen'
)
points(proj2[1:30,1],proj2[1:30,2],col = 'blue')
points(proj2[31:60,1],proj2[31:60,2],col = 'red')

##Figure4(b) sim2_pc2
pc0 = prcomp(X0, center = TRUE, scale. = TRUE) 
pc01 = pc0$rotation[,1]
pc02= pc0$rotation[,2]
pc0s = cbind(pc01,pc02)
proj02 = X0%*%pc0s
proj02 = scale(proj02)
plot(proj02[61:90,1],proj02[61:90,2],xlim=c(min(proj02[,1])-0.01,max(proj02[,1])+0.01 ), 
     ylim=c(min(proj02[,2])-0.01,max(proj02[,2])+0.01 ),
     xlab = '1st Principle Compoment',
     ylab = '2nd Principle Component',
     main = 'First Two Principal Components of X',
     col = 'darkgreen'
)
points(proj02[1:30,1],proj02[1:30,2],col = 'blue')
points(proj02[31:60,1],proj02[31:60,2],col = 'red')

###Figure4(c)sim2_pc3
pc2 = prcomp(X0[,1:50], center = TRUE, scale. = TRUE) 
pcs2 = pc2$rotation[,1:2]
proj2 = X0[,1:50]%*%pcs2
proj2 = scale(proj2)
result1 = fit$best_result
#result1 = c(rep(1,30),rep(2,7),3,rep(2,3),1,rep(2,5),1,rep(2,4),1,rep(2,7),rep(3,8),2,rep(3,3),2,2,rep(3,5),2,rep(3,10))
group1 = which(result1 == 1)
group2 = which(result1 == 2)
group3 = which(result1 == 3)
plot(proj2[group2,1],proj2[group2,2],xlim=c(min(proj2[,1])-0.01,max(proj2[,1])+0.01 ), 
     ylim=c(min(proj2[,2])-0.01,max(proj2[,2])+0.01),
     xlab = '1st Principal Compoment',
     ylab = '2nd Principal Component',
     main = 'Clustering by SAS_gs',
     col = 'red'
)
points(proj2[group1,1],proj2[group1,2],col = 'blue')
points(proj2[group3,1],proj2[group3,2],col = 'darkgreen')

##Figure4(d) sim2_pc4
pc2 = prcomp(X0[,1:50], center = TRUE, scale. = TRUE) 
pcs2 = pc2$rotation[,1:2]
proj2 = X0[,1:50]%*%pcs2
proj2 = scale(proj2)
result2 = fit1$result
#result2 = c(rep(1,31),rep(2,6),3,1,2,2,1,2,2,1,2,2,1,rep(2,4),1,1,rep(2,5),1,2,2,3,3,3,3,2,3,2,3,3,3,2,2,rep(3,5),2,3,3,3,2,rep(3,6))
group1 = which(result2 == 1)
group2 = which(result2 == 2)
group3 = which(result2 == 3)
plot(proj2[group2,1],proj2[group2,2],xlim=c(min(proj2[,1])-0.01,max(proj2[,1])+0.01 ), 
     ylim=c(min(proj2[,2])-0.01,max(proj2[,2])+0.01),
     xlab = '1st Principal Compoment',
     ylab = '2nd Principal Component',
     main = 'Clustering by SAS_gss',
     col = 'red'
)
points(proj2[group1,1],proj2[group1,2],col = 'darkgreen')
points(proj2[group3,1],proj2[group3,2],col = 'blue')


##Figure4(e) sim2_pc5
pc2 = prcomp(X0[,1:50], center = TRUE, scale. = TRUE) 
pcs2 = pc2$rotation[,1:2]
proj2 = X0[,1:50]%*%pcs2
proj2 = scale(proj2)
result3 = c( 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1,2, 1, 2, 1, 2, 2, 1, 2,2, 1,
             1, 1, 1, 3, 1, 3, 1, 2, 2, 3, 2, 2, 3, 1, 3, 3, 1, 1,2, 2, 2, 1, 2, 2,
             2, 2, 1, 2, 3 ,2, 3, 1, rep(3,30))
group1 = which(result3 == 1)
group2 = which(result3 == 2)
group3 = which(result3 == 3)
plot(proj2[group2,1],proj2[group2,2],xlim=c(min(proj2[,1])-0.01,max(proj2[,1])+0.01 ), 
     ylim=c(min(proj2[,2])-0.01,max(proj2[,2])+0.01),
     xlab = '1st Principal Compoment',
     ylab = '2nd Principal Component',
     main = 'Clustering by Sparse K-means',
     col = 'red'
)
points(proj2[group1,1],proj2[group1,2],col = 'darkgreen')
points(proj2[group3,1],proj2[group3,2],col = 'blue')

##Figure4(f) sim2_pc6, the clustering results comes from IF-PCA using MATLAB codes IF-PCA and the data file "X.txt".
pc2 = prcomp(X0[,1:50], center = TRUE, scale. = TRUE) 
pcs2 = pc2$rotation[,1:2]
proj2 = X0[,1:50]%*%pcs2
proj2 = scale(proj2)
result4 = c(3,2,3,3,2,3,3,2,2,2,2,2,2,2,3,2,3,2,3,2,3,2,3,2,1,3,3,2,2,3,2,3,1,3,2,1,2,1,2,1,1,3,1,1,3,1,3,1,3,3,1,1,2,3,1,1,1,3,1,1,3,1,1,1,1,1,1,1,3,3,1,1,3,3,1,1,1,1,3,3,1,1,1,1,1,1,3,1,1,1) 
group1 = which(result4 == 1)
group2 = which(result4 == 2)
group3 = which(result4 == 3)
plot(proj2[group2,1],proj2[group2,2],xlim=c(min(proj2[,1])-0.01,max(proj2[,1])+0.01 ), 
     ylim=c(min(proj2[,2])-0.01,max(proj2[,2])+0.01),
     xlab = '1st Principal Compoment',
     ylab = '2nd Principal Component',
     main = 'Clustering by IF-PCA',
     col = 'blue'
)
points(proj2[group1,1],proj2[group1,2],col = 'darkgreen')
points(proj2[group3,1],proj2[group3,2],col = 'red')

################## Regenerate Figure 4, but instead of using plot, use ggplot ##########
true_result = c(rep(1,30),rep(2,30),rep(3,30))
P = cbind(proj2,result1,result2,result3,result4,true_result)
write.table(P,'proj_result.txt',row.names = FALSE,col.names = FALSE)
P.data = as.data.frame(P)

#Figure4(c) sas_gs
ggplot(P.data, aes(PC1, PC2,shape = factor(result1),color=factor(result1))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)

##Figure4(d) sas_gss
ggplot(P.data, aes(PC1, PC2,shape = factor(result2),color=factor(result2))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)

##Figure4(e) sparse_k
ggplot(P.data, aes(PC1, PC2,shape = factor(result3),color=factor(result3))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)

##Figure4(f) if_pca
ggplot(P.data, aes(PC1, PC2,shape = factor(result4),color=factor(result4))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)

##Figure4(a) true_clust
ggplot(P.data, aes(PC1, PC2,shape = factor(true_result),color=factor(true_result))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)
#Figure4(b): Projection onto the whole matrix components
pc0 = prcomp(X0, center = TRUE, scale. = TRUE) 
pc01 = pc0$rotation[,1]
pc02= pc0$rotation[,2]
pc0s = cbind(pc01,pc02)
proj02 = X0%*%pc0s
proj02 = scale(proj02)
P0 = cbind(proj02,true_result)
P0.data = as.data.frame(P0)
ggplot(P0.data, aes(pc01, pc02,shape = factor(true_result),color=factor(true_result))) + 
  geom_point(size = 3) +
  scale_colour_manual(values=manual_color) + theme(
    panel.border = element_rect(fill = "transparent",colour = "black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=16,face = "bold"),
    axis.text.y = element_text(size=16,face = "bold"),
    panel.background = element_rect(fill = NULL),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+ 
  guides(color=FALSE,shape = FALSE)
