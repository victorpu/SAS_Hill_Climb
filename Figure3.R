#This file generate Figure3 to give a typical example of the weights that sparse K-means returns.
require(sparcl)
require(ggplot2)
set.seed(12)
mu = .7
p = 500
X1 = matrix(rnorm(n=50*20,mean=mu,sd=1),nrow=20)
X2 = matrix(rnorm(n=50*20,mean=0,sd=1),nrow=20)
X3 = matrix(rnorm(n=50*20,mean=-mu,sd=1),nrow=20)
X = rbind(X1,X2,X3)
if(p>50){
  X0 = matrix(rnorm((p-50)*60),nrow = 60)
  X = cbind(X,X0)
}
X = scale(X)

km.perm <- KMeansSparseCluster.permute(X,K=3,nperms=25,silent=TRUE)
km.out <- KMeansSparseCluster(X,K=3,wbounds=km.perm$bestw,silent=TRUE)

#plot(km.out[[1]]$ws)
#write.table(km.out[[1]]$ws,"weights.txt",row.names = FALSE,col.names = FALSE)
#weights = read.table("weights.txt")
weight_plot <- qplot(seq(1,500),km.out[[1]]$ws,xlab = NULL, ylab = NULL,size = I(3))
weight_plot + theme(
                    axis.text.x = element_text(size=16,face = "bold"),
                    axis.text.y = element_text(size=16,face = "bold"))