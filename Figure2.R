#This file runs simulation to generate Figure2 to illustrate the effectiveness of choosing $s$ using gap statistics.
# Need to run the "hill_climb.R" file prior to calling the fuction.
require(ggplot2)
set.seed(586)
mu = 1
p = 500
X1 = matrix(rnorm(n=50*30,mean=mu,sd=1),nrow=30)
X2 = matrix(rnorm(n=50*30,mean=0,sd=1),nrow=30)
X3 = matrix(rnorm(n=50*30,mean=-mu,sd=1),nrow=30)
X = rbind(X1,X2,X3)
if(p>50){
  X0 = matrix(rnorm((p-50)*90),nrow = 90)
  X = cbind(X,X0)
}
X = scale(X)
clust = hill_climb(X,k=3,nbins=500,nperms=20,itermax=10,threshold=0)
write.table(clust$gaps, file = "stats.txt",row.names = F, col.names = F)

gaps = as.matrix(read.csv("stats.txt",header = FALSE))
gaps = gaps[,1]

df <- data.frame(
  x = seq(500,1), y = gaps
)
base <- ggplot(df, aes(x = x, y = y),xlab = NULL,ylab = NULL)
base + geom_line(size = 1,colour = "blue")+theme(
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.x = element_text(size=16,face = "bold"),
  axis.text.y = element_text(size=16,face = "bold"))






