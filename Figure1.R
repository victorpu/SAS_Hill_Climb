#This file runs convergence test on Algorithm 1 using syntetic data generated from three mixture of sparse gaussians.
#We plot the means and confidence intervals of rand indexes and symmetric differences under different regimes (mu = 0.6, 0.7 and 0.8). For simplicity, we only show the regime when mu = 0.6/.
require(phyclust)
require(ggplot2)
require(grid)
require(gtable)
#Basically, the function 'Alternate' implements Algorithm 1 in the paper, but we manually set iteration number to be "itermax", and let itermax = 10 in this simulation.
Alternate = function(X, k,tot, initial_set, s, itermax){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  rand_index = rep(NA,itermax)
  diff = rep(NA,itermax)
  while(iternum< itermax){
    clustering = kmeans(X[,set0], centers=k,nstart=2)
    result = clustering$cluster
    group = list()
    cond = TRUE
    for(j in seq(1,k)){
      group[[j]] = which(result == j)
      cond = cond && length(group[[j]]) > 1 
    }
    center = NULL
    wcss = rep(0,p)
    if(cond){
      for(j in seq(1,k)){
        center = rbind(center,colMeans(X[group[[j]],]))
        Xc =  t(apply( X[group[[j]],], 1, function(x) x-center[j,]))
        wcss = wcss + apply(Xc, 2, function(x){sum(x^2)})
      }
      iternum = iternum + 1
    }
    set1 = set0
    set0 = which(rank((tot-wcss)/tot,ties.method = "random") > p-s)
    rand_index[iternum] = RRand(true_clust, result)$Rand
    diff[iternum] = length(setdiff(set1,SET)) + length(setdiff(SET,set1))
  }
  out = list(diff = diff, rand_index = rand_index,final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}


set.seed(234)
true_clust = c(rep(1,30),rep(2,30),rep(3,30))
mu = 0.6
p = 500
s = 50
k=3
SET = seq(1,50)
rand_matrix = matrix(NA,nrow = 100, ncol = 10)
diff_matrix = matrix(NA,nrow = 100, ncol = 10)
for(i in seq(1:100)){
  X1 = matrix(rnorm(n=50*30,mean=mu,sd=1),nrow=30)
  X2 = matrix(rnorm(n=50*30,mean=0,sd=1),nrow=30)
  X3 = matrix(rnorm(n=50*30,mean=-mu,sd=1),nrow=30)
  X = rbind(X1,X2,X3)
  if(p>50){
    X0 = matrix(rnorm((p-50)*90),nrow = 90)
    X = cbind(X,X0)
  }
  X = scale(X)
  n = dim(X)[1]
  p = dim(X)[2]
  center0 = colMeans(X)
  Xc0 = t(apply(X,1,function(x) x-center0))
  tot = apply(Xc0, 2, function(x) {sum(x^2)})
  wcss = rep(0,p)
  #Initialize the important feature set.
  for(j in 1:p){
    clustering = kmeans(X[,j], centers = k)
    wcss[j] = clustering$tot.withinss
  }
  rank0 = rank((tot-wcss)/tot,ties.method = "random")
  initial_set = which(rank0 > p-s)
  out = Alternate(X, k=3,tot, initial_set, s, itermax=10)
  rand_matrix[i,] = out$rand_index
  diff_matrix[i,] = out$diff
}

iter = seq(1,10)
write.table(rand_matrix, file = "convergence_rand0.txt",row.names = F, col.names = F)
write.table(diff_matrix, file = "convergence_diff0.txt",row.names = F, col.names = F)


data.rand <- read.table('convergence_rand0.txt')[,1:10]
mean.rand <- as.vector(sapply(data.rand, mean))
low.rand <- as.vector(sapply(data.rand, quantile, probs = 0.05))
high.rand <- as.vector(sapply(data.rand, quantile, probs = 0.95))

data.diff <- read.table('convergence_diff0.txt')[,1:10]
mean.diff <- as.vector(sapply(data.diff, mean))
low.diff <- as.vector(sapply(data.diff, quantile, probs = 0.05))
high.diff <- as.vector(sapply(data.diff, quantile, probs = 0.95))

#### Use ggplot2 to generate Figure 1######
# x: iter
# y1: mean.rand
# C1: low.rand, high.rand
# y2: mean.diff
# C2: low.diff, high.diff
grid.newpage()
df <- data.frame(x = iter, y1 = mean.rand, c1min = low.rand, c1max = high.rand, 
                 y2 = mean.diff, c2min = low.diff, c2max = high.diff)
p1 <- ggplot(df,aes(x = x, y = y1),xlab = NULL,ylab = NULL) + 
  scale_y_continuous(limits=c(0.55,1.0),breaks = c(0.6,0.7,0.8,0.9,1.0)) +
  scale_x_discrete(breaks=c("2","4","6","8","10"),limits = c("1","2","3","4","5","6","7","8","9","10"))+
  geom_errorbar(aes(ymin = c1min,ymax = c1max),width =.1,size = 1.0) +
  geom_line(size=1.4) + geom_point(size= 5,shape=21, fill="white") + 
   theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=25,face = "bold"),
    axis.text.y = element_text(size=25,face = "bold"))
p2 <- ggplot(df,aes(x = x, y = y2),xlab = NULL,ylab = NULL,colour = "red") + 
  scale_y_continuous(limits=c(0,100)) + 
  scale_x_discrete(breaks=c("2","4","6","8","10"),limits = c("1","2","3","4","5","6","7","8","9","10")) +
  geom_errorbar(aes(ymin = c2min,ymax = c2max),width =.1,size = 1.0,colour = "red",position="dodge") +
  geom_line(size=1.4,colour = "red") + geom_point(size= 5,shape=21, fill="white",colour = "red") +
 theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.x = element_text(size=25,face = "bold"),
                     axis.text.y = element_text(size=25,face = "bold",colour ="red"),
                     panel.background = element_rect(fill = NA))

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                     pp$l, pp$b, pp$l)

# axis tweaks
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

grid.draw(g)