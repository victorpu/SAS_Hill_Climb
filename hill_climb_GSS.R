#This version is the same as hill_climb_GSS_revised.R in the folder of "Discarded"
require(sparcl)
require(mclust)
require(phyclust)

withinss = function(X,G,K){
  if(is.vector(X) && is.atomic(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + var(X[group[[j]]])*(l-1)
      }
    }
  } else if (is.matrix(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + sum(apply(X[group[[j]],],2,var)*(l-1))
      }   
    }
  } else {
    cat("X is niether a vector nor a matrix! \n")
    return(NULL)
  } 
  
  return(wcss)
}

Alternate= function(X, k,tot, initial_set, s, itermax, threshold){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  while(iternum<= itermax && length(setdiff(set1,set0)) + length(setdiff(set0,set1)) > threshold ){
    clustering = kmeans(X[,set0],iter.max = 20, centers=k,algorithm = "Hartigan-Wong",trace = 0,nstart=2)
    result = clustering$cluster
    
    wcss = apply(X,2,withinss,G = result, K = k)
    iternum = iternum + 1

    set1 = set0
    set0 = which(rank((tot-wcss)/tot,ties.method = "random") > p-s)
  }
  out = list(final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}

#compute within-cluster distance by clustering feature by feature, select S of size s based on this
hill_climb_GSS = function(X,k,nperms=20,itermax,threshold,tolerance){
  X = scale(X)
  n = dim(X)[1]
  p = dim(X)[2]
  tot = apply(X,2,var)*(n-1)
  permx <- list()
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=n, ncol=p)
    for(j in 1:p) permx[[i]][,j] <- sample(X[,j])
  }
  wcss = rep(0,p)
  for(j in 1:p){
    clustering = kmeans(X[,j],iter.max = 10, centers = k,algorithm = "Hartigan-Wong", trace = 0)
    wcss[j] = clustering$tot.withinss
  }
  rank0 = rank((tot-wcss)/tot,ties.method = "random")
  golden.ratio = 2/(sqrt(5) +1)
  iteration = 0
  upper.bound = p
  lower.bound = 1
  p1 = floor(upper.bound - golden.ratio*(upper.bound-lower.bound))
  p2 = floor(lower.bound + golden.ratio*(upper.bound-lower.bound))
  #evaluate the gap statistics using p1 and p2
  initial_set = which(rank0 > p-p1)
  out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap1 = log(out1$betweenss) - mean(log(permtots))
  initial_set = which(rank0 > p-p2)
  out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap2 = log(out2$betweenss) - mean(log(permtots))  
  
  while(abs(upper.bound - lower.bound) > tolerance)
  {
    iteration = iteration + 1
    if(gap2 < gap1) # then the maximum is to the left of x2
    {
      upper.bound = p2
      p2 = p1
      gap2 = gap1
      p1 = floor(upper.bound - golden.ratio*(upper.bound - lower.bound))
      #evaluate gaps for p1
      initial_set = which(rank0 > p-p1)
      out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap1 = log(out1$betweenss) - mean(log(permtots))
    } else {
      # the minimum is to the right of x1
      lower.bound = p1
      p1 = p2
      gap1 = gap2
      p2 = floor(lower.bound + golden.ratio * (upper.bound - lower.bound))
      #evaluate gaps for p2
      initial_set = which(rank0 > p-p2)
      out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap2 = log(out2$betweenss) - mean(log(permtots))
    }   
  }
  s = floor((lower.bound + upper.bound)/2)
  initial_set = which(rank0 > s)
  out = Alternate(X, k,tot, initial_set, s, 2*itermax, threshold)
  output = list(final_set = out$final_set, iternum = iteration, result = out$result, s = s)
  return(output)
}

true_clust = c(rep(1,20),rep(2,20),rep(3,20))
rand1 = rep(0,20)
mu = 0.7
p = 200

  X1 = matrix(rnorm(n=50*20,mean=mu,sd=1),nrow=20)
  X2 = matrix(rnorm(n=50*20,mean=0,sd=1),nrow=20)
  X3 = matrix(rnorm(n=50*20,mean=-mu,sd=1),nrow=20)
  X = rbind(X1,X2,X3)
  if(p>50){
    X0 = matrix(rnorm((p-50)*60),nrow = 60)
    X = cbind(X,X0)
  }
  
 # fit0 = hill_climb(X,k=3,nperms=25,nbins=p/5,itermax = 20,threshold=0)
  fit1 = hill_climb_GSS_revised(X,k=3,nperms=20,itermax=20,threshold=0,tolerance=1)
  fit2 = hill_climb_GSS(X,k=3,nperms=20,itermax=20,threshold=0,tolerance=1)
  RRand(true_clust, fit1$result)$Rand
  RRand(true_clust, fit2$result)$Rand


