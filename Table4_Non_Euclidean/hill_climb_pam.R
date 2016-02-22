#For within cluster distance, (a) our paper and W&T use average within-cluster dissimilarity.
#(b) hill_climb_pam also uses average within-cluster dissimilarity, but the dissimilarity measure is hamming distance.
#(c) Since our simulation in 4.4 is under the regime of binary variables, for simplicity, we only implement hill_climb_pam for binary variables,
# but it is very easy to generalize to multi-class categorical variables (just need to modify hamming.tot and hamming.within).

require(phyclust)
require(cluster)
require(e1071)
hamming.var = function(x){
  #input is a binary vector
  # output is the average hamming distance within the vector: sum_{i,i'}d_{i,i'}/length
  return(sum(x==0)*sum(x==1)/length(x))
}

hamming.tot = function(x){
  #input is a binary matrix
  #output is a vector of average hamming distances for columns of x
  apply(x,2,hamming.var)
}

hamming.within = function(X,k,G,medoids){
  if(is.vector(X) && is.atomic(X)){
    group = list()
    wcds = 0
    for(j in 1:k){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcds = wcds + hamming.var(X[group[[j]]])
      }
    }
  }else if (is.matrix(X)){
    p = dim(X)[2]
    group = list()
    wcds = rep(0,p)
    for(j in 1:k){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
	  wcds = wcds + hamming.tot(X[group[[j]],])
      }
    }
  } else {
    cat("X is niether a vector nor a matrix! \n")
    return(NULL)
  }
  return(wcds)
}

Alternate_pam = function(X, k,tot, initial_set, s, itermax, threshold){
  p = dim(X)[2]
  set0 = initial_set
  set1 = rep(0,s)
  iternum = 0
  while(iternum<= itermax && length(setdiff(set1,set0)) + length(setdiff(set0,set1)) > threshold ){
    D = hamming.distance(X[,set0])
    clustering = pam(D,k,diss=TRUE,cluster.only = FALSE)
    result = clustering$cluster
    medoids = clustering$medoids
    wcds = hamming.within(X,k,result,medoids)
    iternum = iternum + 1
    set1 = set0

    if(any(tot==0)){
      set0 = which(rank((tot-wcds+1)/(tot+1),ties.method = "random") > p-s)
    } else{
      set0 = which(rank((tot-wcds)/tot,ties.method = "random") > p-s)
    }
  }
  info = clustering$clusinfo
  out = list(final_set = set0, iternum = iternum, result = result, withinds = sum(hamming.within(X[,set0],k,result,medoids)))
  return(out)
}

Log1 <- function(x){
  y = x
  y[y==0] <- 1
  return(log(y))
}

hill_climb_pam= function(X,k,nbins=20,nperms=25,itermax=10,threshold=1){
  n = dim(X)[1]
  p = dim(X)[2]
  tot = hamming.tot(X)
  
  permx <- list()
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=n, ncol=p)
    for(j in 1:p) permx[[i]][,j] <- sample(X[,j])
  }  
  stepsize = p/nbins
  tots <- rep(NA,nbins)
  permtots <- matrix(NA, nrow=nbins, ncol=nperms)  
  outs = list()
  for(i in 2:nbins){
    s = floor(p - (i-1)*stepsize)
    initial_set =sample(1:p,s)
    cond = FALSE
    if(s>1) cond = nrow(unique(X[,initial_set]))>k
    if(s==1) cond = length(unique(X[,initial_set]))>k
    if(cond){
      out = Alternate_pam(X, k,tot, initial_set, s, itermax, threshold)
      outs[[i]] = out
      #computing betweenss = tot - withinss
	tots[i] = sum(hamming.tot(X[,out$final_set]))-out$withinds
#	print(out$final_set)
      for(t in 1:nperms){
	  Dperm = hamming.distance(permx[[t]][,out$final_set])
        permresult = pam(Dperm,k,diss=TRUE,cluster.only = FALSE)
        permtots[i,t] <- sum(hamming.tot(permx[[t]][,out$final_set]))- 
        sum(hamming.within(permx[[t]][,out$final_set],k,permresult$clustering,permresult$medoids))
      }
    }
  }
    gaps <- (Log1(tots)-apply(Log1(permtots),1,mean))
    idx = which.max(gaps)[1]
    all_info = list(tots = tots, permtots = permtots, idx= idx, outs = outs, gaps=gaps, feature_set = outs[[idx]]$final_set,
                    best_result = outs[[idx]]$result)
    return(all_info)
}