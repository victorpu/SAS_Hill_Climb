require(sparcl)
require(mclust)
require(phyclust)

Alternate = function(X, k,tot, initial_set, s, itermax, threshold){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  while(iternum<= itermax && length(setdiff(set1,set0)) + length(setdiff(set0,set1)) > threshold ){
    clustering = kmeans(X[,set0],iter.max = 20, centers=k, 
                        algorithm = "Hartigan-Wong",trace = 0,nstart=2)
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
  }
  out = list(final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}

#compute within-cluster distance by clustering feature by feature, initialize S of size s based on this
hill_climb = function(X,k,nbins=50,nperms=25,itermax,threshold){
  n = dim(X)[1]
  p = dim(X)[2]
  center0 = colMeans(X)
  Xc0 = t(apply(X,1,function(x) x-center0))
  tot = apply(Xc0, 2, function(x) {sum(x^2)})
  permx <- list()
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=n, ncol=p)
    for(j in 1:p) permx[[i]][,j] <- sample(X[,j])
  }
  wcss = rep(0,p)
  for(j in 1:p){
    if(length(unique(X[,j]))<k){
      wcss[j] = tot[j]
    }else{
      clustering = kmeans(X[,j],iter.max = 10, centers = k,algorithm = "Hartigan-Wong", trace = 0)
      wcss[j] = clustering$tot.withinss
    } 
  }
  stepsize = p/nbins
  rank0 = rank((tot-wcss)/tot,ties.method = "random")
  tots <- NULL
  permtots <- matrix(NA, nrow=nbins, ncol=nperms)
  tots <- c(tots, sum(tot-wcss))
  outs = list()
  for(i in 2:nbins){
    s = floor(p - (i-1)*stepsize)
    initial_set = which(rank0 > p-s)
    out = Alternate(X, k,tot, initial_set, s, itermax, threshold)
    outs[[i]] = out
    tots <- c(tots, out$betweenss)
    for(t in 1:nperms){
      permresult = kmeans(permx[[t]][,out$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
      permtots[i,t] <- permresult$betweenss
    }
  }
  gaps <- (log(tots)-apply(log(permtots),1,mean))
  idx = which.max(gaps)[1]
  all_info = list(idx= idx, outs = outs, gaps=gaps, feature_set = outs[[idx]]$final_set,
                  best_result = outs[[idx]]$result)
  return(all_info)
}
