#With SparseClustering.R modified to be PamSparseClustering.R, we could replace kmeans in SparseClustering by pam, 
#so that we could perform sparse clustering on categorical data.

PamSparseCluster <- function(x, K=NULL, wbounds=NULL, silent=FALSE, maxiter=6){
    # The criterion is : minimize_{w, C} sum_j w_j (WCSS_j - TSS_j) s.t. ||w||_2=1, ||w||_1<=s, w_j>=0
    # x is the data, nxp
    # K is the number of clusters desired
    # wbounds is a vector of L1 constraints on w, of the form  sum(abs(w))<=wbounds[i]
    if(is.null(K)) stop("Must provide K.")
    if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(ncol(x)), len=20)
    if(min(wbounds)<=1) stop("wbounds should be greater than 1")
    wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
    out <- list()
    if(!is.null(K)){
      clust <- pam(hamming.distance(x),K,diss=TRUE)
      Cs <- clust$clustering
      medoids <- clust$medoids
    } 

    for(i in 1:length(wbounds)){
      if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
      ws <- rep(1/sqrt(ncol(x)), ncol(x)) # Start with equal weights on each feature
      ws.old <- rnorm(ncol(x))
      store.bcss.ws <- NULL
      niter <- 0
      while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
        if(!silent) cat(niter, fill=FALSE)
        niter <- niter+1
        ws.old <- ws
        if(niter>1){
          clust <- UpdateCs(x, K, ws, Cs,medoids) # if niter=1, no need to update!!
          Cs <- clust$Cs
          medoids = clust$medoids
        } 
        ws <- UpdateWs(x, Cs, medoids,wbounds[i])
        store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, Cs,medoids,ws)$bcss.perfeature*ws))
      }
      out[[i]] <- list(ws=ws, Cs=Cs, medoids = medoids,wcss=GetWCSS(x, Cs,medoids, ws),crit=store.bcss.ws, wbound=wbounds[i])
    }
    if(!silent) cat(fill=TRUE)
    #  if(length(wbounds)==1){
    #    out <- out[[1]]
    #    class(out) <- "kmeanssparse"
    #    return(out)
    #  }
    #  class(out) <- "multikmeanssparse"
    #  return(out)
    class(out) <- "PamSparseCluster"
    return(out)
  }
