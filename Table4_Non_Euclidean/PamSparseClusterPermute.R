#use monte carlo simulation to select the tunning parameter $s$.
PamSparseCluster.permute <- function(x, K=NULL,  nperms=25, wbounds=NULL,silent=FALSE, nvals=10){
    if(is.null(wbounds)) wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x))*.9), len=nvals))
    if(min(wbounds) <= 1) stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
    if(length(wbounds)<2) stop("Wbounds should be a vector of at least two elements.")
    # was seq(1.2, sqrt(ncol(x))*.6, len=10)
    if(is.null(K)) stop("Must provide K.")
    permx <- list()
    nnonzerows <- NULL
    for(i in 1:nperms){
      permx[[i]] <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
      for(j in 1:ncol(x)) permx[[i]][,j] <- sample(x[,j])
    }
    tots <- NULL
    #(x, K=NULL, wbounds=NULL, silent=FALSE, maxiter=6)
    out <- PamSparseCluster(x, K, wbounds=wbounds, silent=silent)
    for(i in 1:length(out)){
      nnonzerows <- c(nnonzerows, sum(out[[i]]$ws!=0))
      tots <- c(tots, out[[i]]$wcss$bcss.ws)
    }
    permtots <- matrix(NA, nrow=length(wbounds), ncol=nperms)
    for(k in 1:nperms){
      if(!silent) cat("Permutation ", k, "of ", nperms, fill=TRUE)
      perm.out <- PamSparseCluster(permx[[k]], K, wbounds=wbounds, silent=silent)
      for(i in 1:length(perm.out)){
        permtots[i,k] <- perm.out[[i]]$wcss$bcss.ws
      }
    }
    gaps <- (log(tots)-apply(log(permtots),1,mean))
    out <- list(tots=tots, permtots=permtots, nnonzerows=nnonzerows, gaps=gaps, sdgaps=apply(log(permtots),1,sd), wbounds=wbounds, bestw=wbounds[which.max(gaps)])
    if(!silent) cat(fill=TRUE)
    class(out) <- "PamSparseCluster.permute"
    return(out)
  }
