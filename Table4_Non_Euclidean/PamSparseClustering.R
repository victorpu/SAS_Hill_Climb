#Modification of SparseClustering.R so that it could perform high-dimensional clustering on categorical data
#Basically, we need to modify their dissimilarity measure and functions of GetWCSS and related functions.
GetWCSS <- function(x, Cs, medoids,ws=NULL){
  k = length(unique(Cs))
  wcss.perfeature = hamming.within(X,k,Cs,medoids)
  
  bcss.perfeature <- hamming.tot(x)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature, bcss=sum(bcss.perfeature),bcss.ws=sum(bcss.perfeature*ws)))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}

weighted.hamming.distance <- function(x,y=NULL,ws=NULL){
  z<-NULL
  if(is.vector(x) && is.vector(y) && is.null(ws)){
    z <- sum(x != y)
  }
  else if(is.vector(x) && is.vector(y) && is.vector(ws)){
    z <- sum(ws[x != y])
  }
  else if(is.vector(x) && is.null(y)){
    z <- matrix(0,nrow=length(x),ncol=length(x))
    for(k in 1:(length(x)-1)){
      for(l in (k+1):length(x)){
        z[k,l] <- 1*(x[k]!=x[l])
        z[l,k] <- z[k,l]
      }
    }
  }
  else{
    z <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(k in 1:(nrow(x)-1)){
      for(l in (k+1):nrow(x)){
        z[k,l] <- weighted.hamming.distance(x[k,],x[l,],ws)
        z[l,k] <- z[k,l]
      }
    }
    dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
  }
  return(z)
}

UpdateCs <- function(x, K, ws, Cs,medoids){
  x <- x[,ws!=0]
  D <- weighted.hamming.distance(x=x,ws=ws[ws!=0])
  if(is.null(Cs)){
    km <- pam(D, K,diss=TRUE)
  } else if(!is.null(Cs) && !is.null(medoids)){
    km <- pam(D, K,diss=TRUE,medoids = medoids)
  }
  return(list(Cs = km$clustering,medoids = km$medoids ))
}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

positive.part <- function(x){
  x[which(x<0)] = 0
  return(x)
}

UpdateWs <- function(x, Cs, medoids, l1bound){
  #GetWCSS <- function(x, Cs, medoids,ws=NULL)
  wcss.perfeature <- GetWCSS(x, Cs, medoids)$wcss.perfeature
  tss.perfeature <- hamming.tot(x)
 # aplus = positive.part(-wcss.perfeature+tss.perfeature)
  aplus = -wcss.perfeature+tss.perfeature
  lam <- BinarySearch(aplus, l1bound)
  ws.unscaled <- soft(aplus,lam)
  return(ws.unscaled/l2n(ws.unscaled))
}

### up to here ###
