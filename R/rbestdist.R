rbestdist<-function(dist,r,penalty=NULL,beta=100){
  stopifnot(length(dist$d)==length(r))
  maxr=max(r)
  if (!is.null(penalty)) stopifnot(length(penalty)==maxr)
  if (is.null(penalty)){
    penalty=numeric(maxr)
    if (maxr>1){
      for (j in 2:maxr){
        penalty[j]=(1+penalty[j-1])*beta
      }
    }
  }
  fbpenalty=(1+penalty[maxr])*beta
  for (k in 1:maxr){
    dist$d[which(r==k)]=dist$d[which(r==k)]+penalty[k]
  }
  return(list(newdist=dist, fbpenalty=fbpenalty))
}
