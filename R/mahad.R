mahad<-function(z,X,p=rep(1,length(z)),caliper=1,stdev=FALSE,exact=NULL,nearexact=NULL,penalty=100,rank=TRUE){
  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x) || length(x)==length(z)){
      nocol=1
      x<-matrix(x,nrow=length(z))
    }else{
      nocol=ncol(x)
    }

    if(is.data.frame(x) || is.character(x)){
      if(!is.data.frame(x)) x <- as.data.frame(x)
      X.chars <- which(plyr::laply(x, function(y) 'character' %in% class(y)))
      if(length(X.chars) > 0){
        for(i in X.chars){
          x[,i] <- factor(x[,i])

        }
      }
      #if some variables are factors convert to dummies
      X.factors <-  which(plyr::laply(x, function(y) 'factor' %in% class(y)))

      #handle missing data
      for(i in which(plyr::laply(x, function(y) any(is.na(y))))){
        if(i %in% X.factors){
          #for factors, make NA a new factor level
          x[,i] <- addNA(x[,i])
        }else{
          #for numeric/logical, impute means and add a new indicator for missingness
          x[[paste(colnames(x)[i],'NA', sep = '')]] <- is.na(x[,i])
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }
      for(i in rev(X.factors)){
        dummyXi <- model.matrix(as.formula(
          paste('~',colnames(x)[i], '-1')),data=x)
        x <- cbind(x[,-i], dummyXi)
      }

    }else{
      #handle missing data
      for(i in 1:nocol){
        if(any(is.na(x[,i]))){
          x <- cbind(x,is.na(X[,i]))
          colnames(x)[nocol] <- paste(colnames(X)[i],'NA', sep = '')
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }

    }

    #get rid of columns that do not vary
    varying <- apply(x,2, function(y) length(unique(y)) > 1)
    x <- x[,which(varying),drop = FALSE]

    as.matrix(x)
  }

  #Check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(all((z==1)|(z==0)))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(length(z)==length(p))
  X<-Xmatrix(X)

  if (is.factor(nearexact)){
    levels(nearexact)=1:nlevels(nearexact)
    nearexact<-as.integer(nearexact)
  }

  stopifnot(length(z)==(dim(X)[1]))
  if (length(caliper)==1) stopifnot(caliper>=0)

  # Standardize p using an equally weighted pooled variance
  if (stdev){
    v1<-stats::var(p[z==1])
    v2<-stats::var(p[z==0])
    sp<-sqrt((v1+v2)/2)
    caliper<-caliper*sp
  }

  #Must have treated first
  if(!(min(z[1:(nobs-1)]-z[2:nobs])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    X<-X[o,]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(nearexact)) nearexact<-nearexact[o]
  }

  if (is.vector(X)) X<-matrix(X,ncol=1)

  ids<-1:nobs
  k<-dim(X)[2]

  if (rank){
    for (j in 1:k) X[, j]<-rank(X[, j])
  }
  #get rid of duplicated columns
  X <- X[,!duplicated(t(X))]
  #keep only linearly independent columns
  X<-cbind(rep(1,length(z)),X)
  q <- qr(X)
  X <- X[,q$pivot[seq(q$rank)]]
  X <- X[,-1]
  if (is.vector(X)) X=matrix(X,nrow=nobs)
  cv<-stats::cov(X)
  if (rank){
    vuntied<-stats::var(1:nobs)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
  }
  LL<-chol(cv)

  if (is.vector(X)) X<-matrix(X,ncol=1)
  X0<-X[z==0,]
  X1<-X[z==1,]
  if (is.vector(X0)) X0<-matrix(X0,ncol=1)
  if (is.vector(X1)) X1<-matrix(X1,ncol=1)
  cids <-ids[z==0]
  p0<-p[z==0]
  p1<-p[z==1]
  if (!is.null(exact)){
    exact0<-exact[z==0]
    exact1<-exact[z==1]
  }
  if (!is.null(nearexact)){
    nearex0<-nearexact[z==0]
    nearex1<-nearexact[z==1]
  }

  #  if (length(caliper)==1) caliper<-c(-caliper,caliper)
  #  if (length(caliper)==2) caliper<-sort(caliper)

  edgen<-DiPs::edgenum(z,p,caliper,NULL,exact)
  distance<-numeric(edgen)
  start_node<-numeric(edgen)
  end_node<-numeric(edgen)
  if (!is.null(nearexact)){
    nearex<-numeric(edgen)
  }else{
    nearex<-c()
  }

  current<-0

  for (i in 1:ntreat){
    #use caliper
    d<-abs(p1[i]-p0)
    who<-d<=caliper
    #use exact
    d<-abs(p1[i]-p0)
    if (!is.null(exact)) who<-who&(exact1[i]==exact0)
    num<-sum(who)
    if (any(who)){
      cc<-X0[who,]
      if (num==1) cc<-matrix(cc,nrow=1)
      else if (is.vector(cc)) cc<-matrix(cc,nrow=num)
      tt<-t(as.matrix(X1[i,]))
      distancei<-mvnfast::maha(cc,tt,LL,isChol=TRUE)
    }



    if (num>0){
      distance[(current+1):(current+num)]<-distancei
      start_node[(current+1):(current+num)]<-rep(i, num)
      end_node[(current+1):(current+num)]<-cids[who]

      #use nearexact
      if (!is.null(nearexact)) nearex[(current+1):(current+num)]<-nearex1[i]!=nearex0[who]
    }

    current<-current+num
  }
  if (!is.null(nearexact)){
    out<-list(d=distance+nearex*penalty,start=start_node,end=end_node)
  }else{
    out<-list(d=distance,start=start_node,end=end_node)
  }
  out
}
