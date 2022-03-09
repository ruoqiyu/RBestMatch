rbestmatch<-function(z,data,Xvar=NULL,Pvar=NULL,score=NULL,exact=NULL,fine=NULL,calipers=NULL,cutdist=NULL,nranks=length(calipers)-1,ncontrol=1,rank=TRUE,s.cost=10,rpenalty=50){
  #Check input
  stopifnot(is.data.frame(data))
  stopifnot(is.vector(z))
  if (is.null(fine)) fine=rep(1,length(z))
  if (is.factor(fine)){
    levels(fine)<-1:nlevels(fine)
    fine<-as.integer(fine)
  }
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  n<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  stopifnot(length(z)==length(fine))
  stopifnot(length(z)==nrow(data))
  stopifnot((!is.null(Xvar))|(!is.null(Pvar))|(!is.null(score)))

  if (is.null(Xvar)) stopifnot(!is.null(score))
  else if (is.null(score)) score<-stats::glm(z~Xvar,family=stats::binomial)$fitted.values

  if (is.null(Pvar)){
    if (!is.null(Xvar)) Pvar=Xvar
    else{
      Pvar=score
      nvar=1
    }
  }
  nvar=ncol(Pvar)

  sp<-sqrt((stats::var(score[z==1])+stats::var(score[z==0]))/2)
  if (is.null(calipers)){
    calipers=c(0.1*sp,0.2*sp,0.3*sp,0.4*sp,1)
    nranks=4
  }
  if (is.null(cutdist)) cutdist=c(0.1*2*nvar,0.2*2*nvar,0.5*2*nvar,2*nvar,Inf)
  dist=mahad(z,Pvar,score,caliper=calipers[nranks],stdev=FALSE,exact=exact,nearexact=NULL,penalty=100,rank=rank)
  score1=score[z==1]
  score0=score[z==0]
  scored=abs(score1[dist$start]-score0[dist$end-sum(z)])
  r=rep(nranks+1,length(dist$d))
  for (ri in nranks:1){
    r[which((scored<=calipers[ri]) & (dist$d<=cutdist[ri]))]=ri
  }
  whichr=which(r<=nranks)
  dist$d=dist$d[whichr]
  dist$start=dist$start[whichr]
  dist$end=dist$end[whichr]
  rr=r[whichr]
  rbdist=rbestdist(dist,rr,penalty=NULL,beta=rpenalty)
  dist=rbdist$newdist
  penalty=rbdist$fbpenalty

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    data<-data[o,]
    fine<-fine[o]
  }

  if (!requireNamespace("optmatch", quietly=TRUE)) {
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run match command, you must install optmatch first and agree to the terms of its license.")
  }

  ##building network
  #create basic treated-vs-control bipartite graph
  fine1=fine[z==1]
  fine0=fine[z==0]
  startn<-dist$start
  endn<-dist$end
  cost<-dist$d
  tcarcs<-length(startn) # number of treatment-control arcs
  ucap<-rep(1,tcarcs)
  b<-rep(ncontrol,ntreat) #supply for treated nodes
  b<-c(b,rep(0,ncontr)) #flow conservation at control nodes
  #Make costs integer
  if (any(cost<0)){
    lb=min(cost)
    cost<-cost-lb
    penalty<-penalty-lb
  }

  #create a duplicate for each control to make sure each control is only used once
  startn=c(startn,(ntreat+1):n)
  endn=c(endn,(n+1):(n+ncontr))
  cost=c(cost,rep(0,ncontr))
  ucap=c(ucap,rep(1,ncontr))
  b<-c(b,rep(0,ncontr))

  #Add structure to the bipartite graph for near fine balance
  tb<-table(z,fine)
  nt<-as.vector(tb[2,])
  nwant<-nt*ncontrol #desired number
  finelevels<-as.vector(as.numeric(colnames(tb)))

  bypassix=c()
  useix=length(startn)

  #Add a node for fine balance category k
  sinks<-NULL
  for (k in 1:(dim(tb)[2])){
    sinkk<-length(b)+1
    sinks<-c(sinks,sinkk)
    b<-c(b,-nt[k]*ncontrol)
    who0<-fine0==finelevels[k]
    if (sum(who0)>0){
      startn<-c(startn,rep(n,sum(who0))+which(who0))
      endn<-c(endn,rep(sinkk,sum(who0)))
      ucap<-c(ucap,rep(1,sum(who0)))
      cost<-c(cost,rep(0,sum(who0)))
    }

    who1<-fine1==finelevels[k]
    if (sum(who1)>0){
      startn<-c(startn,which(who1))
      endn<-c(endn,rep(sinkk,sum(who1)))
      ucap<-c(ucap,rep(ncontrol,sum(who1)))
      cost<-c(cost,rep(penalty,sum(who1)))
    }
    bypassix=c(bypassix,useix+sum(who0)+(1:sum(who1)))
    useix=useix+sum(who0)+sum(who1)
  }

  #Make costs integer
  cost<-round(cost*s.cost)
  net<-list(startn=startn,endn=endn,ucap=ucap,b=b,cost=cost,tcarcs=tcarcs)

  if (any(net$cost==Inf)) net$cost[net$cost==Inf]<-2*max(net$cost[net$cost!=Inf])

  #do match
  callrelax <- function(net){
    if (!requireNamespace("optmatch", quietly = TRUE)){
      stop('Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.')
    }
    startn <- net$startn
    endn <- net$endn
    ucap <- net$ucap
    b <- net$b
    cost <- net$cost
    nnodes <- length(b)
    my.expr <- parse(text = '.Fortran("relaxalg", nnodes, as.integer(length(startn)),
                     as.integer(startn), as.integer(endn), as.integer(cost),
                     as.integer(ucap), as.integer(b), x1 = integer(length(startn)),
                     crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4),
                     feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = "optmatch")')
    fop <- eval(my.expr)
    x <- fop$x1
    feasible <- fop$feasible1
    crash <- fop$crash1
    list(crash = crash, feasible = feasible, x = x)
  }

  #output<-rcbalance::callrelax(net)
  output<-callrelax(net)

  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    bypassflow=output$x[bypassix]
    if (ncontrol>1 & sum(bypassflow)>0){
      m<-list(feasible=output$feasible,d=NULL)
      warning("Fine balance is not feasible. Try increasing the number of edge sets to use (nranks) or decreasing the number of controls matched to each treated unit (ncontrol)")
    }
    else{
      x<-output$x[1:net$tcarcs]
      treated<-net$startn[1:net$tcarcs]
      control<-net$endn[1:net$tcarcs]
      treated=treated[which(x==1)]
      control=control[which(x==1)]
      match.df=data.frame('treat'=treated,'control'=control)
      match.df$treat<-as.factor(as.character(match.df$treat))
      matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                     function(treat.edges) treat.edges$control,.drop_o=FALSE))
      id1<-(1:n)[z==1]
      id0<-(1:n)[z==0]
      matchid<-matrix(c(id1[as.numeric(row.names(matches))], id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
      matchid<-as.vector(t(matchid))
      dat1<-data[matchid,]
      zm<-z[matchid]
      mset<-rep(1:nrow(matches),each=ncontrol+1)
      dat1$mset<-mset
      nmatchedpairs=dim(matches)[1]
      m<-list(data=dat1,nmatchedpairs=nmatchedpairs,rank_before=r,rank_after=rr[x==1])
    }
  }
  m
}
