\name{rbestmatch}
\alias{rbestmatch}
\title{Constructs a r-best matching.}
\usage{
rbestmatch(z,data,Xvar=NULL,Pvar=NULL,score=NULL,exact=NULL,fine=NULL,
calipers=NULL,cutdist=NULL,nranks=length(calipers)-1,ncontrol=1,
rank=TRUE,s.cost=10,rpenalty=50)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{data}{ A dataframe with length(z) rows containing the treatment indicator, all observed covariates and outcomes of interest. }
\item{Xvar}{ A matrix with length(z) rows giving all covariates to estimate the propensity score. If NULL, user must supply a score. }
\item{Pvar}{ A matrix with length(z) rows giving the priority covariates to calculate the Mahalanobis distance. X should be of full column rank. If NULL, Pvar=Xvar; if Xvar is also NULL, Pvar=score. }
\item{score}{ A vector of score for which calipers are applied for, e.g., score can be chosen as the propensity score. }
\item{exact}{ If not NULL, then a vector of length length(z) giving the nominal variable to be exactly matched. }
\item{fine}{ If not NULL, then a vector of length length(z) giving the nominal variable to be finely balanced. }
\item{calipers}{ A vector of calipers for the score, defining the ranks.  If calipers=NULL and let psd denote the equally weighted pooled standard deviation of score, use (.1*psd,.2*psd,.3*psd,0.4*psd,1). }
\item{cutdist}{ A vector of distance values that define the ranks. If cutdist=NULL and there are k=ncol(Pvar) covariates to calculate the distance, use (0.1*2*k,0.2*2*k,0.5*2*k,2*k,Inf). }
\item{nranks}{ Use up to nranks of the best edge sets. }
\item{ncontrol}{ Number of controls matched to each treated individual.}
\item{rank}{ A boolean denoting whether the distance is the standard Mahalanobis distance or its robust version. The default choice is TRUE to calculate the rank-based Mahalanobis distance.}
\item{s.cost}{ Scale of all distances (cost) for the network flow problem.}
\item{rpenalty}{ The tuning parameter to penalize the distances for r-best matching.}
}

\value{
\item{data}{ A dataframe containing selective rows of the original dataframe, corresponding to the matched individuals. The matched dataframe contains one additional column mset denoting which matched set each individual belongs to.}
\item{nmatchedpairs}{ Number of matched sets in the matched data.}
\item{rank_before}{ Rank for all edges before matching.}
\item{rank_after}{ Rank for the edges corresponding to the matched sets.}
}
\description{
Computes a Mahalanobis distance (or its robust version) list for use in matching. In this case, we compute the distance for all possible pairs of treated and control within caliper for score p.

This function and its use are discussed in Rosenbaum (2010).  The robust
Mahalanobis distance in described in Chapter 8 of Rosenbaum (2010).

}
\details{
The usual Mahalanobis distance works well for multivariate Normal covariates, but can exhibit odd behavior with typical covariates. Long tails or an outlier in a covariate can yield a large estimated
variance, so the usual Mahalanobis distance pays little attention to large differences in this covariate. Rare binary covariates have a small variance, so a mismatch on a rare binary covariate is
viewed by the usual Mahalanobis distance as extremely important. If you were matching for binary covariates indicating US state of residence, the usual Mahalanobis distance would regard a
mismatch for Wyoming as much worse than a mismatch for California.

The robust Mahalanobis distance uses ranks of covariates rather than the covariates themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a variable more important. Binary covariates are, of course, heavily tied.
}
\examples{
# To run this example, you MUST install the optmatch package.
# The optmatch match is available in the archived repository of CRAN and the Github page <https://github.com/markmfredrickson/optmatch>
\donttest{
install.packages("devtools")

library(devtools)
install.packages("optmatch", repos = "https://cran.microsoft.com/snapshot/2022-02-15/")
install_github("ruoqiyu/DiPs")
install_github("ruoqiyu/RBestMatch")

library(optmatch)
library(RBestMatch)

### RHC example in the manuscript

##data cleaning
data(rhc)
z=as.numeric(rhc$swang1=='RHC')
attach(rhc)
female=as.numeric(sex=="Female")
race_black=as.numeric(race=="black")
race_other=as.numeric(race=="other")
income1=as.numeric(income=='$11-$25k')
income2=as.numeric(income=='$25-$50k')
income3=as.numeric(income=='> $50k')
ins_care=as.numeric(ninsclas=='Medicare')
ins_pcare=as.numeric(ninsclas=='Private & Medicare')
ins_caid=as.numeric(ninsclas=='Medicaid')
ins_no=as.numeric(ninsclas=='No insurance')
ins_carecaid=as.numeric(ninsclas=='Medicare & Medicaid')
cat1_copd=as.numeric(cat1=='COPD')
cat1_mosfsep=as.numeric(cat1=='MOSF w/Sepsis')
cat1_mosfmal=as.numeric(cat1=='MOSF w/Malignancy')
cat1_chf=as.numeric(cat1=='CHF')
cat1_coma=as.numeric(cat1=='Coma')
cat1_cirr=as.numeric(cat1=='Cirrhosis')
cat1_lung=as.numeric(cat1=='Lung Cancer')
cat1_colon=as.numeric(cat1=='Colon Cancer')
cat2_mosfsep=as.numeric(match(cat2,'MOSF w/Sepsis',nomatch = 0)>0)
cat2_coma=as.numeric(match(cat2,'Coma',nomatch = 0)>0)
cat2_mosfmal=as.numeric(match(cat2,'MOSF w/Malignancy',nomatch = 0)>0)
cat2_lung=as.numeric(match(cat2,'Lung Cancer',nomatch = 0)>0)
cat2_cirr=as.numeric(match(cat2,'Cirrhosis',nomatch = 0)>0)
cat2_colon=as.numeric(match(cat2,'Colon Cancer',nomatch = 0)>0)
adld3p_na=as.numeric(is.na(adld3p))
adld3p_impute=adld3p
adld3p_impute[is.na(adld3p)]=mean(adld3p,na.rm=TRUE)
ca_yes=as.numeric(ca=='Yes')
ca_meta=as.numeric(ca=='Metastatic')
wt0=as.numeric(wtkilo1==0)
urin1_na=as.numeric(is.na(urin1))
urin1_impute=urin1
urin1_impute[is.na(urin1)]=mean(urin1,na.rm=TRUE)
NumComorbid=cardiohx+chfhx+dementhx+psychhx+chrpulhx+renalhx+liverhx+
  gibledhx+malighx+immunhx+transhx+amihx
Resp=as.numeric(resp=='Yes')
Card=as.numeric(card=='Yes')
Neuro=as.numeric(neuro=='Yes')
Gastr=as.numeric(gastr=='Yes')
Renal=as.numeric(renal=='Yes')
Meta=as.numeric(meta=='Yes')
Hema=as.numeric(hema=='Yes')
Seps=as.numeric(seps=='Yes')
Trauma=as.numeric(trauma=='Yes')
Ortho=as.numeric(ortho=='Yes')
Dnr1=as.numeric(dnr1=='Yes')
#estimate propensity score
pr<-glm(z~age+female+race_black+race_other+edu+income1+income2+income3+
          ins_care+ins_pcare+ins_caid+ins_no+ins_carecaid+
          cat1_copd+cat1_mosfsep+cat1_mosfmal+cat1_chf+cat1_coma+cat1_cirr+cat1_lung+cat1_colon+
          cat2_mosfsep+cat2_coma+cat2_mosfmal+cat2_lung+cat2_cirr+cat2_colon+
          Resp+Card+Neuro+Gastr+Renal+Meta+Hema+Seps+Trauma+Ortho+
          adld3p_impute+adld3p_na+das2d3pc+Dnr1+ca_yes+ca_meta+surv2md1+aps1+scoma1+
          wtkilo1+wt0+temp1+meanbp1+resp1+hrt1+pafi1+paco21+ph1+wblc1+hema1+
          sod1+pot1+crea1+bili1+alb1+urin1_impute+urin1_na+
          cardiohx+chfhx+dementhx+psychhx+chrpulhx+renalhx+liverhx+
          gibledhx+malighx+immunhx+transhx+amihx,family=binomial)$fitted.values
X<-cbind(aps1,surv2md1,age,NumComorbid,adld3p_impute,adld3p_na,das2d3pc,temp1,hrt1,meanbp1,
         resp1,wblc1,pafi1,paco21,ph1,crea1,alb1,scoma1,
         cat1_copd,cat1_mosfsep,cat1_mosfmal,cat1_chf,cat1_coma,cat1_cirr,cat1_lung,cat1_colon)
XX=cbind(X,pr)
Xfull<-cbind(age,female,race_black,race_other,edu,income1,income2,income3,
             ins_care,ins_pcare,ins_caid,ins_no,ins_carecaid,
             cat1_copd,cat1_mosfsep,cat1_mosfmal,cat1_chf,cat1_coma,cat1_cirr,cat1_lung,cat1_colon,
             cat2_mosfsep,cat2_coma,cat2_mosfmal,cat2_lung,cat2_cirr,cat2_colon,
             Resp,Card,Neuro,Gastr,Renal,Meta,Hema,Seps,Trauma,Ortho,
             adld3p_impute,adld3p_na,das2d3pc,Dnr1,ca_yes,ca_meta,surv2md1,aps1,scoma1,
             wtkilo1,wt0,temp1,meanbp1,resp1,hrt1,pafi1,paco21,ph1,wblc1,hema1,
             sod1,pot1,crea1,bili1,alb1,urin1_impute,urin1_na,
             cardiohx,chfhx,dementhx,psychhx,chrpulhx,renalhx,liverhx,
             gibledhx,malighx,immunhx,transhx,amihx)
XXf=cbind(Xfull,pr)
#create propensity score strata
pr_st4=as.integer(cut(pr,quantile(pr,c(0,.16,.5,.84,1)),include.lowest = TRUE))
detach(rhc)
data<-cbind(rhc,cbind(pr,pr_st4,z,female,race_black,race_other,income1,income2,income3,
                      ins_care,ins_pcare,ins_caid,ins_no,ins_carecaid,cat1_copd,cat1_mosfsep,
                      cat1_mosfmal,cat1_chf,cat1_coma,cat1_cirr,cat1_lung,cat1_colon,
                      cat2_mosfsep,cat2_coma,cat2_mosfmal,cat2_lung,cat2_cirr,cat2_colon,
                      adld3p_na,adld3p_impute,ca_yes,ca_meta,wt0,urin1_na,urin1_impute,NumComorbid,
                      Resp,Card,Neuro,Gastr,Renal,Meta,Hema,Seps,Trauma,Ortho,Dnr1))
rm(pr,pr_st4,z,female,race_black,race_other,income1,income2,income3,
   ins_care,ins_pcare,ins_caid,ins_no,ins_carecaid,cat1_copd,cat1_mosfsep,
   cat1_mosfmal,cat1_chf,cat1_coma,cat1_cirr,cat1_lung,cat1_colon,
   cat2_mosfsep,cat2_coma,cat2_mosfmal,cat2_lung,cat2_cirr,cat2_colon,
   adld3p_na,adld3p_impute,ca_yes,ca_meta,wt0,urin1_na,urin1_impute,NumComorbid,
   Resp,Card,Neuro,Gastr,Renal,Meta,Hema,Seps,Trauma,Ortho,Dnr1)
data$id=numeric(nrow(data))
data$id[which(data$z==1)]=1:sum(data$z)
data$id[which(data$z==0)]=(sum(data$z)+1):nrow(data)
data$cat2[is.na(data$cat2)]='Missing'

##match in Connors et al (1996)
data1=data[data$z==1,]
data0=data[data$z==0,]
pdm=abs(outer(data1$pr,data0$pr,"-"))
edm1=outer(data1$cat1,data0$cat1,"==")
edm2=outer(data1$cat2,data0$cat2,"==")
set.seed(1)
o=sample(1:sum(data$z),sum(data$z))
matched_or_not=rep(F,sum(1-data$z))
#control_num=numeric(sum(data$z))
greedy=c()
for (i in 1:sum(data$z)){
  ix=o[i]
  select=(pdm[ix,]<=0.03) & (edm1[ix,]==T) & (!matched_or_not)
  if (any(select)){
    cix=which(pdm[ix,]==min(pdm[ix,select]))[1]
    greedy=rbind(greedy,data1[ix,],data0[cix,])
    matched_or_not[cix]=T
  }
}

npairs=nrow(greedy)/2
npairs #1472
greedy$mset=rep(1:npairs,each=2)

Xmf_greedy=greedy[,c('age','female','race_black','race_other','edu','income1','income2','income3',
               'ins_care','ins_pcare','ins_caid','ins_no','ins_carecaid',
               'cat1_copd','cat1_mosfsep','cat1_mosfmal','cat1_chf','cat1_coma','cat1_cirr','cat1_lung','cat1_colon',
               'cat2_mosfsep','cat2_coma','cat2_mosfmal','cat2_lung','cat2_cirr','cat2_colon',
               'Resp','Card','Neuro','Gastr','Renal','Meta','Hema','Seps','Trauma','Ortho',
               'adld3p_impute','adld3p_na','das2d3pc','Dnr1','ca_yes','ca_meta','surv2md1','aps1','scoma1',
               'wtkilo1','wt0','temp1','meanbp1','resp1','hrt1','pafi1','paco21','ph1','wblc1','hema1',
               'sod1','pot1','crea1','bili1','alb1','urin1_impute','urin1_na',
               'cardiohx','chfhx','dementhx','psychhx','chrpulhx','renalhx','liverhx',
               'gibledhx','malighx','immunhx','transhx','amihx','pr')]
btb_greedy=balance(XXf,Xmf_greedy,data$z,greedy$z)
btb_greedy

##r best matching
fbv=as.factor(data$pr_st4):as.factor(data$cat1)
dm=mahad(data$z,X)
dist_quantile=c(0.01,0.1,0.25,0.5,1)
caliper_cutoff=c(0.01,0.02,0.03,0.03,1)
dist_cutoff=quantile(dm$d,dist_quantile)
dist_cutoff[length(dist_cutoff)]=Inf
result=rbestmatch(z=data$z,data=data,Xvar=Xfull,Pvar=X,score=data$pr,fine=fbv,
                  calipers=caliper_cutoff,cutdist=dist_cutoff,
                  nranks=4,ncontrol=1,rank=TRUE,s.cost=10,rpenalty=50)
result$nmatchedpairs
rbest=result$data
Xmf_rbest=rbest[,c('age','female','race_black','race_other','edu','income1','income2','income3',
             'ins_care','ins_pcare','ins_caid','ins_no','ins_carecaid','cat1_copd',
             'cat1_mosfsep','cat1_mosfmal','cat1_chf','cat1_coma','cat1_cirr','cat1_lung',
             'cat1_colon','cat2_mosfsep','cat2_coma','cat2_mosfmal','cat2_lung','cat2_cirr',
             'cat2_colon','Resp','Card','Neuro','Gastr','Renal','Meta','Hema','Seps',
             'Trauma','Ortho','adld3p_impute','adld3p_na','das2d3pc','Dnr1','ca_yes',
             'ca_meta','surv2md1','aps1','scoma1','wtkilo1','wt0','temp1','meanbp1','resp1',
             'hrt1','pafi1','paco21','ph1','wblc1','hema1','sod1','pot1','crea1','bili1',
             'alb1','urin1_impute','urin1_na','cardiohx','chfhx','dementhx','psychhx',
             'chrpulhx','renalhx','liverhx','gibledhx','malighx','immunhx','transhx','amihx',
             'pr')]
btb_rbest=balance(XXf,Xmf_rbest,data$z,rbest$z)
btb_rbest


##Table 1: 4 largest absolute SMD
tb1=rbind(sort(abs(btb_greedy[,6]),decreasing = TRUE)[1:4],
           sort(abs(btb_greedy[,5]),decreasing = TRUE)[1:4],
           sort(abs(btb_rbest[,5]),decreasing = TRUE)[1:4])
colnames(tb1)<-c(77,76,75,74)
rownames(tb1)<-c('Before','Greedy','R best')
tb1

##Table 2: fine balance
tbfb=table(rbest$z,rbest$pr_st4,rbest$cat1)
tb2=cbind(
  rbind(tbfb[2,,1],tbfb[1,,1],tbfb[2,,2],tbfb[1,,2],tbfb[2,,3],tbfb[1,,3],
        tbfb[2,,4],tbfb[1,,4],tbfb[2,,5],tbfb[1,,5],tbfb[2,,6],tbfb[1,,6],
        tbfb[2,,7],tbfb[1,,7],tbfb[2,,8],tbfb[1,,8],tbfb[2,,9],tbfb[1,,9]))
tb2

##Table 3: rbest match
grade=1:5
tb_match=c(as.numeric(table(result$rank_after)),0)
percent_match=tb_match/sum(tb_match)*100
tb_edge=as.numeric(table(result$rank_before))
tb_edge[5]=sum(data$z)*sum(1-data$z)-sum(tb_edge[1:4])
percent_edge=tb_edge/sum(tb_edge)*100
tb3=cbind(grade,caliper_cutoff,dist_quantile,dist_cutoff,
          tb_edge,percent_edge,tb_match,percent_match)
tb3

##Figure 1: propensity score
boxplot(
  data$pr[data$z==1],
  data$pr[data$z==0],
  greedy$pr[greedy$z==1],
  greedy$pr[greedy$z==0],
  rbest$pr[rbest$z==1],
  rbest$pr[rbest$z!=0],
  names=c("T","C","T","C","T","C"),
  ylim=c(0,1.15),cex.axis=.8,cex.lab=.9,
  ylab="Propensity Score",
  xlab="Three matching methods: G=greedy, O=optimal, R=rank maximal",
  main="Propensity Scores, Before and After Matching",cex.main=.9,
  axes=T,boxwex=rep(.4,8)
)

abline(v=c(2.5,4.5))
text(1,1.1,paste("n=2184"),cex=.7)
text(2,1.1,paste("n=3551"),cex=.7)
text(3.5,1.1,paste("G:",nrow(greedy)/2,"pairs"),cex=.7)
text(5.5,1.1,paste("R:",nrow(rbest)/2,"pairs"),cex=.7)
abline(h=1,lwd=2)

##Figure 2: distance
mdist_greedy=numeric(nrow(greedy)/2)
for (k in 1:length(mdist_greedy)){
  id1=greedy$id[2*k-1]
  id0=greedy$id[2*k]
  ix=which((dm$start==id1) & (dm$end==id0))
  mdist_greedy[k]=dm$d[ix]
}

mdist_rbest=numeric(nrow(rbest)/2)
for (k in 1:length(mdist_rbest)){
  id1=rbest$id[2*k-1]
  id0=rbest$id[2*k]
  ix=which((dm$start==id1) & (dm$end==id0))
  mdist_rbest[k]=dm$d[ix]
}

boxplot(
  mdist_greedy,
  mdist_rbest,
  names=c("G","R"),ylab="Distance",main="Within Pair Covariate Distances",
  xlab="Two matching methods: G=greedy, R=r-best",
  ylim=c(0,96),cex.lab=.85,cex.axes=.9

)

text(1,93,"1472 Pairs",cex=.75)
text(2,93,"1567 Pairs",cex=.75)
}

}

