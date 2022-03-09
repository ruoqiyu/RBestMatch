\name{rbestdist}
\alias{rbestdist}
\title{Creates a penalized distance list for r-best matching.}
\usage{
rbestdist(dist,r,penalty=NULL,beta=100)
}
\arguments{
\item{dist}{ A distance list used for matching. dist contains three elements, d: a distance object for each pair of treated and control; start: the treated subject for each distance; end: the control subject for each distance.}
\item{r}{ A vector of rank with length length(dist$d), denoting the preference for each treated-control pair considered in dist (1: best; 2: second best; etc).}
\item{penalty}{ A vector of penalty with length max(r) for different edge sets to conduct r-best matching. If NULL, use the tuning parameter beta to construct the penalty}.
\item{beta}{ If penalty is not specified, then it is calculated based on the tuning parameter beta.}
}

\value{
\item{newdist}{ A new distance object with d being the penalized distance and start and end same as the original distance object dist.}
\item{fbpenalty}{ The penalty for fine balance in r-best matching.}
}
\description{
Computes a penalized Mahalanobis distance (or its robust version) list for use in r-best matching.

}
\details{
The r-best matching problem is solved as a network flow problem with penalized distances, with preferance to edges that rank high. This function is used to construct the new distance object.
}
\examples{
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

dm0=mahad(z,X,pr,0.1)
pr1=pr[z==1]
pr0=pr[z==0]
prd=abs(pr1[dm0$start]-pr0[dm0$end-sum(z)])
calip1<-prd<=(.01)&(dm0$d<=quantile(dm0$d,0.01))
calip2<-prd<=(.02)&(dm0$d<=quantile(dm0$d,0.1))
calip3<-prd<=(.03)&(dm0$d<=quantile(dm0$d,0.25))
calip4<-prd<=(.03) &(dm0$d<=quantile(dm0$d,0.5))
r=rep(5,length(dm0$d))
r[calip4]=4
r[calip3]=3
r[calip2]=2
r[calip1]=1
which4=which(r<=4)
dm0$d=dm0$d[which4]
dm0$start=dm0$start[which4]
dm0$end=dm0$end[which4]
r4=r[which4]
dm=rbestdist(dm0,r4,beta=50)
head(dm$newdist$d)
dm$fbpenalty
detach(rhc)
}

