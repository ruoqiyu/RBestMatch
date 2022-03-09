\name{mahad}
\alias{mahad}
\title{Creates a Mahalanobis distance or its robust version for matching.}
\usage{
mahad(z,X,p=rep(1,length(z)),caliper=1,stdev=FALSE,exact=NULL,nearexact=NULL,
penalty=100,rank=TRUE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{X}{ A matrix with length(z) rows giving the covariates to calculate the Mahalanobis distance. X should be of full column rank.}
\item{p}{ A vector of score for which caliper is applied for, e.g., p can be chosen as the propensity score. }
\item{caliper}{ A caliper applied to score p, requiring the treated and control pairs have score distance no more than caliper. The caliper is applied to remove edges from the network for matching, i.e., the reduced network only contains treated-control edges satisfying the caliper.}
\item{stdev}{ A boolean denoting whether the caliper is on pooled standard deviation scale. The default choice is FALSE so that the caliper is on the original scale of score p.}
\item{exact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{penalty}{ The penalty for a mismatch on nearexact.}
\item{rank}{ A boolean denoting whether the distance is the standard Mahalanobis distance or its robust version. The default choice is TRUE to calculate the rank-based Mahalanobis distance.}
}

\value{
\item{d}{ A distance object for each pair of treated and control.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
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
detach(rhc)

dm1=mahad(z,X,pr,0.1)
length(dm1$d)
\donttest{
dm2=mahad(z,X,pr)
length(dm2$d)
}
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
