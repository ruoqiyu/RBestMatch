\name{rhc}
\alias{rhc}
\docType{data}
\title{
Right Heart Catheterization Dataset
}
\description{
This dataset was used in Connors
et al. (1996), JAMA: The effectiveness of RHC in the initial care of critically III patients.
276:889-897, The data and description are available from Harrell's Hmisc R package, using getHdata() function, obtaining information from https://hbiostat.org/data.
}
\usage{data("rhc")}
\format{
  A data frame with 5735 observations on the following 62 variables.
  \describe{
    \item{\code{cat1}}{Primary disease category}
    \item{\code{cat2}}{Secondary disease category}
    \item{\code{ca}}{Cancer}
    \item{\code{sadmdte}}{Study admission date}
    \item{\code{dschdte}}{Hospital discharge date}
    \item{\code{dthdte}}{Date of death}
    \item{\code{lstctdte}}{Date of last contact}
    \item{\code{death}}{Death at any time up to 180 Days}
    \item{\code{cardiohx}}{Indicator for acute MI, peripheral vascular disease, severe cardiovascular symptoms (NYHA-Class III), very severe cardiovascular symptoms (NYHA-Class IV)}
    \item{\code{chfhx}}{Indicator for congestive heart failure}
    \item{\code{dementhx}}{Indicator for dementia, stroke or cerebral infarct, parkinson’s disease}
    \item{\code{psychhx}}{Indicator for psychiatric history, active psychosis or severe depression}
    \item{\code{chrpulhx}}{Indicator for chronic pulmonary disease, severe pulmonary disease, very severe pulmonary disease}
    \item{\code{renalhx}}{Indicator for chronic renal disease, chronic hemodialysis or peritoneal dialysis}
    \item{\code{liverhx}}{Indicator for cirrhosis, hepatic failure}
    \item{\code{gibledhx}}{Indicator for upper GI bleeding}
    \item{\code{malighx}}{Indicator for solid tumor, metastatic disease, chronic leukemia/myeloma, acute leukemia, lymphoma}
    \item{\code{immunhx}}{Indicator for immunosupperssion, organ transplant, HIV positivity, diabetes mellitus without end organ damage, diabetes mellitus with end organ damage, connective tissue disease}
    \item{\code{transhx}}{Indicator for transfer (> 24 Hours) from another hospital}
    \item{\code{amihx}}{Indicator for definite myocardial infarction}
    \item{\code{age}}{Age}
    \item{\code{sex}}{Sex}
    \item{\code{edu}}{Years of education}
    \item{\code{surv2md1}}{Support model estimate of the prob. of surviving 2 months at day 1}
    \item{\code{das2d3pc}}{DASI (Duke Activity Status Index)}
    \item{\code{t3d30}}{min{dthdte - sadmdte, 30}}
    \item{\code{dth30}}{Death within 30 days of admission}
    \item{\code{aps1}}{APACHE score at day 1}
    \item{\code{scoma1}}{Glasgow Coma Score at day 1}
    \item{\code{meanbp1}}{Mean blood pressure at day 1}
    \item{\code{wblc1}}{White blood cell count at day 1}
    \item{\code{hrt1}}{Heart rate at day 1}
    \item{\code{resp1}}{Respiratory rate at day 1}
    \item{\code{temp1}}{Temperature (celcius) at day 1}
    \item{\code{pafi1}}{PaO2/FIO2 ratio at day 1}
    \item{\code{alb1}}{Albumin at day 1}
    \item{\code{hema1}}{Hematocrit at day 1}
    \item{\code{bili1}}{Bilirubin at day 1}
    \item{\code{crea1}}{Serum creatinineat at day 1}
    \item{\code{sod1}}{Serum sodium at day 1}
    \item{\code{pot1}}{Serum potassium at day 1}
    \item{\code{paco21}}{PaCo2 at day 1}
    \item{\code{ph1}}{Serum PH (arterial) at day 1}
    \item{\code{swang1}}{Whether or not a patient received a RHC (also called the Swan-Ganz catheter) onthe first day in which the patient qualified for the SUPPORT study}
    \item{\code{wtkilo1}}{Weight in kilograms}
    \item{\code{dnr1}}{DNR status on day 1}
    \item{\code{ninsclas}}{Medical insurance}
    \item{\code{resp}}{Categories of admission diagnosis: Respiratory Diagnosis}
    \item{\code{card}}{Categories of admission diagnosis: Cardiovascular Diagnosis}
    \item{\code{neuro}}{Categories of admission diagnosis: Neurological Diagnosis}
    \item{\code{gastr}}{Categories of admission diagnosis: Gastrointestinal Diagnosis}
    \item{\code{renal}}{Categories of admission diagnosis: Renal Diagnosis}
    \item{\code{meta}}{Categories of admission diagnosis: Metabolic Diagnosis}
    \item{\code{hema}}{Categories of admission diagnosis: Hematologic Diagnosis}
    \item{\code{seps}}{Categories of admission diagnosis: Sepsis Diagnosis}
    \item{\code{trauma}}{Categories of admission diagnosis: Trauma Diagnosis}
    \item{\code{ortho}}{Categories of admission diagnosis: Orthopedic Diagnosis}
    \item{\code{adld3p}}{Activities of daily living at day 3}
    \item{\code{urin1}}{Urine output at day 1}
    \item{\code{race}}{Race}
    \item{\code{income}}{Income}
    \item{\code{ptid}}{Patient ID}
  }
}
\details{
The following code constructed the data as used here.

library(Hmisc)
getHdata(rhc)
save(rhc,file='rhc.rda', compress='xz')

}
\source{
The data and description are available from Harrell's Hmisc R package, using getHdata() function, obtaining information from https://hbiostat.org/data.
}
\references{
Harrell Jr, F. E., & Dupont, M. C. (2021). The Hmisc package. R package version, 4.6-0.
}
\examples{
data(rhc)
summary(rhc)
}
\keyword{datasets}
