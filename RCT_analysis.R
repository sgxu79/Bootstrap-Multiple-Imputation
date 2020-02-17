# Analysis with RCT data
#
# Steven Xu, 2020
#
# Estimate hazard ratio (HR), incidence rate (IR) and incidence rate differ-
# ence (IRD) using trial data only. Parametric confidence intervals constru- 
# cted for HR and IR, and bootstrap confidence interval constructed for IRD.

library(haven)
library(tidyr)
library(tibble)
library(survival)
library(tictoc)
library(boot)
library(gridExtra)
library(Cairo)
library(survey)
library(parallel)

gen <- data.frame(read_sas("..."))

group = gen$study

group[which(group=="")] = "CNICS"

study = gen$study

num_study= factor(study,c("A5095","A5142","A5175","A5202"),labels=1:4)

#Reconstruct the dataset to include only variables of interest
gen_s <- data.frame(age = as.numeric(gen$AGE),cd4 = as.numeric(gen$bcd4),lrna = as.numeric(gen$lrna),
                    sex = factor(gen$SEX,1:2,labels = c("Male","Female")),
                    race = factor(gen$raceth2,c(1:3,9),c("White","Black","Hispanic","Other")),
                    idu = factor(gen$ivdrugb,0:1,c("No","Yes")), aids = factor(gen$hxaids,0:1,c("No","Yes")),
                    hep = factor(gen$hep,0:1,c("Negative/Missing/Indeterminate","Positive")),
                    dep = factor(gen$hxdeprx,0:1,c("No","Yes")),rand_efv = factor(gen$rand_efv,0:1,c("EFV-free","EFV-containing")),
                    s = factor(gen$S,0:1,c("No","Selected")),study = factor(gen$study,c("A5095","A5142","A5175","A5202")), z = gen$Zobs,
                    isuicyrd = as.numeric(gen$isuicyr_d),ln_suicyrd = gen$ln_suicyrd,isuicwk_d=as.numeric(gen$isuicwk_d),isuic_d=as.numeric(gen$isuic_d),
                    datsuic = as.numeric(gen$datsuic),datsuicw = as.numeric(gen$datsuicw), datsuicy=as.numeric(gen$datsuicy),
                    ln_datsuicy = gen$ln_datsuicy,efv_indicator = factor(gen$efv_indicator,0:1,c("EFV-free","EFV-containing")))

gen_s$ln_datsuicy[gen_s$datsuicy == 0 ] = log(1/365.25)
gen_s$ln_suicyrd[gen_s$isuicyrd == 0 ] = log(1/365.25)


#Unweighted analysis using RCTs
RCT.df = subset(gen_s,s=="Selected")

design.ps = svydesign(ids=~1,weights=~1,strata =~study,data=RCT.df)

poi_fit = svyglm(isuic_d~offset(ln_suicyrd)+rand_efv,design=design.ps,family=quasipoisson())

cph_fit = coxph(Surv(isuicwk_d,isuic_d) ~ rand_efv+strata(study),ties = "efron",data = RCT.df)

v.ir.t = matrix(c(1,1),nrow=1)

v.ir.c = matrix(c(1,0),nrow=1)

log_IR_t = (v.ir.t%*%poi_fit$coefficients)[1]

log_IR_t_se = sqrt(v.ir.t%*%poi_fit$cov.unscaled%*%t(v.ir.t))[1]

round(exp(log_IR_t)*1000, 1) # Incidence Rate, EFV-cont group
round(exp(c(log_IR_t-qnorm(0.975)*log_IR_t_se,log_IR_t+qnorm(0.975)*log_IR_t_se))*1000, 1)

log_IR_c = (v.ir.c%*%poi_fit$coefficients)[1]

log_IR_c_se = sqrt(v.ir.c%*%poi_fit$cov.unscaled%*%t(v.ir.c))[1]

round(exp(log_IR_c)*1000, 1) # Incidence rate, EFV-free group
round(exp(c(log_IR_c-qnorm(0.975)*log_IR_c_se,log_IR_c+qnorm(0.975)*log_IR_c_se))*1000, 1)

log_HR = cph_fit$coefficients #estimated log Hazard Ratio
log_HR_se = sqrt(cph_fit$var)[1]

round(exp(log_HR),1)
exp(c(log_HR-qnorm(0.975)*log_HR_se,log_HR+qnorm(0.975)*log_HR_se))

complete.analysis = function(data,indices){
  
  library(survey)
  
  current.df = data[indices,]
  
  design.ps = svydesign(ids=~1, weights=~1, strata=~study, data=subset(current.df, s == "Selected"))
  
  poi_fit = svyglm(isuic_d ~ rand_efv + offset(ln_suicyrd), design=design.ps, family=quasipoisson())

  v.ir.t = matrix(c(1,1),nrow=1)
  
  v.ir.c = matrix(c(1,0),nrow=1)
  
  log_IR_t = v.ir.t%*%poi_fit$coefficients
  
  log_IR_c = v.ir.c%*%poi_fit$coefficients
  
  result = exp(log_IR_t[1])-exp(log_IR_c[1])
  
  return(result)
}


set.seed(1234,kind = "L'Ecuyer-CMRG")

mycpu = detectCores()

tic("IRD bootstrap CI computing")

#Stratified bootstrap wrapper
boot.out <- boot(data=RCT.df,complete.analysis,R=1000,strata=na.omit(num_study),parallel = "snow",ncpus=mycpu) 

toc()

round((exp(log_IR_t)-exp(log_IR_c))*1000, 1)         #Estimated IRD
round(quantile(boot.out$t,c(0.025,0.975))*1000, 1)   #Bootstrap 95% CI


