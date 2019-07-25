#############################################
#############################################
###### 
###### Rosiglitazone data analysis
###### Author: Hwanhee Hong (hwanhee.hong@duke.edu)
###### Last update: July 18, 2019
###### Related paper: Meta-Analysis of Rare Adverse Events in Randomized Clinical Trials: Bayesian and Frequentist Methods (2019+)
######                by Hwanhee Hong, Chenguang Wang, and Gary L. Rosner
###### Data source: Nissen SE and Wolski K. 
######              Rosiglitazone revisited: an updated meta-analysis of risk for myocardial infarction and cardiovascular mortality.
######              Archives of internal medicine. 2010;170(14):1191-1201.
###### 
###### This code is to fit meta-analysis models (frequentist and Bayesian) to the rosiglitazone data.
###### This code provides 15 log odds ratio (LORs), 13 log relative risks (LRRs), and 17 risk differences estimates.
###### Bayesian models are fitted using rstan. Stan needs to be installed.
###### For Bayesian model comparison, we estimate Watanabe-Akaike or widely applicable information criterion (WAIC). 
###### 
###### This code produces 4 files:
######    (1) LOR, LRR, RD estimates (.csv)
######    (2) WAIC.txt
######    (3) Bayesian model results from Stan (.txt)
######    (4) Bayesian model diagnostic plots (.pdf)
######
#############################################
#############################################

###### Load packages
library(dplyr)
library(metafor)
library(data.table)  ## rbindlist
library(ggplot2)
library(gridExtra)
library(cowplot)  ## get_legend
library(loo) ## To calculate waic
#library(matrixStats)  ## fit$BUGSoutput$sims.list
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###### Set working directory
setwd("...")
source("AnalyzeData_Stan_Source.R")
source("Stan_models.R")

#### DATA should have the following variable names:
## Trial_ID: Actual trial ID
## n_rosi, n_ctl, y_rosi, y_ctl (e.g., y_rosi=MI_rosi or CVdeath_rosi)
#### outcome = MI or CV
#### dataset:
## 1=All studies (56 studies)
## 2=All studies except RECORD trial (55 studies)
## 3=Trials comparing rosiglitazone+X vs. X alone (42 studies)
## 4=Trials comparing rosiglitazone vs. placebo (13 studies)

DATA.raw <- read.table("RosiglitazoneData.csv",sep=",",header=T)
DATA1 <- DATA.raw
DATA2 <- DATA.raw[DATA.raw$Trial_ID!="RECORD trial",]
DATA3 <- DATA.raw[DATA.raw$Subgroup==1,]
DATA4 <- DATA.raw[DATA.raw$placebo==1,]

clean.data <- function(DATA){
  MI <- DATA[,c("Trial_ID","n_rosi","n_ctl","MI_rosi","MI_ctl")]
  CV <- DATA[,c("Trial_ID","n_rosi","n_ctl","CVdeath_rosi","CVdeath_ctl")]
  colnames(MI) <- colnames(CV) <- c("Trial_ID","n_rosi","n_ctl","y_rosi","y_ctl")
  
  output <- list(MI, CV)
  names(output) <- c("MI", "CV")
  return(output)
}

data1 <- clean.data(DATA1)
data2 <- clean.data(DATA2)
data3 <- clean.data(DATA3)
data4 <- clean.data(DATA4)

######
###### Run models
###### 
res.MI1 <- AnalyzeData(data1$MI, outcome="MI", dataset=1, niter=50000, nburnin=20000)
#res.CV1 <- AnalyzeData(data1$CV, outcome="CV", dataset=1, niter=50000, nburnin=20000)

#res.MI2 <- AnalyzeData(data2$MI, outcome="MI", dataset=2, niter=50000, nburnin=20000)
#res.CV2 <- AnalyzeData(data2$CV, outcome="CV", dataset=2, niter=50000, nburnin=20000)

#res.MI3 <- AnalyzeData(data3$MI, outcome="MI", dataset=3, niter=50000, nburnin=20000)
#res.CV3 <- AnalyzeData(data3$CV, outcome="CV", dataset=3, niter=50000, nburnin=20000)

#res.MI4 <- AnalyzeData(data4$MI, outcome="MI", dataset=4, niter=50000, nburnin=20000)
#res.CV4 <- AnalyzeData(data4$CV, outcome="CV", dataset=4, niter=50000, nburnin=20000)

######
###### WAIC Table
######

sink("WAIC.txt", append=FALSE)
cat("\n","\n","#### Outcome=MI; Data1; Logit ####","\n","\n")
format(res.MI1$LOR[res.MI1$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data2; Logit ####","\n","\n")
#format(res.MI2$LOR[res.MI2$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data3; Logit ####","\n","\n")
#format(res.MI3$LOR[res.MI3$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data4; Logit ####","\n","\n")
#format(res.MI4$LOR[res.MI4$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n")
#cat("\n","\n","#### Outcome=CV; Data1; Logit ####","\n","\n")
#format(res.CV1$LOR[res.CV1$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data2; Logit ####","\n","\n")
#format(res.CV2$LOR[res.CV2$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data3; Logit ####","\n","\n")
#format(res.CV3$LOR[res.CV3$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data4; Logit ####","\n","\n")
#format(res.CV4$LOR[res.CV4$LOR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
cat("\n","\n")
cat("\n","\n","#### Outcome=MI; Data1; Log ####","\n","\n")
format(res.MI1$LRR[res.MI1$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data2; Log ####","\n","\n")
#format(res.MI2$LRR[res.MI2$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data3; Log ####","\n","\n")
#format(res.MI3$LRR[res.MI3$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=MI; Data4; Log ####","\n","\n")
#format(res.MI4$LRR[res.MI4$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n")
#cat("\n","\n","#### Outcome=CV; Data1; Log ####","\n","\n")
#format(res.CV1$LRR[res.CV1$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data2; Log ####","\n","\n")
#format(res.CV2$LRR[res.CV2$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data3; Log ####","\n","\n")
#format(res.CV3$LRR[res.CV3$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
#cat("\n","\n","#### Outcome=CV; Data4; Log ####","\n","\n")
#format(res.CV4$LRR[res.CV4$LRR$col.group=="Bayesian", c("elpd_waic", "p_waic", "waic")], nsmall=1)
sink()











