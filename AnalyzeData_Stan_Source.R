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
###### Source file associated with MetaAnalysis_Rosiglitazone.R
###### This source file includes user-written function to fit frequentist and Bayesian meta-analysis models.
######
#############################################
#############################################

###### Basic calculation functions
## pc = probability (risk) in control group
## pa = probability (risk) in active group
OR.cal <- function(pc, pa) {
  OR <- (pa/(1-pa)) / (pc/(1-pc))
  return(OR)
}

RR.cal <- function(pc, pa) {
  RR <- pa/pc
  return(RR)
}

RD.cal <- function(pc, pa) {
  RD <- pa-pc
  return(RD)
}

expit <- function(p) { exp(p)/(1+exp(p)) }
q.025<-function(x){quantile(x, probs = 0.025)}
q.975<-function(x){quantile(x, probs = 0.975)}

#### Save results for naive estimator
## For Data Analysis (DA)
output.pool.DA <- function(mod, measure, const, est, se){
  low <- est - const*se
  up <- est + const*se
  width <- up - low
  tau <- NA
  res <- c(est, se, low, up, width, tau, rep(NA,6))
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}
output.pool.RE.DA <- function(mod, measure, const, est, se, tau){
  low <- est - const*se
  up <- est + const*se
  width <- up - low
  tau <- tau
  res <- c(est, se, low, up, width, tau, rep(NA,6))
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

#### Save results for frequentist models (fixed effect)
## For Data Analysis (DA)
output.FE.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, NA, rep(NA,6))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,12)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}
#### Save results for frequentist models (random effect)
output.RE.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, sqrt(res$tau2), rep(NA,6))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,12)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}

#### Save results for Bayesian models with waic (fixed effect)
####
output.Bayes.cte.DA <- function(mod, measure, mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}
#### Save rsults for Bayesian models with waic (random effect)
output.Bayes.hte.DA <- function(mod, measure, mcmc, tau.mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- median(tau.mcmc)
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

#### Analyze NW data for CT revision
#### DATA should have the following variable names:
## Trial_ID: Actual trial ID
## n_rosi, n_ctl, y_rosi, y_ctl (e.g., y_rosi=MI_rosi or CVdeath_rosi)
#### outcome = MI or CV
#### dataset:
## 1=All studies (56 studies)
## 2=All studies except RECORD trial (55 studies)
## 3=Trials comparing rosiglitazone+X vs. X alone (42 studies)
## 4=Trials comparing rosiglitazone vs. placebo (13 studies)

AnalyzeData <- function(DATA, outcome, dataset, niter, nburnin){
  
  ###### Modify data structure
  data <- DATA %>% 
    rename(nc=n_ctl, na=n_rosi, ea=y_rosi, ec=y_ctl) %>%
    mutate(study = 1:n(),
           ec_cc1 = ifelse(ec==0 | ea==0, ec+0.5, ec),
           ea_cc1 = ifelse(ec==0 | ea==0, ea+0.5, ea),
           nc_cc1 = ifelse(ec==0 | ea==0, nc+1, nc),
           na_cc1 = ifelse(ec==0 | ea==0, na+1, na))
  data.l <- data.frame(rbindlist( list(cbind(data[,c("Trial_ID","study","nc","ec")],1), cbind(data[,c("Trial_ID","study","na","ea")],2))))
  colnames(data.l) <- c("Trial_ID","study","n","y","t")
  
  ###### Frequentist models ######
  ###### 1. Naive
  a <- sum(data$ec)
  b <- sum(data$nc) - a
  c <- sum(data$ea)
  d <- sum(data$na) - c
  #### 1a. LOR
  lor.e.naive <- log(OR.cal( pc=a/(a+b), pa=c/(c+d) ))
  lor.se.naive <- sqrt( 1/a + 1/b + 1/c + 1/d )
  lor.res.naive <- output.pool.DA(mod="naive", measure="lor", const=1.96, est=lor.e.naive, se=lor.se.naive)
  
  #### 1b. LRR
  lrr.e.naive <- log(RR.cal( pc=a/(a+b), pa=c/(c+d) ))
  lrr.se.naive <- sqrt( 1/c + 1/a - 1/(a+b) - 1/(c+d) )
  lrr.res.naive <- output.pool.DA(mod="naive", measure="lrr", const=1.96, est=lrr.e.naive, se=lrr.se.naive)
  
  #### 1c. RD
  rd.e.naive <- c/(c+d) - a/(a+b)
  rd.se.naive <- sqrt( a*b/(a+b)^3 + c*d/(c+d)^3 )
  rd.res.naive <- output.pool.DA(mod="naive", measure="rd", const=1.96, est=rd.e.naive, se=rd.se.naive)
  
  ###### 2. Peto
  #### 2a. LOR
  lor.res.peto <- output.FE.DA(mod="peto",measure="lor",
                               res=try(rma.peto(ai=ea, ci=ec, n1i=na, n2i=nc, data=data)))
  
  ###### 3. MH (no data modification)
  #### 3a. LOR
  lor.res.mh <- output.FE.DA(mod="mh",measure="lor",
                             res=try(rma.mh(ai=ea, ci=ec, n1i=na, n2i=nc, measure="OR",data=data,correct=F)))
  #### 3b. LRR
  lrr.res.mh <- output.FE.DA(mod="mh",measure="lrr",
                             res=try(rma.mh(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RR",data=data,correct=F)))
  #### 3c. RD
  rd.res.mh <- output.FE.DA(mod="mh",measure="rd",
                            res=try(rma.mh(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RD",data=data,correct=F)))
  
  ###### 4. MH (data modification; add 0.5 to zero cells)
  #### 4a. LOR
  lor.res.mhdm <- output.FE.DA(mod="mhdm",measure="lor",
                               res=try(rma.mh(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="OR",data=data,correct=F)))
  #### 4b. LRR
  lrr.res.mhdm <- output.FE.DA(mod="mhdm",measure="lrr",
                               res=try(rma.mh(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="RR",data=data,correct=F)))
  #### 4c. RD
  rd.res.mhdm <- output.FE.DA(mod="mhdm",measure="rd",
                              res=try(rma.mh(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="RD",data=data,correct=F)))
  
  ###### 5. IV, fixed effect (data modification has been implemented in the code; add 0.5 to zero cells)
  #### 5a. LOR
  lor.res.ivfe <- output.FE.DA(mod="ivfe",measure="lor",
                               res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="OR",data=data,method="FE")))
  #### 5b. LRR
  lrr.res.ivfe <- output.FE.DA(mod="ivfe",measure="lrr",
                               res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RR",data=data,method="FE")))
  #### 5c. RD
  rd.res.ivfe <- output.FE.DA(mod="ivfe",measure="rd",
                              res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RD",data=data,method="FE")))
  
  ###### 6. IV, random effect DL (data modification has been implemented in the code; add 0.5 to zero cells)
  #### 6a. LOR
  lor.res.ivre <- output.RE.DA(mod="ivre",measure="lor",
                               res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="OR",data=data,method="DL")))
  #### 6b. LRR
  lrr.res.ivre <- output.RE.DA(mod="ivre",measure="lrr",
                               res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RR",data=data,method="DL")))
  #### 6c. RD
  rd.res.ivre <- output.RE.DA(mod="ivre",measure="rd",
                              res=try(rma(ai=ea, ci=ec, n1i=na, n2i=nc, measure="RD",data=data,method="DL")))
  
  ###### 7. Shuster unweighted (Shuster et al. 2012 RSM)
  #### 1 = control; 2= active; p1j=risk in control in study j (hat); p1=mean of p1j
  M <- nrow(data)
  p1j <- data$ec/data$nc
  p2j <- data$ea/data$na
  p1 <- sum(p1j)/M
  p2 <- sum(p2j)/M
  c11 <- sum((p1j-p1)^2)/(M-1)
  c22 <- sum((p2j-p2)^2)/(M-1)
  c12 <- sum((p1j-p1)*(p2j-p2))/(M-1)
  
  #### 7a. LOR
  Q <- c11/((p1*(1-p1))^2)
  R <- c22/((p2*(1-p2))^2)
  S <- c12/(p1*p2*(1-p1)*(1-p2))
  tval <- qt(0.975, M-1)
  lor.e.sgsunwgt <- log(OR.cal( pc=p1, pa=p2 ))
  lor.se.sgsunwgt <- sqrt( (Q+R-2*S)/M )
  lor.res.sgsunwgt <- output.pool.DA(mod="sgsunwgt", measure="lor", const=tval, est=lor.e.sgsunwgt, se=lor.se.sgsunwgt)
  
  #### 7b. LRR
  Q <- c11/(p1^2)
  R <- c22/(p2^2)
  S <- c12/(p1*p2)
  tval <- qt(0.975, M-1)
  lrr.e.sgsunwgt <- log(RR.cal( pc=p1, pa=p2 ))
  lrr.se.sgsunwgt <- sqrt( (Q+R-2*S)/M )
  lrr.res.sgsunwgt <- output.pool.DA(mod="sgsunwgt", measure="lrr", const=tval, est=lrr.e.sgsunwgt, se=lrr.se.sgsunwgt)
  
  #### 7c. RD
  Q <- c11
  R <- c22
  S <- c12
  tval <- qt(0.975, M-1)
  rd.e.sgsunwgt <- p2-p1
  rd.se.sgsunwgt <- sqrt( (Q+R-2*S)/M )
  rd.res.sgsunwgt <- output.pool.DA(mod="sgsunwgt", measure="rd", const=tval, est=rd.e.sgsunwgt, se=rd.se.sgsunwgt)
  
  ###### 8. Shuster weighted (Shuster et al. 2012 RSM)
  M <- nrow(data)
  p1j <- data$ec/data$nc
  p2j <- data$ea/data$na
  uj <- (data$nc + data$na)/2
  u <- mean(uj)
  
  A1j <- p1j*uj; B1j <- (1-p1j)*uj
  A2j <- p2j*uj; B2j <- (1-p2j)*uj
  A1 <- mean(A1j); B1 <- mean(B1j)
  A2 <- mean(A2j); B2 <- mean(B2j)
  
  sA1 <-sd(A1j); sB1 <- sd(B1j)
  sA2 <-sd(A2j); sB2 <- sd(B2j)
  cA1B1 <- cov(A1j, B1j)
  cA1B2 <- cov(A1j, B2j)
  cA1A2 <- cov(A1j, A2j)
  cA2B1 <- cov(A2j, B1j)
  cA2B2 <- cov(A2j, B2j)
  cB1B2 <- cov(B1j, B2j)
  
  #### 8a. LOR
  lor.e.sgswgt <- log((A2*B1)/(A1*B2))
  lor.var.sgswgt <- ( (sA1/A1)^2 + (sA2/A2)^2 + (sB1/B1)^2 + (sB2/B2)^2 + 2*(cA1B2/(A1*B2)) + 2*(cA2B1/(A2*B1)) - 2*(cA1A2/(A1*A2)) - 2*(cB1B2/(B1*B2)) - 2*(cA1B1/(A1*B1)) - 2*(cA2B2/(A2*B2)) )/M
  tval <- qt(0.975, M-2)
  lor.res.sgswgt <- output.pool.DA(mod="sgswgt", measure="lor", const=tval, est=lor.e.sgswgt, se=sqrt(lor.var.sgswgt))
  
  #### 8b. LRR
  lrr.e.sgswgt <- log(A2/A1)
  lrr.var.sgswgt <- ( (sA1/A1)^2 + (sA2/A2)^2 - 2*(cA1A2/(A1*A2)) )/M
  tval <- qt(0.975, M-2)
  lrr.res.sgswgt <- output.pool.DA(mod="sgswgt", measure="lrr", const=tval, est=lrr.e.sgswgt, se=sqrt(lrr.var.sgswgt))
  
  #### 8c. RD
  Dj <- A2j - A1j
  D <- mean(Dj)
  tval <- qt(0.975, M-2)
  rd.e.sgswgt <- D/u
  rd.var.sgswgt <- ( (var(Dj)/u^2) + (var(uj)*D^2/u^4) - 2*(D*cov(Dj, uj)/u^3))/M
  rd.res.sgswgt <- output.pool.DA(mod="sgswgt", measure="rd", const=tval, est=rd.e.sgswgt, se=sqrt(rd.var.sgswgt))
  
  ###### 9. Simple average (Bhaumik et al. 2012 JASA) with tau_DSL (add 0.5 to every cells!!)
  #### 9a. LOR
  k <- nrow(data)
  pcihat <- (data$ec+0.5)/(data$nc+1)
  paihat <- (data$ea+0.5)/(data$na+1)
  lori <- log(OR.cal(pcihat, paihat))
  lor.e.sa <- mean(lori)
  ## Calculate tau_DSL
  var0i <- 1/(data$na*paihat*(1-paihat)) + 1/(data$nc*pcihat*(1-pcihat))
  w0i <- 1/var0i
  sa.w0 <- sum(w0i*lori)/sum(w0i)
  Q <- sum(w0i*(lori-sa.w0)^2)
  var.dsl <- (Q-(k-1))/(sum(w0i) - (sum(w0i^2))/(sum(w0i)) )
  tau.dsl <- ifelse(var.dsl>0, sqrt(var.dsl), 0)
  lor.vari <- var0i + tau.dsl
  lor.se.sa <- sqrt(sum(lor.vari)/k^2)
  tval <- qt(0.975, k-1)
  lor.res.sa <- output.pool.RE.DA(mod="sa", measure="lor", const=tval, est=lor.e.sa, se=lor.se.sa, tau=tau.dsl)
  
  ###### Bayesian models using Stan ######
  ###### Bayesian models don't use data modification ######
  #niter <- 50000
  #nburnin <- 20000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  s <- data$study 
  ec <- data$ec 
  ea <- data$ea 
  nc <- data$nc 
  na <- data$na
  
  ###### 1. CTE-vague
  #### 1a. LOR & RD.logit
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lor", "rd", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=ctelogit.stan, data=stan.data, pars=stan.params, seed=752346,
                  chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: ctevaguelogit ####","\n","\n")
  print(stan.fit, par=c("lor", "rd"), digits=4)
  sink()
  pdf(file=paste("diag_ctevaguelogit_",outcome,dataset,".pdf",sep=""), onefile=TRUE)
  pairs(stan.fit, pars=c("lor", "lp__"))
  print(traceplot(stan.fit, pars=c("lor", "rd"), nrow=2, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  rd.mcmc <- fit.mcmc$rd
  #loglik.mcmc <- extract_log_lik(stan.fit)
  loglik.mcmc <- fit.mcmc$log_lik
  lor.res.ctevaguelogit <- output.Bayes.cte.DA(mod="ctevaguelogit",measure="lor",mcmc=lor.mcmc, loglik=loglik.mcmc)
  rd.res.ctevaguelogit <- output.Bayes.cte.DA(mod="ctevaguelogit",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  #### 1b. LRR & RD.log
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lrr", "rd", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=ctelog.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: ctevaguelog ####","\n","\n")
  print(stan.fit, par=c("lrr", "rd"), digits=4)
  sink()
  pdf(file=paste("diag_ctevaguelog_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("lrr", "lp__"))
  print(traceplot(stan.fit, pars=c("lrr", "rd"), nrow=2, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  loglik.mcmc <- fit.mcmc$log_lik
  lrr.res.ctevaguelog <- output.Bayes.cte.DA(mod="ctevaguelog",measure="lrr", mcmc=lrr.mcmc, loglik=loglik.mcmc)
  rd.res.ctevaguelog <- output.Bayes.cte.DA(mod="ctevaguelog",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  ###### 2. HTE-vague
  #### 2a. LOR & RD.logit
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lor", "rd", "tau","log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=htelogit.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: htevaguelogit ####","\n","\n")
  print(stan.fit, par=c("lor", "rd", "tau"), digits=4)
  sink()
  pdf(file=paste("diag_htevaguelogit_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("lor", "tau", "lp__"))
  print(traceplot(stan.fit, pars=c("lor", "rd", "tau"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  rd.mcmc <- fit.mcmc$rd
  tau.mcmc <- fit.mcmc$tau
  loglik.mcmc <- fit.mcmc$log_lik
  lor.res.htevaguelogit <- output.Bayes.hte.DA(mod="htevaguelogit",measure="lor", mcmc=lor.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.htevaguelogit <- output.Bayes.cte.DA(mod="htevaguelogit",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  #### 2b. LRR & RD.log
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lrr", "rd", "tau","log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=htelog.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: htevaguelog ####","\n","\n")
  print(stan.fit, par=c("lrr", "rd", "tau"), digits=4)
  sink()
  pdf(file=paste("diag_htevaguelog_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("lrr", "tau", "lp__"))
  print(traceplot(stan.fit, pars=c("lrr", "rd", "tau"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  tau.mcmc <- fit.mcmc$tau
  loglik.mcmc <- fit.mcmc$log_lik
  lrr.res.htevaguelog <- output.Bayes.hte.DA(mod="htevaguelog",measure="lrr", mcmc=lrr.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.htevaguelog <- output.Bayes.cte.DA(mod="htevaguelog",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)

  ###### 3. HTE-shirinkage
  #### 3a. LOR & RD.logit
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lor", "rd", "tau", "mP", "tauP", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=hteshlogit.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: hteshlogit ####","\n","\n")
  print(stan.fit, par=c("lor", "rd", "tau", "mP", "tauP"), digits=4)
  sink()
  pdf(file=paste("diag_hteshlogit_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("lor", "tau", "mP", "tauP", "lp__"))
  print(traceplot(stan.fit, pars=c("lor", "rd", "tau", "mP", "tauP"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  rd.mcmc <- fit.mcmc$rd
  tau.mcmc <- fit.mcmc$tau
  loglik.mcmc <- fit.mcmc$log_lik
  lor.res.hteshlogit <- output.Bayes.hte.DA(mod="hteshlogit",measure="lor", mcmc=lor.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.hteshlogit <- output.Bayes.cte.DA(mod="hteshlogit",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)

  #### 3b. LRR & RD.log
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lrr", "rd", "tau", "mP", "tauP", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=hteshlog.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: hteshlog ####","\n","\n")
  print(stan.fit, par=c("lrr", "rd", "tau", "mP", "tauP"), digits=4)
  sink()
  pdf(file=paste("diag_hteshlog_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("lrr", "tau", "mP", "tauP", "lp__"))
  print(traceplot(stan.fit, pars=c("lrr", "rd", "tau", "mP", "tauP"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  tau.mcmc <- fit.mcmc$tau
  loglik.mcmc <- fit.mcmc$log_lik
  lrr.res.hteshlog <- output.Bayes.hte.DA(mod="hteshlog",measure="lrr", mcmc=lrr.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.hteshlog <- output.Bayes.cte.DA(mod="hteshlog",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  ###### 4. HTE-AB
  #### 4a. LOR & RD.logit
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na, "NT"=2, "mean0"=c(0,0), "nu"=2, "L_inv"=diag(sqrt(4),2))
  stan.params <- c("lor", "rd", "phat_c", "phat_a", "logit_risk_c", "logit_risk_a", "tau_c", "tau_a", "A_inv_L_inv", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=ablogit.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: ablogit ####","\n","\n")
  print(stan.fit, par=c("lor", "rd", "phat_c", "phat_a", "logit_risk_c", "logit_risk_a", "tau_c", "tau_a", "A_inv_L_inv"), digits=4)
  sink()
  pdf(file=paste("diag_ablogit_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("logit_risk_c", "logit_risk_a", "tau_c", "tau_a", "lp__"))
  print(traceplot(stan.fit, pars=c("lor", "rd", "phat_c", "phat_a", "tau_c", "tau_a"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  rd.mcmc <- fit.mcmc$rd
  loglik.mcmc <- fit.mcmc$log_lik
  # AB models' taus are not for variability of lor, so we don't record them.
  lor.res.ablogit<- output.Bayes.cte.DA(mod="ablogit",measure="lor", mcmc=lor.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.ablogit <- output.Bayes.cte.DA(mod="ablogit",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  #### 4b. LRR & RD.log
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na, "NT"=2, "mean0"=c(0,0), "nu"=2, "L_inv"=diag(sqrt(4),2))
  stan.params <- c("lrr", "rd", "phat_c", "phat_a", "log_risk_c", "log_risk_a", "tau_c", "tau_a", "A_inv_L_inv", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=ablog.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: ablog ####","\n","\n")
  print(stan.fit, par=c("lrr", "rd", "phat_c", "phat_a", "tau_c", "tau_a", "log_risk_c", "log_risk_a", "A_inv_L_inv"), digits=4)
  sink()
  pdf(file=paste("diag_ablog_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("log_risk_c", "log_risk_a", "tau_c", "tau_a", "lp__"))
  print(traceplot(stan.fit, pars=c("lrr", "rd", "phat_c", "phat_a", "tau_c", "tau_a"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  loglik.mcmc <- fit.mcmc$log_lik
  # AB models' taus are not for variability of lor, so we don't record them.
  lrr.res.ablog<- output.Bayes.cte.DA(mod="ablog",measure="lrr", mcmc=lrr.mcmc, loglik=loglik.mcmc)
  # we don't estimate tau for RD
  rd.res.ablog <- output.Bayes.cte.DA(mod="ablog",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  ###### 5. CTE-Beta
  #### 5a. LOR, LRR, RD
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("pc", "pa", "lor", "lrr", "rd", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=ctebeta.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: ctebeta ####","\n","\n")
  print(stan.fit, par=c("pc", "pa", "lor", "lrr", "rd"), digits=4)
  sink()
  pdf(file=paste("diag_ctebeta_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("pc", "pa", "lp__"))
  print(traceplot(stan.fit, pars=c("pc", "pa", "lor", "lrr", "rd"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  loglik.mcmc <- fit.mcmc$log_lik
  # Beta prior models' taus are not for variability of lor, so we don't record them.
  lor.res.ctebeta<- output.Bayes.cte.DA(mod="ctebeta",measure="lor", mcmc=lor.mcmc, loglik=loglik.mcmc)
  lrr.res.ctebeta<- output.Bayes.cte.DA(mod="ctebeta",measure="lrr", mcmc=lrr.mcmc, loglik=loglik.mcmc)
  rd.res.ctebeta <- output.Bayes.cte.DA(mod="ctebeta",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)

  ###### 6. HTE-Beta
  #### 6a. LOR, LRR, RD
  ## sensitive to the prior for Vc and Va!!
  stan.data <- list("NS"=NS, "s"=s, "ec"=ec, "ea"=ea, "nc"=nc, "na"=na)
  stan.params <- c("lor", "lrr", "rd", "Uc", "Ua", "Vc", "Va", "log_lik")
  stan.fit <- NULL
  stan.fit <- stan(model_code=htebeta.stan, data=stan.data, pars=stan.params, seed=752346,
                   chains=nchains, iter=niter, warmup=nburnin, control = list(adapt_delta = 0.99))
  sink(paste("Print_stan_fit_",outcome,dataset,".txt",sep=""), append=TRUE)
  cat("\n","\n","#### Model: htebeta ####","\n","\n")
  print(stan.fit, par=c("lor", "lrr", "rd", "Uc", "Ua", "Vc", "Va"), digits=4)
  sink()
  pdf(file=paste("diag_htebeta_",outcome,dataset,".pdf",sep=""))
  pairs(stan.fit, pars=c("Uc", "Ua", "Vc", "Va", "lp__"))
  print(traceplot(stan.fit, pars=c("lor", "lrr", "rd", "Uc", "Ua", "Vc", "Va"), nrow=3, inc_warmup=TRUE))
  dev.off()
  fit.mcmc <- extract(stan.fit, permuted=TRUE)
  lor.mcmc <- fit.mcmc$lor
  lrr.mcmc <- fit.mcmc$lrr
  rd.mcmc <- fit.mcmc$rd
  loglik.mcmc <- fit.mcmc$log_lik
  # Beta prior models' taus are not for variability of lor, so we don't record them.
  lor.res.htebeta<- output.Bayes.cte.DA(mod="htebeta",measure="lor", mcmc=lor.mcmc, loglik=loglik.mcmc)
  lrr.res.htebeta<- output.Bayes.cte.DA(mod="htebeta",measure="lrr", mcmc=lrr.mcmc, loglik=loglik.mcmc)
  rd.res.htebeta <- output.Bayes.cte.DA(mod="htebeta",measure="rd", mcmc=rd.mcmc, loglik=loglik.mcmc)
  
  ######
  ###### Combine results
  #### 1. LOR
  lor <- data.frame(rbind(lor.res.naive, lor.res.peto, lor.res.mh, lor.res.mhdm, lor.res.ivfe, 
                          lor.res.sgsunwgt, lor.res.sgswgt, lor.res.ivre, lor.res.sa, lor.res.ctevaguelogit,
                          lor.res.ctebeta, lor.res.htevaguelogit, lor.res.hteshlogit, lor.res.ablogit, lor.res.htebeta))
  lor$modelname <- c("naive", "peto", "mh", "mhdm", "ivfe",
                     "sgsunwgt", "sgswgt", "ivre", "sa", "ctevaguelogit",
                     "ctebeta", "htevaguelogit", "hteshlogit", "ablogit", "htebeta")
  lor$col.group <- factor(c(rep(1,9), rep(2,6))); levels(lor$col.group) <- c("Frequentist", "Bayesian")
  lor$pch.group <- factor(c(1,1,1,1,1, 1,1,2,2,1, 1,2,2,2,2)); levels(lor$pch.group) <- c("CTE model", "HTE model")
  lor$outcome <- outcome
  lor$dataset <- dataset
  
  #### 2. LRR
  lrr <- data.frame(rbind(lrr.res.naive, lrr.res.mh, lrr.res.mhdm, lrr.res.ivfe, lrr.res.sgsunwgt,
                          lrr.res.sgswgt, lrr.res.ivre, lrr.res.ctevaguelog, lrr.res.ctebeta, lrr.res.htevaguelog,
                          lrr.res.hteshlog, lrr.res.ablog, lrr.res.htebeta))
  lrr$modelname <- c("naive", "mh", "mhdm", "ivfe", "sgsunwgt",
                     "sgswgt", "ivre", "ctevaguelog", "ctebeta", "htevaguelog",
                     "hteshlog", "ablog", "htebeta")
  lrr$col.group <- factor(c(rep(1,7), rep(2,6))); levels(lrr$col.group) <- c("Frequentist", "Bayesian")
  lrr$pch.group <- factor(c(1,1,1,1,1, 1,2,1,1,2, 2,2,2)); levels(lrr$pch.group) <- c("CTE model", "HTE model")
  lrr$outcome <- outcome
  lrr$dataset <- dataset
  
  #### 3. RD
  rd <- data.frame(rbind(rd.res.naive, rd.res.mh, rd.res.mhdm, rd.res.ivfe, rd.res.sgsunwgt,
                         rd.res.sgswgt, rd.res.ivre, rd.res.ctevaguelogit, rd.res.ctevaguelog, rd.res.ctebeta,
                         rd.res.htevaguelogit, rd.res.htevaguelog, rd.res.hteshlogit, rd.res.hteshlog, rd.res.ablogit,
                         rd.res.ablog, rd.res.htebeta))
  rd$modelname <- c("naive", "mh", "mhdm", "ivfe", "sgsunwgt",
                    "sgswgt", "ivre", "ctevaguelogit", "ctevaguelog", "ctebeta",
                    "htevaguelogit", "htevaguelog", "hteshlogit", "hteshlog", "ablogit",
                    "ablog", "htebeta")
  rd$col.group <- factor(c(rep(1,7), rep(2,10))); levels(rd$col.group) <- c("Frequentist", "Bayesian")
  rd$pch.group <- factor(c(1,1,1,1,1, 1,2,1,1,1, 2,2,2,2,2, 2,2)); levels(rd$pch.group) <- c("CTE model", "HTE model")
  rd$outcome <- outcome
  rd$dataset <- dataset
  
  colnames(lor) <- colnames(lrr) <- colnames(rd) <- c("est","se","low","up","width","tau",
                                                      "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic",
                                                      "modelname","col.group","pch.group","outcome","dataset")
  
  write.csv(lor, file=paste(outcome,dataset,"_lor.csv",sep=""))
  write.csv(lrr, file=paste(outcome,dataset,"_lrr.csv",sep=""))
  write.csv(rd, file=paste(outcome,dataset,"_rd.csv",sep=""))
  
  output <- list(lor, lrr, rd)
  names(output) <- c("LOR", "LRR", "RD")
  return(output)
}

