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
###### This source file includes Stan model code to fit Bayesian meta-analysis methods.
######
#############################################
#############################################

ctelogit.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lor;
  vector[NS] mu;
}
transformed parameters{
  vector[NS] logit_pc;
  vector[NS] logit_pa;
  vector[NS] pc;
  vector[NS] pa;
  //
  for(n in 1:NS) {
    logit_pc[n] = mu[n];
    logit_pa[n] = mu[n] + lor;
    //
    pc[n] = inv_logit(logit_pc[n]);
    pa[n] = inv_logit(logit_pa[n]);
  }
}
model{
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  lor ~ normal(0, 100);
  mu ~ normal(0, 100);
}
generated quantities{
  vector[2*NS] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc[n]);
    log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa[n]);
  }
  for(n in 1:NS) {
    rdi[n] = pa[n]-pc[n];
  }
  rd = mean(rdi);
}"

ctelog.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lrr;
  vector[NS] mu;
}
transformed parameters{
  vector[NS] log_pc;
  vector[NS] log_pa;
  vector[NS] lambda_c;
  vector[NS] lambda_a;
  vector[NS] offset_c;
  vector[NS] offset_a;
  //
  for(n in 1:NS) {
    offset_c[n] = log(nc[n]);
    offset_a[n] = log(na[n]);
    //
    log_pc[n] = mu[n];
    log_pa[n] = mu[n] + lrr;
    //
    lambda_c[n] = exp(log_pc[n] + offset_c[n]);
    lambda_a[n] = exp(log_pa[n] + offset_a[n]);
  }
}
model{
  ec ~ poisson(lambda_c);
  ea ~ poisson(lambda_a);
  lrr ~ normal(0, 100);
  mu ~ normal(0, 100);
}
generated quantities{
  vector[2*NS] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = poisson_lpmf(ec[n] | lambda_c[n]);
    log_lik[NS+n] = poisson_lpmf(ea[n] | lambda_a[n]);
    //
    rdi[n] =  exp(log_pa[n]) - exp(log_pc[n]);
  }
  rd = mean(rdi);
}"

htelogit.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lor;
  real mu[NS];
  real<lower=0> tau;  // sd of delta 
  real theta_tilde[NS];
}
transformed parameters{
  vector[NS] logit_pc;
  vector[NS] logit_pa;
  vector[NS] delta;
  vector[NS] pc;
  vector[NS] pa;
  //
  for(n in 1:NS) {
    delta[n] = lor + tau*tau*theta_tilde[n];
    logit_pc[n] = mu[n];
    logit_pa[n] = mu[n] + delta[n];
    //
    pc[n] = inv_logit(logit_pc[n]);
    pa[n] = inv_logit(logit_pa[n]);
  }
}
model{
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  lor ~ normal(0, 100);
  mu ~ normal(0, 100);
  theta_tilde ~ normal(0, 1);
  tau ~ uniform(0,2);
}
generated quantities{
  vector[NS*2] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc[n]);
    log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa[n]);
  }
  for(n in 1:NS) {
    rdi[n] = pa[n]-pc[n];
  }
  rd = mean(rdi);
}"

htelog.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lrr;
  real mu[NS];
  real<lower=0> tau;  // sd of delta 
  real theta_tilde[NS];
}
transformed parameters{
  vector[NS] log_pc;
  vector[NS] log_pa;
  vector[NS] lambda_c;
  vector[NS] lambda_a;
  vector[NS] offset_c;
  vector[NS] offset_a;
  vector[NS] delta;
  //
  for(n in 1:NS) {
    offset_c[n] = log(nc[n]);
    offset_a[n] = log(na[n]);
    //
    delta[n] = lrr + tau*tau*theta_tilde[n];
    log_pc[n] = mu[n];
    log_pa[n] = mu[n] + delta[n];
    //
    lambda_c[n] = exp(log_pc[n] + offset_c[n]);
    lambda_a[n] = exp(log_pa[n] + offset_a[n]);
  }
}
model{
  ec ~ poisson(lambda_c);
  ea ~ poisson(lambda_a);
  //
  lrr ~ normal(0, 100);
  mu ~ normal(0, 100);
  theta_tilde ~ normal(0, 1);
  tau ~ uniform(0,2);
}
generated quantities{
  vector[NS*2] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = poisson_lpmf(ec[n] | lambda_c[n]);
    log_lik[NS+n] = poisson_lpmf(ea[n] | lambda_a[n]);
    //
    rdi[n] =  exp(log_pa[n]) - exp(log_pc[n]);
  }
  rd = mean(rdi);
}"

hteshlogit.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lor;
  real<lower=0> tau;  // sd of delta 
  real theta_tilde[NS];
  real mu_tilde[NS];
  real mP;  // mean of mu
  real<lower=0> tauP;  // sd of mu
}
transformed parameters{
  vector[NS] logit_pc;
  vector[NS] logit_pa;
  vector[NS] delta;
  vector[NS] mu;
  vector[NS] pc;
  vector[NS] pa;
  //
  for(n in 1:NS) {
    delta[n] = lor + tau*tau*theta_tilde[n];
    mu[n] = mP + tauP*tauP*mu_tilde[n];
    //
    logit_pc[n] = mu[n];
    logit_pa[n] = mu[n] + delta[n];
    //
    pc[n] = inv_logit(logit_pc[n]);
    pa[n] = inv_logit(logit_pa[n]);
  }
}
model{
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  mu_tilde ~ normal(0, 1);
  theta_tilde ~ normal(0, 1);
  //
  lor ~ normal(0, 100);
  mP ~ normal(0, 100);
  tau ~ uniform(0,2);
  tauP ~ uniform(0,2);
}
generated quantities{
  vector[NS*2] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc[n]);
    log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa[n]);
  }
  for(n in 1:NS) {
    rdi[n] = pa[n]-pc[n];
  }
  rd = mean(rdi);
}"

hteshlog.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real lrr;
  real<lower=0> tau;  // sd of delta 
  real theta_tilde[NS];
  real mu_tilde[NS];
  real mP;  // mean of mu
  real<lower=0> tauP;  // sd of mu
}
transformed parameters{
  vector[NS] log_pc;
  vector[NS] log_pa;
  vector[NS] lambda_c;
  vector[NS] lambda_a;
  vector[NS] offset_c;
  vector[NS] offset_a;
  vector[NS] delta;
  vector[NS] mu;
  //
  for(n in 1:NS) {
    offset_c[n] = log(nc[n]);
    offset_a[n] = log(na[n]);
    //
    delta[n] = lrr + tau*tau*theta_tilde[n];
    mu[n] = mP + tauP*tauP*mu_tilde[n];
    //
    log_pc[n] = mu[n];
    log_pa[n] = mu[n] + delta[n];
    //
    lambda_c[n] = exp(log_pc[n] + offset_c[n]);
    lambda_a[n] = exp(log_pa[n] + offset_a[n]);
  }
}
model{
  ec ~ poisson(lambda_c);
  ea ~ poisson(lambda_a);
  //
  mu_tilde ~ normal(0, 1);
  theta_tilde ~ normal(0, 1);
  //
  lrr ~ normal(0, 100);
  mP ~ normal(0, 100);
  tau ~ uniform(0,2);
  tauP ~ uniform(0,2);
}
generated quantities{
  vector[NS*2] log_lik;
  vector[NS] rdi;  // study-specific risk difference
  real rd;
  for(n in 1:NS) {
    log_lik[n] = poisson_lpmf(ec[n] | lambda_c[n]);
    log_lik[NS+n] = poisson_lpmf(ea[n] | lambda_a[n]);
    //
    rdi[n] =  exp(log_pa[n]) - exp(log_pc[n]);
  }
  rd = mean(rdi);
}"

#### AB logit using Cholesky decomposition
ablogit.stan = "
data{
  int NS;  // number of studies 
  int NT;  // number of treatments
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
  matrix[NT, NT] L_inv;  // Cholesky factor of scale matrix (inversed); Sigma ~ inv_Wishart(nu, L_inv^T*L_inv)
  int nu;
  vector[NT] mean0;  // mean (0,0) for eta 
}
parameters{
  real logit_risk_c;  // logit(risk) in control
  real logit_risk_a;  // logit(risk) in active
  row_vector[NT] eta[NS];  // eta is a NS by NT matrix 
  vector<lower=0>[NT] c;   // For Cholesky factor
  real z;  // This should be $vector[0.5*NT*(NT-1)] z$ but it is not vector when NT=2
}
transformed parameters{
  vector[NS] logit_pc;
  vector[NS] logit_pa;
  vector[NS] pc;
  vector[NS] pa;
  matrix[NT, NT] A;
  matrix[NT, NT] A_inv_L_inv;
  matrix[NT, NT] R;
  //
  for(n in 1:NS) {
    logit_pc[n] = logit_risk_c + eta[n,1];
    logit_pa[n] = logit_risk_a + eta[n,2];
    //
    pc[n] = inv_logit(logit_pc[n]);
    pa[n] = inv_logit(logit_pa[n]);
  }
  // Define matrix A element by element because it is 2 by 2
  A[1,1]=sqrt(c[1]);
  A[1,2]=0;
  A[2,1]=z;
  A[2,2]=sqrt(c[2]);
  A_inv_L_inv = mdivide_left_tri_low(A, L_inv);
  R = crossprod(A_inv_L_inv);
}
model{
  for(i in 1:NT) {
    c[i] ~ chi_square(nu-i+1);
  }
  z ~ normal(0,1);  // Sigma=crossprod(A_inv_L_inv) ~ inv_wishart(nu, L_inv` * L_inv)
  eta ~ multi_normal(mean0, R); 
  //
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  logit_risk_c ~ normal(0, 100);
  logit_risk_a ~ normal(0, 100);
}
generated quantities{
  vector[NS*2] log_lik;
  real rd;
  real lor;
  real phat_c;
  real phat_a;
  real tau_c;
  real tau_a;
  //
  for(n in 1:NS) {
      log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc[n]);
      log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa[n]);
  }
  lor = logit_risk_a - logit_risk_c;
  phat_c = inv_logit(logit_risk_c);
  phat_a = inv_logit(logit_risk_a);
  rd = phat_a - phat_c;
  tau_c = sqrt(R[1,1]);
  tau_a = sqrt(R[2,2]);
}"

ablog.stan = "
data{
  int NS;  // number of studies 
  int NT;  // number of treatments
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
  matrix[NT, NT] L_inv;  // Cholesky factor of scale matrix (inversed); Sigma ~ inv_Wishart(nu, L_inv^T*L_inv)
  int nu;
  vector[NT] mean0;  // mean (0,0) for eta 
}
parameters{
  real log_risk_c;  // log(risk) in control
  real log_risk_a;  // log(risk) in active
  row_vector[NT] eta[NS];  // eta is a NS by NT matrix 
  vector<lower=0>[NT] c;   // For Cholesky factor
  real z;  // This should be $vector[0.5*NT*(NT-1)] z$ but it is not vector when NT=2
}
transformed parameters{
  vector[NS] log_pc;
  vector[NS] log_pa;
  vector[NS] lambda_c;
  vector[NS] lambda_a;
  vector[NS] offset_c;
  vector[NS] offset_a;
  matrix[NT, NT] A;
  matrix[NT, NT] A_inv_L_inv;
  matrix[NT, NT] R;
  //
  for(n in 1:NS) {
    offset_c[n] = log(nc[n]);
    offset_a[n] = log(na[n]);
    //
    log_pc[n] = log_risk_c + eta[n,1];
    log_pa[n] = log_risk_a + eta[n,2];
    //
    lambda_c[n] = exp(log_pc[n] + offset_c[n]);
    lambda_a[n] = exp(log_pa[n] + offset_a[n]);
  }
  // Define matrix A element by element because it is 2 by 2
  A[1,1]=sqrt(c[1]);
  A[1,2]=0;
  A[2,1]=z;
  A[2,2]=sqrt(c[2]);
  A_inv_L_inv = mdivide_left_tri_low(A, L_inv);
  R = crossprod(A_inv_L_inv);
}
model{
  for(i in 1:NT) {
    c[i] ~ chi_square(nu-i+1);
  }
  z ~ normal(0,1);  // Sigma=crossprod(A_inv_L_inv) ~ inv_wishart(nu, L_inv` * L_inv)
  eta ~ multi_normal(mean0, R); 
  //
  ec ~ poisson(lambda_c);
  ea ~ poisson(lambda_a);
  //
  log_risk_c ~ normal(0, 100);
  log_risk_a ~ normal(0, 100);
}
generated quantities{
  vector[NS*2] log_lik;
  real rd;
  real lrr;
  real phat_c;
  real phat_a;
  real tau_c;
  real tau_a;
  //
  for(n in 1:NS) {
      log_lik[n] = poisson_lpmf(ec[n] | lambda_c[n]);
      log_lik[NS+n] = poisson_lpmf(ea[n] | lambda_a[n]);
  }
  phat_c = exp(log_risk_c);
  phat_a = exp(log_risk_a);
  lrr = log_risk_a - log_risk_c;
  rd = phat_a - phat_c;
  tau_c = sqrt(R[1,1]);
  tau_a = sqrt(R[2,2]);
}"

ctebeta.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real pc;
  real pa;
}
model{
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  pc ~ beta(1,1);
  pa ~ beta(1,1);
}
generated quantities{
  vector[2*NS] log_lik;
  real lor;
  real lrr;
  real rd;
  for(n in 1:NS) {
    log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc);
    log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa);
  }
  lor=logit(pa)-logit(pc);
  lrr=log(pa/pc);
  rd=pa-pc;
}"

htebeta.stan = "
data{
  int NS;  // number of studies 
  int s[NS];  // study id
  int ec[NS];  // number of events in control
  int ea[NS];  // number of events in active
  int nc[NS];  // sample size in control
  int na[NS];  // sample size in active
}
parameters{
  real<lower=0,upper=1> Uc;  // mean of risk in control
  real<lower=0,upper=1> Ua;  // mean of risk in active
  real<lower=0> Vc;
  real<lower=0> Va;
  vector<lower=0,upper=1>[NS] pc;
  vector<lower=0,upper=1>[NS] pa;
}
transformed parameters{
  real<lower=0> alpha_c;
  real<lower=0> alpha_a;
  real<lower=0> beta_c;
  real<lower=0> beta_a;
  //
  alpha_c = Uc*Vc;
  beta_c = (1-Uc)*Vc;
  alpha_a = Ua*Va;
  beta_a = (1-Ua)*Va;
}
model{
  Uc ~ beta(1,1);
  Ua ~ beta(1,1);
  Vc ~ gamma(1, 0.01);  // mean=100 (alpha/beta)
  Va ~ gamma(1, 0.01);  // mean=100 (alpha/beta)
  ec ~ binomial(nc, pc);
  ea ~ binomial(na, pa);
  //
  pc ~ beta(alpha_c, beta_c);
  pa ~ beta(alpha_a, beta_a);
}
generated quantities{
  vector[2*NS] log_lik;
  real lor;
  real lrr;
  real rd;
  for(n in 1:NS) {
    log_lik[n] = binomial_lpmf(ec[n] | nc[n], pc[n]);
    log_lik[NS+n] = binomial_lpmf(ea[n] | na[n], pa[n]);
  }
  lor=logit(Ua)-logit(Uc);
  lrr=log(Ua/Uc);
  rd=Ua-Uc;
}"




