# MetaAnalysisRareEvents
Hwanhee Hong (hwanhee.hong@duke.edu)

R code to run meta-analysis with binary (rare) events

Analysis of rosiglitazone data

* Related paper

Hwanhee Hong, Chenguang Wang, and Gary L. Rosner (2019+) Meta-Analysis of Rare Adverse Events in Randomized Clinical Trials: Bayesian and Frequentist Methods  

* Data source

Nissen SE and Wolski K (2010) Rosiglitazone revisited: an updated meta-analysis of risk for myocardial infarction and cardiovascular mortality. *Archives of internal medicine*. 170(14):1191-1201.

This code is to fit meta-analysis models (frequentist and Bayesian) to the rosiglitazone data.
This code provides 15 log odds ratio (LORs), 13 log relative risks (LRRs), and 17 risk differences estimates.
Bayesian models are fitted using rstan. Stan needs to be installed.
For Bayesian model comparison, we estimate Watanabe-Akaike or widely applicable information criterion (WAIC).

This code produces 4 files:
    1. LOR, LRR, RD estimates (.csv)
    2. WAIC.txt
    3. Bayesian model results from Stan (.txt)
    4. Bayesian model diagnostic plots (.pdf)
