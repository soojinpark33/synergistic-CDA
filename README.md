# Causal Decomposition Analysis with Synergistic Interventions: A Triply-Robust Machine Learning Approach to Addressing Multiple Dimensions of Social Disparities
R codes for "Causal Decomposition Analysis with Synergistic Interventions: A Triply-Robust Machine Learning Approach to Addressing Multiple Dimensions of Social Disparities"

Soojin Park<sup>1</sup>, Suyeon Kim<sup>1</sup>, and Chioun Lee<sup>2</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Sociology, University of California, Riverside


## Overview

Educational disparities are rooted in, and perpetuate, social inequalities across multiple dimensions such as race, socioeconomic status, and geography. To reduce disparities, most intervention strategies focus on a single domain and frequently evaluate their effectiveness by using causal decomposition analysis. However, a growing body of research suggests that single-domain interventions may be insufficient for individuals marginalized on multiple fronts. While interventions across multiple domains are increasingly proposed, there is limited guidance on appropriate methods for evaluating their effectiveness. To address this gap, we develop an extended causal decomposition analysis that simultaneously targets multiple causally ordered intervening factors, allowing for the assessment of their synergistic effects. These scenarios often involve challenges related to model misspecification due to complex interactions among group categories, intervening factors, and their confounders with the outcome. To mitigate these challenges, we introduce a triply-robust estimator that leverages machine learning techniques to address potential model misspecification. We apply our method to a cohort of students from the High School Longitudinal Study (HSLS:2009), focusing on math achievement disparities between Black, Hispanic, and White high schoolers. Specifically, we examine how two sequential interventions—equalizing the proportion of students who attend high-\textcolor{black}{performing} schools and equalizing enrollment in Algebra I by 9th grade across racial groups—may reduce these disparities. 

For more details of our proposed methods, see [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis. 

## Case Study

* `data2.Rdata` 
  
  For our case study, we used data from the High School Longitudinal Study 2009 (HSLS:09) study. The HSLS09 data is restricted from circulation, and the original data can be downloaded from the NCES portal by clicking [here]([https://nces.ed.gov/surveys/hsls09/hsls09_data.asp]). 

* `casestudy_main_code.R` 
 
   This `R` file replicates Tables 1 and 2 of our study.

* `casestudy_source_code.R` 
 
   This `R` file contains functions for 'casestudy_main_code.R'.

## Simulation Study

* `main_code_sim.R`  

   This `R` file contains the simulation codes for our propposed CDA with synergistic interventions. This code replicates Figure 2 of our paper.

* `source_code_park_cv.R` 
 
   This `R` file includes source functions required to run our simulation codes. 

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.
