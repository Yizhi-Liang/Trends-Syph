# Trends in syphilis prevalence by race and ethnicity among people who are pregnant  in the United States 2016–2023.

The data and code for paper titled "Trends in syphilis prevalence by race and ethnicity among people who are pregnant  in the United States 2016–2023
."

The outline:

- data: data input for parameterization
  - natality: natality data
  - fetal_death: fetal deaths by race and year
  - parameter: prior data sources for each parameter in the model
  - test_sens_spec: Sensitivity and Specificity for the syphilis test used for meta-analysis
  
- code:
  - 0_functions: pre-defined functions
  - 1_parameterization: the determination of each parameter in the model
  - 2_analysis:
    - 1_scenario1: estimating the posterior of syphilis prevalence by year and race with a lower screening coverage from Medicaid for each race and ethnicity.
    - 2_scenario2: estimating the posterior of syphilis prevalence by year and race with a higher screening coverage, where the lower bound is the screening coverage from Medicaid and the upper bound is the proportion of women receiving at lease one prenatal care, for each race and ethnicity.
    - year_race_prevalence.stan: The Stan code for MCMC sampling and Bayesian inference.
  - 3_plot: the generation of tables and figures in the manuscript
- output: 
  - stan_output: simulated draws from Stan programs
  - figure:
    - params: plots for priors for some parameters
    - est_prevalence: plots for posteriors for syphilis prevalence by year and race
  - table: posterior summary of the model
