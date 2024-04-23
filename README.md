# Evolving trends in racial and ethnic disparities in syphilis prevalence among pregnant women in the United States from 2014 to 2022
The data and code for paper title "Evolving trends in racial and ethnic disparities in syphilis prevalence among pregnant women in the United States from 2014 to 2022"

The outline:

- data: data input for parameterization
  - natality: natality data
    - raw: raw datasets
    - selected: by year dataset with only selected variables
  - fetal_death: fetal deaths by race and year
  - parameter: prior info for each parameter in the model
  - test_sens_spec: Sensitivity and Specificity for the syphilis test using for meta-analysis

- code:
  - 0_functions: pre-defined functions
  - 1_parameterization: the determination of each parameter in the model
  - 2_analysis:
    - 1_main: the main analysis
    - 2_SA1: the sensitivity analysis with an alternative prior for proportion of women screened by race/ethnicity
    - 3_SA2: the sensitivity analysis with alternative weekly informative priors for syphilis prevalence
  - 3_plot: the creation of tables and figures in the manuscript
- output: 
  - stan_output: simulated draws from Stan programs
  - figure:
    - params: plots for priors for some parameters
    - est_prevalence: plots for posteriors for syphilis prevalence by year and race
  - table: posterior summary of the model
