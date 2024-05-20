

rm(list = ls())

pacman::p_load(tidyverse,
               tidybayes,
               here,
               rio,
               rstan)

# model in parallel
options(mc.cores = parallel::detectCores())
# says to excute
rstan_options(auto_write = TRUE)
# helper functions
source(here("code", "0_functions", "order_levels.R"))


# 1 parameter load ------------------------------------------------------------------------------------------------

## positivities
positivities_yr_race = import(here("data",
                                   "parameter",
                                   "0_positivities.csv")) |>
  order_levels(c("year", "race"))

## sensitivity
sens = import(here("data",
                   "parameter",
                   "1_sens.csv"))

## specificity
spec = import(here("data",
                   "parameter",
                   "2_spec.csv"))

## prop of screened
screened_race = import(here("data",
                            "parameter",
                            "3_prop_screened.csv")) |>
  order_levels("race")

## pregnancies
pregnancies_yr_race = import(here("data",
                                  "parameter",
                                  "4_pregnancies.csv")) |>
  order_levels(c("year", "race")) |>
  select(year, race, Pop = births)

# 2 Data preparation for Stan -------------------------------------------------------------------------------------

df = positivities_yr_race |>
  left_join(pregnancies_yr_race, by = c("year", "race")) |>
  left_join(screened_race, by = "race") |>
  mutate(
    Sens_alpha = sens$Sens_alpha,
    Sens_beta = sens$Sens_beta,
    Spec_alpha = spec$Spec_alpha,
    Spec_beta = spec$Spec_beta,
    P_alpha = 1,
    P_beta = 5
  )

data_ls = df |>
  compose_data(
    Sens_alpha = Sens_alpha[1],
    Sens_beta = Sens_beta[1],
    Spec_alpha = Spec_alpha[1],
    Spec_beta = Spec_beta[1],
    P_alpha = P_alpha[1],
    P_beta = P_beta[1]
  )

# 3 Fit Bayesian model --------------------------------------------------------------------------------------------

fit = stan(
  file = here("code", "2_analysis", "year_race_prevalence.stan"),
  data = data_ls,
  iter = 12000,
  chains = 4,
  warmup = 8000,
  control = list(adapt_delta = 0.995, max_treedepth = 12),
  seed = 123 #123 #547
)

fit_tidy = fit |>
  recover_types(df) |>
  spread_draws(P[year, race],
               Sens, Spec,
               Test[year, race],
               theta[year, race],
               pred_Positivities[year, race]) |>
  ungroup()

export(fit_tidy,
       here("output", "stan_output", "1_main.rds"))
