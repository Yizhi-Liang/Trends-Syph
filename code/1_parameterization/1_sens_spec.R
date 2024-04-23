rm(list = ls())

pacman::p_load(rriskDistributions, tidyverse, rio, here, meta)

# 0 functions defined -----------------------------------------------------------

# mean to beta params
source(here("code", "0_functions", "get_beta_normal.R"))

# 1 load data -----------------------------------------------------------------------------------------------------

sens = import(here(
  "data",
  "test_sens_spec",
  "1_sens_park.csv"
)) |> as_tibble() |>
  mutate(Stage = factor(
    Stage,
    levels = c(
      "primary",
      "secondary",
      "early latent",
      "late latent",
      "tertiary",
      "uncertain",
      "all"
    )
  ),
  # stage: primary and secondary, early latent
  stage_1 = factor(
    case_when(
      Stage %in% c("primary", "secondary") ~ "primary_secondary",
      Stage %in% c("early latent") ~ "early_latent",
      .default = NA_character_
    )
  ))

spec = import(here(
  "data",
  "test_sens_spec",
  "2_spec_park.csv"
)) |> as_tibble()

# 2 sensitivity ---------------------------------------------------------------------------------------------------

sens_prop <- metaprop(
  event = `True Pos`,
  n = `Actual Pos`,
  studlab = Author,
  data = sens[sens$Stage %in% c("primary", "secondary", "early latent"), ],
  method = "GLMM",
  sm = "PLOGIT",
  fixed = FALSE,
  random = TRUE,
  title = "Sensitivity for Syphilis Test",
  prediction = TRUE,
  prediction.subgroup = TRUE,
  subgroup = Stage,
  tau.common = FALSE,
  overall = TRUE
)

summary(sens_prop)

Sens_all = tibble(sens = 0.9769,
       sens_lower = 0.9635,
       sens_upper = 0.9854) |>
  mutate(params = pmap(list(sens_lower, sens, sens_upper),
                       function(sens_lower, sens, sens_upper) {
                         get.beta.par(
                           q = c(sens_lower, sens, sens_upper),
                           plot = FALSE,
                           show.output = FALSE
                         )
                       })) |>
  unnest_wider(params) |> 
  rename(Sens_alpha = shape1, Sens_beta = shape2) |> 
  mutate(
    simu = map2(Sens_alpha, Sens_beta,
                ~ rbeta(10000, .x, .y)),
    mean = map_dbl(simu, mean),
    sd = map_dbl(simu, sd),
    lower = map_dbl(simu, quantile, probs = 0.025),
    upper = map_dbl(simu, quantile, probs = 0.975)
  ) |> 
  select(-simu)

export(Sens_all |> select(Sens_alpha, Sens_beta), here(
  "data",
  "parameter",
  "1_sens.csv"
))

forest(
  sens_prop,
  # layout = "subgroup",
  calcwidth.hetstat = TRUE,
  prediction = FALSE,
  prediction.subgroup = FALSE,
  colgap.forest.left = unit(5, "line"),
  width = 12,
  file = here(
    "output",
    "figure",
    "params",
    "1_sens.pdf"
  )
)

# sens_prop_2 <- metaprop(
#   event = `True Pos`,
#   n = `Actual Pos`,
#   studlab = Author,
#   data = sens[sens$stage_1 %in% c("primary_secondary", "early_latent"), ],
#   method = "GLMM",
#   sm = "PLOGIT",
#   fixed = FALSE,
#   random = TRUE,
#   title = "Sensitivity for Syphilis Test",
#   prediction = TRUE,
#   prediction.subgroup = TRUE,
#   subgroup = stage_1,
#   tau.common = FALSE,
#   overall = TRUE
# )
# 
# forest(
#   sens_prop_2,
#   # layout = "subgroup",
#   calcwidth.hetstat = TRUE,
#   prediction = FALSE,
#   prediction.subgroup = FALSE,
#   colgap.forest.left = unit(5, "line"),
#   width = 12,
#   file = here(
#     "04_plot",
#     "0_params_sets",
#     "1_sensitivity_specificity",
#     "1_sens_2.pdf"
#   )
# )

# 3 specificity ---------------------------------------------------------------------------------------------------

spec_prop <- metaprop(
  event = `True Neg`,
  n = `Actual Neg`,
  studlab = Author,
  data = spec,
  method = "GLMM",
  sm = "PLOGIT",
  fixed = FALSE,
  random = TRUE,
  title = "Sensitivity for Syphilis Test",
  prediction = TRUE,
  prediction.subgroup = TRUE,
  subgroup = Type,
  tau.common = FALSE,
  overall = TRUE
)

summary(spec_prop)

Spec_all = tibble(spec = 0.9903,
                  spec_lower = 0.9833,
                  spec_upper = 0.9944) |>
  mutate(params = pmap(list(spec_lower, spec, spec_upper),
                       function(spec_lower, spec, spec_upper) {
                         get.beta.par(
                           q = c(spec_lower, spec, spec_upper),
                           plot = FALSE,
                           show.output = FALSE
                         )
                       })) |>
  unnest_wider(params) |> 
  rename(Spec_alpha = shape1, Spec_beta = shape2) |> 
  mutate(
    simu = map2(Spec_alpha, Spec_beta,
                ~ rbeta(10000, .x, .y)),
    mean = map_dbl(simu, mean),
    sd = map_dbl(simu, sd),
    lower = map_dbl(simu, quantile, probs = 0.025),
    upper = map_dbl(simu, quantile, probs = 0.975)
  ) |> 
  select(-simu)

export(Spec_all |> select(Spec_alpha, Spec_beta), here(
  "data",
  "parameter",
  "2_spec.csv"
))


forest(
  spec_prop,
  layout = "jama",
  calcwidth.hetstat = TRUE,
  prediction = FALSE,
  prediction.subgroup = FALSE,
  colgap.forest.left = unit(2, "line"),
  width = 12,
  file = here(
    "output",
    "figure",
    "params",
    "1_spec.pdf"
  )
)

