rm(list = ls())

pacman::p_load(rriskDistributions, tidyverse, rio, here)

source(here("code", "0_functions", "get_mean_beta.R"))
source(here("code", "0_functions", "get_range_beta.R"))
source(here("code", "0_functions", "get_beta_normal.R"))

test_prop =
  list(
    c(62.2, 73.6, 86.3, 60.9, 59.2, 74.4) / 100,
    c(74.1, 75.7, 91.4, 64.2, 63.2, 78.8) / 100,
    c(20.7, 58.7, 91.2, 33.6, 25.4, 60.0) / 100,
    c(53.6, 66.5, 91.2, 55, 38.7, 71.3) / 100,
    c(62.1, 73.5, 89.1, 61.3, 60.8, 80) / 100,
    c(74, 77, 92.6, 64.9, 63.4, 82.7) / 100,
    c(35.8, 57.3, 91.8, 36.9, 26.6, 80) / 100,
    c(55.5, 64.9, 91.3, 56.7, 42.1, 69.4) / 100
  )

test_num = list(
  c(17020, 15367, 8484, 15381, 7165, 14786),
  c(20344, 2795, 11104, 13042, 8298, 7939),
  c(330, 910, 1185, 3072, 449, 12),
  c(2399, 2056, 10180, 2968, 1785, 6970),
  c(16094, 14596, 7479, 15084, 7354, 15892),
  c(19439, 2730, 10090, 12953, 8321, 8335),
  c(625, 884, 927, 3357, 470, 16),
  c(2901, 2174, 11248, 2937, 1939, 6614)
)


# 1 main analysis -------------------------------------------------------------------------------------------------

# white, black, hispanic

test_prop1 = tibble(
  year = rep(2017:2018, each = 4),
  race = rep(c(
    "Non-H White", "Non-H Black", "Hispanic", "OthersMissings"
  ), 2),
  test_prop = test_prop,
  test_num = test_num
) |> 
  mutate(
    denominator = map2_dbl(test_prop, test_num,
                           ~ sum(.y / .x)),
    numerator = map_dbl(test_num, sum)
  ) |> 
  select(-test_prop,-test_num) |>
  group_by(race) |>
  summarise(denominator = sum(denominator),
            numerator = sum(numerator)) |>
  mutate(pool_test_prop = numerator / denominator) |>
  mutate(
    pool_test_prop_sd = map2_dbl(pool_test_prop, denominator,
                                 ~ get_range_beta(.x, .y)$sd),
    pool_test_prop_lower = map2_dbl(pool_test_prop, denominator,
                                    ~ get_range_beta(.x, .y)$lower),
    pool_test_prop_upper = map2_dbl(pool_test_prop, denominator,
                                    ~ get_range_beta(.x, .y)$upper),
    params = map2(
      pool_test_prop,
      pool_test_prop_sd,
      ~ get_beta_normal(mu = .x, sd = .y)
    )
  ) |>
  unnest_wider(params) |>
  rename(Test_alpha = alpha, Test_beta = beta) |>
  mutate(
    simu = map2(Test_alpha, Test_beta, ~ rbeta(
      n = 1000, shape1 = .x, shape2 = .y
    )),
    mean_simu = map_dbl(simu, mean),
    sd_simu = map_dbl(simu, sd),
    lower_simu = map_dbl(simu, ~ quantile(., probs = 0.025)),
    upper_simu = map_dbl(simu, ~ quantile(., probs = 0.975))
  ) |> 
  select(race, Test_alpha, Test_beta)

# AIAN, Asian, Multiracial, NHPI

test_prop2 = tibble(
  race = c("Non-H AIAN", "Non-H Asian", "Non-H Multiracial", "Non-H NHPI"),
  test_prop_lower = 0.4004618,
  test_prop_upper = 0.7346180
) |>
  mutate(
    params = map2(
      test_prop_lower,
      test_prop_upper,
      ~ get.beta.par(
        p = c(0.025, 0.975),
        q = c(.x, .y),
        show.output = FALSE,
        plot = FALSE
      )
    )
  ) |>
  unnest_wider(params) |>
  rename(Test_alpha = shape1, Test_beta = shape2) |>
  mutate(
    simu = map2(Test_alpha, Test_beta, ~ rbeta(
      n = 10000, shape1 = .x, shape2 = .y
    )),
    mean = map_dbl(simu, mean),
    sd = map_dbl(simu, sd),
    lower = map_dbl(simu, ~ quantile(., probs = 0.025)),
    upper = map_dbl(simu, ~ quantile(., probs = 0.975))
  ) |>
  select(race, Test_alpha, Test_beta)

# save
test_prop1 |> 
  add_row(test_prop2) |> 
  filter(race != "OthersMissings") |> 
  mutate(race = str_remove(race, "Non-H ")) |> 
  arrange(race) |> 
  export(here(
      "data",
      "parameter",
      "3_prop_screened.csv"
  ))

# 2 sensitivity analysis ------------------------------------------------------------------------------------------

rm(list = ls())
pacman::p_load(rriskDistributions, tidyverse, rio, here)

births = import(here("data", "natality", "birth_all_coded.rds")) |> as_tibble()

test_prop_yr_race = tibble(
  year = rep(2017:2018, each = 4),
  race = rep(c(
    "White", "Black", "Hispanic", "OthersMissings"
  ), 2),
  test_prop =
    list(
      c(62.2, 73.6, 86.3, 60.9, 59.2, 74.4) / 100,
      c(74.1, 75.7, 91.4, 64.2, 63.2, 78.8) / 100,
      c(20.7, 58.7, 91.2, 33.6, 25.4, 60.0) / 100,
      c(53.6, 66.5, 91.2, 55, 38.7, 71.3) / 100,
      c(62.1, 73.5, 89.1, 61.3, 60.8, 80) / 100,
      c(74, 77, 92.6, 64.9, 63.4, 82.7) / 100,
      c(35.8, 57.3, 91.8, 36.9, 26.6, 80) / 100,
      c(55.5, 64.9, 91.3, 56.7, 42.1, 69.4) / 100
    ),
  test_num =
    list(
      c(17020, 15367, 8484, 15381, 7165, 14786),
      c(20344, 2795, 11104, 13042, 8298, 7939),
      c(330, 910, 1185, 3072, 449, 12),
      c(2399, 2056, 10180, 2968, 1785, 6970),
      c(16094, 14596, 7479, 15084, 7354, 15892),
      c(19439, 2730, 10090, 12953, 8321, 8335),
      c(625, 884, 927, 3357, 470, 16),
      c(2901, 2174, 11248, 2937, 1939, 6614)
    )
) |>
  mutate(
    numerator = map_dbl(test_num, ~ sum(.x)),
    denominator = map2_dbl(test_prop, test_num,
                           ~ sum(.y / .x)),
    pool_test_prop = map2_dbl(test_prop, test_num,
                              ~ sum(.y) / sum(.y / .x))
  )

test_prop_race <- test_prop_yr_race |>
  select(year, race, numerator, denominator) |>
  group_by(race) |>
  summarise(numerator = sum(numerator),
            denominator = sum(denominator)) |>
  mutate(pool_test_prop = numerator / denominator) |>
  add_row(race = c("AIAN", "Asian", "NHPI", "Multiracial")) |>
  mutate(pool_test_prop = ifelse(!is.na(pool_test_prop), 
                                 pool_test_prop, 
                                 pool_test_prop[race == "OthersMissings"])) |> 
  filter(race != "OthersMissings") |> 
  select(race, lower = pool_test_prop)

pnc_prop_race = births |>
  select(year, race, No_PNC) |>
  drop_na() |>
  mutate(pnc = 1 - No_PNC) |>
  group_by(race) |>
  count(pnc) |>
  mutate(pnc_prop = n / sum(n)) |>
  filter(pnc == 1) |>
  select(race, upper = pnc_prop)

test_prop_race |> 
  left_join(pnc_prop_race, by = "race") |> 
  mutate(
    params = map2(lower, upper, ~ get.beta.par(
      p = c(0.025, 0.975),
      q = c(.x, .y),
      show.output = FALSE,
      plot = FALSE
    ))
  ) |> 
  unnest_wider(params) |> 
  select(race, Test_alpha = shape1, Test_beta = shape2) |> 
  arrange(race) |> 
  export(here(
    "data",
    "parameter",
    "3_prop_screened_SA1.csv"
  ))

