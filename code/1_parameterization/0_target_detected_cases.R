rm(list = ls())
pacman::p_load(tidyverse, rio, here)

syph_yr_race = import(here("data", "natality", "birth_all_coded.rds"), 
                      trust = TRUE) |> 
  drop_na(race) |> 
  group_by(year, race) |> 
  count(syph, name = "Detected") |> 
  filter(syph == 1) |> 
  select(-syph)

syph_yr_race |> 
  export(here(
    "data",
    "parameter",
    "0_detected_cases.csv"
  ))

# import(here("data", "natality", "birth_all_coded.rds")) |> 
#   filter(insurance == "medicaid") |> 
#   drop_na(race) |> 
#   group_by(year, race) |> 
#   count(syph, name = "Detected") |> 
#   filter(syph == 1) |> 
#   select(-syph) |> 
#   export(here(
#     "data",
#     "parameter",
#     "0_positivities_medicaid.csv"
#   ))
