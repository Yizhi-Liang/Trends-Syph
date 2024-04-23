rm(list = ls())
pacman::p_load(tidyverse, rio, here)

syph_yr_race = import(here("data", "natality", "birth_all_coded.rds")) |> 
  drop_na(race) |> 
  group_by(year, race) |> 
  count(syph, name = "Positivities") |> 
  filter(syph == 1) |> 
  select(-syph)

syph_yr_race |> 
  export(here(
    "data",
    "parameter",
    "0_positivities.csv"
  ))
  