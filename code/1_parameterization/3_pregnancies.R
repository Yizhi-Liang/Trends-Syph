
# Notes: 3 parts in this script
# Part 1: process the raw natality data and combine them together
# Part 2: recode values to make variables readable
# Part 3: get the number of pregnancies by year and race/ethnicity

# 1 pre-prosessing natality data ----------------------------------------------------------------------------------

rm(list = ls())
pacman::p_load(data.table, here, rio)

## select variables
selected_vars <- c(
  "dob_yy", "mager", "mracehisp", "bfacil",
  "meduc", "precare", "previs", "pay_rec", "dbwt", "oegest_comb",
  "ip_gon", "ip_syph", "ip_chlam", "dplural",
  "ab_anti", "ld_antb", "dmar", "mbstate_rec",
  "cig0_r", "cig1_r", "cig2_r", "cig3_r",
  "bmi", "priorlive", "priordead", "priorterm"
)

## process data
process_year <- function(year) {
  # Construct the file name and path
  file_name <- paste0("natality", year, "us", ".csv")
  file_path <- here("data", "natality", "raw", file_name)
  
  # Read the data
  data <- fread(file_path, select = selected_vars)
  
  export_file_name <- paste0("birth_", year, ".rds")
  export_file_path <- here("data", "natality", "selected", export_file_name)
  
  export(data, export_file_path)
  
  return(paste0("Processed year ", year))
}

# export each year
lapply(2014:2022, process_year)

# read each year's data
read_parquet_year <- function(year) {

  file_name <- paste0("birth_", year, ".rds")
  file_path <- here("data", "natality", "selected", file_name)
  
  # Read the data
  data <- import(file_path)
  
  return(data)
}

all_dat <- rbindlist(lapply(2014:2022, read_parquet_year))

export(all_dat, here("data", "natality", "birth_all.rds"))

# 2 Recode variables ----------------------------------------------------------------------------------------------

rm(list = ls())
pacman::p_load(data.table, rio, here)

full_data = as.data.table(import(here("data", 
                                      "natality", 
                                      "birth_all.rds")))

## year
clean_dat <- copy(full_data)
clean_dat[, year := as.numeric(dob_yy)]

## age
clean_dat[, age := as.numeric(mager)] 

## race/ethnicity
clean_dat[, race := factor(
  dplyr::case_when(
    mracehisp == 1 ~ "White",
    mracehisp == 2 ~ "Black",
    mracehisp == 3 ~ "AIAN",
    mracehisp == 4 ~ "Asian",
    mracehisp == 5 ~ "NHPI",
    mracehisp == 6 ~ "Multiracial",
    mracehisp == 7 ~ "Hispanic",
    mracehisp == 8 ~ NA_character_,
    .default = NA_character_
  )
)]

## education
clean_dat[, educ := factor(
  dplyr::case_when(
    meduc %in% 1:2 ~ "<HS",
    meduc %in% 3:4 ~ "HS",
    meduc %in% 5:6 ~ "College",
    meduc %in% 7:8 ~ ">College",
    .default = NA_character_
  ),
  levels = c("<HS", "HS", "College", ">College")
)]

## insurance
clean_dat[, insurance := factor(
  dplyr::case_when(
    pay_rec == 1 ~ "medicaid",
    pay_rec == 2 ~ "private",
    pay_rec == 3 ~ "self",
    pay_rec == 4 ~ "other",
    .default = NA_character_
  ),
  levels = c("medicaid", "private", "self", "other")
)]

## birthplace
clean_dat[, birthplace := case_when(
  bfacil == 1 ~ "hospital",
  bfacil == 2 ~ "freestanding",
  bfacil == 3 ~ "home_intended",
  bfacil == 4 ~ "home_notintended",
  bfacil == 5 ~ "home_unknownintended",
  bfacil == 6 ~ "clinic_doctor_office",
  bfacil == 7 ~ "other",
  .default = NA_character_
)]

## plural
clean_dat[, birth_num :=
            dplyr::case_when(
              dplural %in% 1:5 ~ as.numeric(dplural),
              .default = NA_integer_
            )]

## infectious
clean_dat[, c("gon", "syph", "chlam") := list(
  ifelse(ip_gon == "Y", 1, ifelse(ip_gon == "N", 0, NA_integer_)),
  ifelse(ip_syph == "Y", 1, ifelse(ip_syph == "N", 0, NA_integer_)),
  ifelse(ip_chlam == "Y", 1, ifelse(ip_chlam == "N", 0, NA_integer_))
)]

## marital
clean_dat[, marital := dplyr::case_when(
  dmar == 9 ~ NA_integer_,
  dmar == 1 ~ 1,
  dmar == 2 ~ 0,
  .default = NA_integer_
)]

## BMI
clean_dat[, BMI := dplyr::case_when(
  bmi %in% 13:70 ~ as.numeric(bmi),
  .default = NA
)]

## smoke
clean_dat[, c("smoke_bf", "smoke_preg") := list(
  dplyr::case_when(
    cig0_r == 6 ~ NA_integer_,
    cig0_r == 0 ~ 0,
    .default = 1
  ),
  dplyr::case_when(
    cig1_r == 6 | cig2_r == 6 | cig3_r == 6 ~ NA_integer_,
    cig1_r == 0 | cig2_r == 0 | cig3_r == 0 ~ 0,
    .default = 1
  )
)]

## PNC
clean_dat[, c("initial_month", "visit_num", "gest_wks", "birth_wt") := list(
  ifelse(precare %in% 0:10, precare, NA_integer_), # initial month
  ifelse(previs %in% 0:98, previs, NA_integer_), # visit number
  ifelse(oegest_comb %in% 17:47, oegest_comb, NA_integer_), # gestational age
  ifelse(dbwt %in% 227:8165, dbwt, NA_integer_) # birth weight (grams)
)]

clean_dat[, No_PNC := dplyr::case_when(
  initial_month == 0 | visit_num == 0 ~ 1,
  is.na(initial_month) | is.na(visit_num) ~ NA_integer_,
  .default = 0
)]

## save the data
cols_to_keep <- names(clean_dat)[which(names(clean_dat) == "year"):which(names(clean_dat) == "No_PNC")]
final <- clean_dat[, ..cols_to_keep]
export(final, here("data", "natality", "birth_all_coded.rds"))

# 3 Num. of pregnancies -------------------------------------------------------------------------------------------

rm(list = ls())
pacman::p_load(tidyverse, rio, here)

pregnancies_yr_race = import(here("data", "natality", "birth_all_coded.rds")) |> 
  drop_na(race) |> 
  group_by(year, race) |> 
  summarise(births = n()) |> 
  ungroup() |> 
  left_join(import(here("data", "fetal_death", "fetal_death_yr_race.csv")) |> 
              mutate(race = str_remove(race, "Non-H ")) |> 
              rename(fetal_deaths = fetal_death),
            by = c("year", "race")) |> 
  mutate(
    death_rate = 1000 * fetal_deaths / births,
    death_rate = ifelse(year == 2022, death_rate[year==2021], death_rate),
    fetal_deaths = ifelse(year == 2022, births*death_rate/1000, fetal_deaths),
    total = births + fetal_deaths
  ) |> 
  select(year, race, births, fetal_deaths, total)

pregnancies_yr_race |> 
  export(here(
    "data",
    "parameter",
    "4_pregnancies.csv"
  ))

import(here("data", "natality", "birth_all_coded.rds")) |> 
  filter(insurance == "medicaid") |> 
  group_by(year) |> 
  summarise(births = n())
