
rm(list = ls())

pacman::p_load(tidyverse, ggplot2, here, rio, tidybayes, cowplot)

theme_set(theme_half_open())
source(here("code", "0_functions", "order_levels.R"))

size_fig = 18

# color set
detected_color = "#40a9ff"
prevalence_color_S1 = "#a4262c"
prevalence_color_SA = "#ff80ab"
prevalence_color_S2 = "#ef6c00"
diagnoses_color = "#ffb900"

# 0 load data -----------------------------------------------------------------------------------------------------

pregnancies_yr_race = import(here("data", "parameter", "4_pregnancies.csv")) |>
  order_levels(c("year", "race")) |>
  select(year, race, Pop_mean = births, fetal_deaths, total) |>
  group_by(year) |>
  mutate(wt = Pop_mean / sum(Pop_mean)) |>
  ungroup()

## estimated prevalence in scenario 1
fit_yr_race_S1 = import(here("output", "stan_output", "1_S1.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

P_pregnant_yr_S1 <- fit_yr_race_S1 |>
  group_by(year, race) |>
  summarise(
    P_mean = mean(P),
    P_lower = quantile(P, probs = 0.025),
    P_upper = quantile(P, probs = 0.975)
  ) |>
  ungroup() |>
  left_join(pregnancies_yr_race, by = c("year", "race")) |>
  mutate(
    numerator = P_mean * Pop_mean,
    numerator_lower = P_lower * Pop_mean,
    numerator_upper = P_upper * Pop_mean
  ) |>
  group_by(year) |>
  summarise(
    numerator = sum(numerator),
    numerator_lower = sum(numerator_lower),
    numerator_upper = sum(numerator_upper),
    denominator = sum(Pop_mean)
  ) |>
  mutate(
    P_mean = numerator / denominator,
    P_lower = numerator_lower / denominator,
    P_upper = numerator_upper / denominator
  ) |>
  select(year, P_mean, P_lower, P_upper)

## estimated prevalence in scenario 2
fit_yr_race_S2 = import(here("output", "stan_output", "2_S2.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

P_pregnant_yr_S2 <- fit_yr_race_S2 |>
  group_by(year, race) |>
  summarise(
    P_mean = mean(P),
    P_lower = quantile(P, probs = 0.025),
    P_upper = quantile(P, probs = 0.975)
  ) |>
  ungroup() |>
  left_join(pregnancies_yr_race, by = c("year", "race")) |>
  mutate(
    numerator = P_mean * Pop_mean,
    numerator_lower = P_lower * Pop_mean,
    numerator_upper = P_upper * Pop_mean
  ) |>
  group_by(year) |>
  summarise(
    numerator = sum(numerator),
    numerator_lower = sum(numerator_lower),
    numerator_upper = sum(numerator_upper),
    denominator = sum(Pop_mean)
  ) |>
  mutate(
    P_mean = numerator / denominator,
    P_lower = numerator_lower / denominator,
    P_upper = numerator_upper / denominator
  ) |>
  select(year, P_mean, P_lower, P_upper)

## detected syphilis rate
detected_yr_race = import(here("data", "parameter", "0_detected_cases.csv")) |>
  order_levels(c("year", "race"))

detected_rate_yr_race = pregnancies_yr_race |>
  left_join(detected_yr_race, by = c('year', 'race')) |>
  mutate(detected_rate = Detected / Pop_mean) |>
  select(year, race, detected_rate)

## CDC reported diagnosis rate among women at reproductive ages
syph_diag_yr_race <- import(here("data", "surveillance", "atlas_syphilis_yr_race.csv")) |>
  as_tibble() |>
  select(
    type = Indicator,
    year = Year,
    age_group = `Age Group`,
    race = `Race/Ethnicity`,
    cases = Cases,
    rate = `Rate per 100000`
  ) |>
  mutate(
    type = case_when(str_detect(type, "Primary and Secondary") ~ "PS", TRUE ~ "Early"),
    year = factor(as.numeric(ifelse(
      str_detect(year, "2020"), 2020, year
    ))),
    race = ifelse(race == "Unknown", NA_character_, race),
    cases = as.numeric(str_replace(cases, ",", "")),
    rate = ifelse(str_detect(rate, "Data"), NA_real_, as.numeric(rate))
  ) |>
  drop_na(race) |>
  mutate(
    race = case_when(
      race == "American Indian/Alaska Native" ~ "AIAN",
      race == "Asian" ~ "Asian",
      race == "Black/African American" ~ "Black",
      race == "Hispanic/Latino" ~ "Hispanic",
      race == "Multiracial" ~ "Multiracial",
      race == "Native Hawaiian/Other Pacific Islander" ~ "NHPI",
      race == "White" ~ "White"
    )
  ) |>
  order_levels(c("year", "race")) |>
  mutate(denominator = ifelse(rate == 0, 0, 100000 * cases / rate)) |>
  group_by(type, year, race) |>
  summarise(cases = sum(cases),
            denominator = sum(denominator)) |>
  ungroup() |>
  mutate(syph_diag_rate = ifelse(denominator == 0, 0, cases / denominator))

## CDC reported diagnosis rate among pregnant women
syph_pregnant_yr <- import(here("data", "surveillance", "syphilis_yr_pregnant.csv")) |>
  as_tibble() |>
  mutate(year = as.numeric(Year), cases = as.numeric(Cases)) |>
  select(year, cases) |>
  mutate(year = factor(year, levels = c(
    2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023
  ))) |>
  left_join(pregnancies_yr_race |>
              group_by(year) |>
              summarise(Pop = sum(total)),
            by = "year") |>
  mutate(diag_rate = cases / Pop) |>
  select(year, diag_rate)

## detected syphilis rate among pregnant women
detected_pregnant_yr <- pregnancies_yr_race |>
  left_join(detected_yr_race, by = c('year', 'race')) |>
  group_by(year) |>
  summarise(Pop = sum(Pop_mean), Detected = sum(Detected)) |>
  mutate(detected_rate = Detected / Pop) |>
  select(year, detected_rate)

## Reported syphilis diagnosis rates among stillbirths
cs_cause_yr <- import(here("data", "surveillance", "CS_stillbirth.csv")) |>
  as_tibble() |>
  filter(`Vital Status` == "Stillbirth") |>
  select(year = Year, cases = Cases) |>
  order_levels(c("year")) |>
  drop_na() |>
  left_join(pregnancies_yr_race |>
              group_by(year) |>
              summarise(fetal_deaths = sum(fetal_deaths)),
            by = "year") |>
  mutate(CS_prop = cases / fetal_deaths) |>
  select(year, CS_prop)


# 1 Fig 1: YEAR: estimated in S1 and S2, diagnoses ------------------------

pregnant_compare = P_pregnant_yr_S1 |>
  full_join(P_pregnant_yr_S2,
            by = "year",
            suffix = c("_S1", "_S2")) |>
  full_join(detected_pregnant_yr, by = 'year') |>
  full_join(syph_pregnant_yr, by = "year") |>
  full_join(
    syph_diag_yr_race |>
      group_by(year) |>
      summarise(
        cases = sum(cases),
        denominator = sum(denominator)
      ) |>
      mutate(syph_diag_rate = cases / denominator),
    by = "year"
  ) |>
  select(-cases, -denominator) |>
  mutate(across(-year, ~ .x * 100000)) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses among women of reproductive age")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses among women of reproductive age", group = 1)) +
  geom_point(aes(y = diag_rate, color = "Diagnoses among pregnant women")) +
  geom_line(aes(y = diag_rate, color = "Diagnoses among pregnant women", group = 1)) +
  geom_point(aes(y = detected_rate, color = "Detected syphilis in birth certificate")) +
  geom_line(aes(y = detected_rate, color = "Detected syphilis in birth certificate", group = 1)) +
  geom_point(aes(y = P_mean_S1, color = "Estimated prevalence in Scenario 1")) +
  geom_line(aes(y = P_mean_S1, color = "Estimated prevalence in Scenario 1", group = 1)) +
  geom_linerange(
    aes(ymin = P_lower_S1, ymax = P_upper_S1, color = "Estimated prevalence in Scenario 1"),
    linewidth = 0.5
  ) +
  geom_point(aes(
    x = as.numeric(year) + 0.05,
    y = P_mean_S2,
    color = "Estimated prevalence in Scenario 2"
  )) +
  geom_line(
    aes(
      x = as.numeric(year) + 0.05,
      y = P_mean_S2,
      color = "Estimated prevalence in Scenario 2",
      group = 1
    )
  ) +
  geom_linerange(
    aes(
      x = as.numeric(year) + 0.05,
      ymin = P_lower_S2,
      ymax = P_upper_S2,
      color = "Estimated prevalence in Scenario 2"
    ),
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Diagnoses among women of reproductive age" = "#107c10",
      "Diagnoses among pregnant women" = diagnoses_color,
      "Detected syphilis in birth certificate" = detected_color,
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2
    ),
    breaks = c(
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2",
      "Detected syphilis in birth certificate",
      "Diagnoses among pregnant women",
      "Diagnoses among women of reproductive age"
    )
  ) +
  labs(x = "Year", y = "Syphilis cases per 100,000 women", color = NULL) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    strip.background = element_rect(color = NULL, fill = NA),
    axis.text.x = element_text(
      size = size_fig - 2.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = size_fig),
    axis.title.x = element_text(size = size_fig),
    axis.title.y = element_text(size = size_fig),
    strip.text = element_text(size = size_fig),
    legend.text = element_text(size = size_fig)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

pregnant_compare

ggsave(
  plot = pregnant_compare,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "1_diagnoses_compare_main.svg"
  ),
  width = 12,
  height = 9,
  dpi = 800
)

# 2 Prevalence in Two S and Detected --------------------------------------

P_yr_race_S1 = fit_yr_race_S1 |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  ungroup()

P_yr_race_S2 = fit_yr_race_S2 |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  ungroup()

P_detected_yr_race = P_yr_race_S1 |>
  bind_rows(P_yr_race_S2, .id = "type") |>
  left_join(detected_rate_yr_race, by = c("year", "race")) |>
  mutate(across(c(mean:detected_rate), ~ . * 100000)) |>
  mutate(
    type = ifelse(
      type == 1,
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2"
    )
  ) |>
  ggplot(aes(x = year, color = type)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.5,
                 position = position_dodge2(width = 0.5)) +
  geom_point(aes(y = detected_rate, color = "Detected syphilis in birth certificate")) +
  
  geom_line(aes(y = detected_rate, color = "Detected syphilis in birth certificate", group = 1)) +
  geom_line(aes(y = mean, color = type, group = type), position = position_dodge(width = 0.5)) +
  
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap( ~ race, ncol = 3, scales = "free") +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2,
      "Detected syphilis in birth certificate" = detected_color
    ),
    breaks = c(
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2",
      "Detected syphilis in birth certificate"
    )
  ) +
  labs(x = "Year", y = "Syphilis cases per 100,000 live births", color = NULL) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    strip.background = element_rect(color = NULL, fill = NA),
    axis.text.x = element_text(
      size = size_fig - 2.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = size_fig),
    axis.title.x = element_text(size = size_fig),
    axis.title.y = element_text(size = size_fig),
    strip.text = element_text(size = size_fig),
    legend.text = element_text(size = size_fig),
    panel.spacing.x = unit(2, "lines"),
    panel.spacing.y = unit(2, "lines")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

P_detected_yr_race

# index calculation
index_df_S1 <- fit_yr_race_S1 |>
  left_join(P_pregnant_yr_S1, by = "year") |>
  select(year, race, P, .chain, .iteration, wt, P_mean) |>
  group_nest(year, .chain, .iteration) |>
  mutate(index = map_dbl(data, ~ {
    .x |>
      mutate(P_diff = abs(P - P_mean),
             index = 100 * sum(P_diff * wt) / P_mean) |>
      pull(index) |>
      mean()
  })) |>
  select(year, index) |>
  group_by(year) |>
  summarise(
    mean_index = mean(index),
    lower_index = quantile(index, probs = 0.025),
    upper_index = quantile(index, probs = 0.975)
  ) |>
  ungroup()

index_df_S2 <- fit_yr_race_S2 |>
  left_join(P_pregnant_yr_S2, by = "year") |>
  select(year, race, P, .chain, .iteration, wt, P_mean) |>
  group_nest(year, .chain, .iteration) |>
  mutate(index = map_dbl(data, ~ {
    .x |>
      mutate(P_diff = abs(P - P_mean),
             index = 100 * sum(P_diff * wt) / P_mean) |>
      pull(index) |>
      mean()
  })) |>
  select(year, index) |>
  group_by(year) |>
  summarise(
    mean_index = mean(index),
    lower_index = quantile(index, probs = 0.025),
    upper_index = quantile(index, probs = 0.975)
  ) |>
  ungroup()

index_df_detected <- detected_rate_yr_race |>
  left_join(fit_yr_race_S1 |> select(year, race, wt) |> distinct(),
            by = c("year", "race")) |>
  left_join(
    detected_rate_yr_race |>
      group_by(year) |>
      summarise(rate_mean = mean(detected_rate)) |>
      ungroup(),
    by = "year"
  ) |>
  group_nest(year) |>
  mutate(index = map_dbl(data, ~ {
    .x |>
      mutate(
        detected_diff = abs(detected_rate - rate_mean),
        index = 100 * sum(detected_diff * wt) / rate_mean
      ) |>
      pull(index) |>
      mean()
  })) |>
  select(year, index)

## index plot
index_df <- index_df_S1 |>
  bind_rows(index_df_S2,
            index_df_detected |> rename(mean_index = index),
            .id = "type") |>
  mutate(type = factor(
    case_when(
      type == "1" ~ "Estimated prevalence in Scenario 1",
      type == "2" ~ "Estimated prevalence in Scenario 2",
      type == "3" ~ "Detected syphilis in birth certificate"
    ),
    levels = c(
      "Estimated prevalence in Scenario 1",
      "Detected syphilis in birth certificate",
      "Estimated prevalence in Scenario 2"
    )
  ))

disparity_index = index_df |>
  ggplot(aes(x = year, y = mean_index, color = type)) +
  geom_point(position = position_dodge2(width = 0.25)) +
  geom_line(aes(group = type), position = position_dodge2(width = 0.25)) +
  geom_linerange(aes(ymin = lower_index, ymax = upper_index),
                 position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2,
      "Detected syphilis in birth certificate" = detected_color
    ),
    breaks = c(
      "Estimated prevalence in Scenario 1",
      "Detected syphilis in birth certificate",
      "Estimated prevalence in Scenario 2"
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = size_fig - 2.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = size_fig),
    axis.title.x = element_text(size = size_fig),
    axis.title.y = element_text(size = size_fig),
    plot.title = element_text(size = size_fig)
  )

disparity_index

ggsave(
  plot = ggdraw(P_detected_yr_race) +
    draw_plot(disparity_index + panel_border(), .58, .09, .4, .3),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "2_pre_yr_race_detected.svg"
  ),
  width = 12,
  height = 9,
  dpi = 800
)


# 3 Infer stillbirths -----------------------------------------------------

stillbirth_plot = P_pregnant_yr_S1 |>
  left_join(P_pregnant_yr_S2,
            by = "year",
            suffix = c("_S1", "_S2")) |>
  left_join(cs_cause_yr, by = "year") |>
  mutate(across(-year, ~ .x * 100000)) |>
  ggplot(aes()) +
  
  geom_point(aes(x = year, y = P_mean_S1, color = "Estimated prevalence in Scenario 1")) +
  geom_line(aes(
    x = year,
    y = P_mean_S1,
    color = "Estimated prevalence in Scenario 1",
    group = 1
  )) +
  geom_linerange(
    aes(
      x = year,
      ymin = P_lower_S1,
      ymax = P_upper_S1,
      color = "Estimated prevalence in Scenario 1"
    ),
    linewidth = 0.5
  ) +
  
  geom_point(aes(x = as.numeric(year)+0.05, y = P_mean_S2, color = "Estimated prevalence in Scenario 2")) +
  geom_line(aes(
    x = as.numeric(year)+0.05,
    y = P_mean_S2,
    color = "Estimated prevalence in Scenario 2",
    group = 1
  )) +
  geom_linerange(
    aes(
      x = as.numeric(year)+0.05,
      ymin = P_lower_S2,
      ymax = P_upper_S1,
      color = "Estimated prevalence in Scenario 2"
    ),
    linewidth = 0.5
  ) +
  
  geom_point(aes(x = year, y = CS_prop, color = "Stillbirths")) +
  geom_line(aes(
    x = year,
    y = CS_prop,
    color = "Stillbirths",
    group = 1
  )) +
  
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Syphilis cases per 100,000 births", color = NULL) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    strip.background = element_rect(color = NULL, fill = NA),
    axis.text.x = element_text(
      size = size_fig - 2.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = size_fig),
    axis.title.x = element_text(size = size_fig),
    axis.title.y = element_text(size = size_fig),
    strip.text = element_text(size = size_fig),
    legend.text = element_text(size = size_fig)
  ) +
  # scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 2" = prevalence_color_S2,
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Stillbirths" = "#5c005c"
    ),
    breaks = c(
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2",
      "Stillbirths"
    ),
    labels = c(
      "Estimated prevalence in Scenario 1" = "Estimated prevalence among live births in Scenario 1",
      "Estimated prevalence in Scenario 2" = "Estimated prevalence among live births in Scenario 2",
      "Stillbirths" = "Prevalence based on reported stillbirths attributable to congenital syphilis among stillbirths"
    )
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

stillbirth_plot

ggsave(
  plot = stillbirth_plot,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3_stillbirths.svg"
  ),
  width = 12,
  height = 9,
  dpi = 300
)
