


rm(list = ls())

pacman::p_load(tidyverse, ggplot2, here, rio, tidybayes, cowplot)

theme_set(theme_half_open())
source(here("code", "0_functions", "order_levels.R"))

size_fig = 18

# color set
positives_color = "#40a9ff"
prevalence_color_main = "#a4262c"
prevalence_color_SA = "#ff80ab"
diagnoses_color = "#ffb900"

# 0 load data -----------------------------------------------------------------------------------------------------

pregnancies_yr_race = import(here("data", "parameter", "4_pregnancies.csv")) |>
  order_levels(c("year", "race")) |>
  select(year, race, Pop_mean = births, fetal_deaths, total) |>
  group_by(year) |>
  mutate(wt = Pop_mean / sum(Pop_mean)) |>
  ungroup()

fit_yr_race = import(here("output", "stan_output", "1_main.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

positivities_yr_race = import(here("data", "parameter", "0_positivities.csv")) |>
  order_levels(c("year", "race"))

positive_rate = pregnancies_yr_race |>
  left_join(positivities_yr_race, by = c('year', 'race')) |>
  mutate(positivity_rate = Positivities / Pop_mean) |>
  select(year, race, positivity_rate)

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

syph_pregnant <- import(here("data", "surveillance", "syphilis_yr_pregnant.csv")) |>
  as_tibble() |>
  mutate(year = as.numeric(Year), cases = as.numeric(Cases)) |>
  select(year, cases) |>
  mutate(year = factor(year, levels = c(2016, 2017, 2018, 2019, 2020, 2021, 2022))) |>
  left_join(pregnancies_yr_race |>
              group_by(year) |>
              summarise(Pop = sum(total)),
            by = "year") |>
  mutate(diag_rate = cases / Pop) |>
  select(year, diag_rate)

positivies_pregnant_yr <- pregnancies_yr_race |>
  left_join(positivities_yr_race, by = c('year', 'race')) |>
  group_by(year) |>
  summarise(Pop = sum(Pop_mean),
            Positivities = sum(Positivities)) |>
  mutate(pos_rate = Positivities / Pop) |>
  select(year, pos_rate)

P_pregnant_yr <- fit_yr_race |>
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

pregnancies_yr = pregnancies_yr_race |>
  group_by(year) |>
  summarise(fetal_deaths = sum(fetal_deaths))

cs_cause_yr <- import(here("data", "surveillance", "CS_stillbirth.csv")) |>
  as_tibble() |>
  filter(`Vital Status` == "Stillbirth") |>
  select(year = Year, cases = Cases) |>
  order_levels(c("year")) |>
  drop_na() |>
  left_join(pregnancies_yr, by = "year") |>
  mutate(CS_prop = cases / fetal_deaths) |>
  select(year, CS_prop)

prevalence_yr = fit_yr_race |>
  group_by(year, race) |>
  summarise(
    P_mean = mean(P),
    P_lower = quantile(P, probs = 0.025),
    P_upper = quantile(P, probs = 0.975)
  ) |>
  left_join(pregnancies_yr_race |>
              select(year, race, live_births = Pop_mean),
            by = c("year", "race")) |>
  mutate(
    numerator = P_mean * live_births,
    lower_numerator = P_lower * live_births,
    upper_numerator = P_upper * live_births
  ) |>
  summarise(
    numerator = sum(numerator),
    lower_numerator = sum(lower_numerator),
    upper_numerator = sum(upper_numerator),
    denominator = sum(live_births)
  ) |>
  mutate(
    P = numerator / denominator,
    P_lower = lower_numerator / denominator,
    P_upper = upper_numerator / denominator
  ) |>
  select(year, P:P_upper)


# 1 Estimated prevalence and observed positives -------------------------------------------------------------------

pos_compare = positive_rate |>
  left_join(
    fit_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P),
        lower = quantile(P, probs = 0.025),
        upper = quantile(P, probs = 0.975)
      ) |>
      ungroup(),
    by = c("year", "race")
  ) |>
  mutate(
    positivity_rate = positivity_rate * 100000,
    mean = mean * 100000,
    lower = lower * 100000,
    upper = upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = mean, color = "Estimated prevalence")) +
  geom_point(aes(y = positivity_rate, color = "Birth certificate positivity")) +
  
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated prevalence"),
                 linewidth = 0.5) +
  geom_line(aes(y = mean, color = "Estimated prevalence", group = 1)) +
  geom_line(aes(y = positivity_rate, color = "Birth certificate positivity", group = 1)) +
  scale_color_manual(
    values = c(
      "Birth certificate positivity" = positives_color,
      "Estimated prevalence" = prevalence_color_main
    ),
    breaks = c("Estimated prevalence", "Birth certificate positivity")
  ) +
  facet_wrap(~ race, scales = "free") +
  labs(x = "Year", y = "Cases per 100,000 live births", color = NULL) +
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
  # order of the color legend
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

pos_compare

# disparity index
P_mean_yr <- fit_yr_race |>
  group_by(year) |>
  summarise(P_mean = mean(P)) |>
  ungroup()

## index calculation
index_df <- fit_yr_race |>
  left_join(P_mean_yr, by = "year") |>
  select(year, race, P, .chain, .iteration, wt, P_mean) |>
  group_nest(year, .chain, .iteration) |>
  mutate(index = map_dbl(data, ~ {
    .x |>
      mutate(P_diff = abs(P - P_mean),
             index = 100 * sum(P_diff * wt) / P_mean) |>
      pull(index) |>
      mean()
  })) |>
  select(year, index)

## index plot
disparity_index = index_df |>
  group_by(year) |>
  summarise(
    mean_index = mean(index),
    lower_index = quantile(index, probs = 0.025),
    upper_index = quantile(index, probs = 0.975)
  ) |>
  ungroup() |>
  ggplot(aes(x = year, y = mean_index, color = "black")) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_linerange(aes(ymin = lower_index, ymax = upper_index), linewidth = 0.5) +
  scale_color_manual(
    values = c("black" = "black"),
    labels = c("black" = "Index of disparity")
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(
      size = size_fig - 2.5,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = size_fig),
    axis.title.x = element_text(size = size_fig),
    axis.title.y = element_text(size = size_fig),
    plot.title = element_text(size = size_fig),
    legend.text = element_text(size = size_fig),
    legend.position = c(0.45, 0.1)
  ) +
  guides(color = guide_legend(title = NULL))

disparity_index

ggsave(
  plot = ggdraw(pos_compare) +
    draw_plot(disparity_index + panel_border(), .48, .09, .4, .3),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "1_positives_compare_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)


# 2 Estimated Prevalence, Positives, and Diagnoses among Pregnant women -------------------------------------------

pregnant_compare = P_pregnant_yr |>
  full_join(positivies_pregnant_yr, by = 'year') |>
  full_join(syph_pregnant, by = "year") |>
  full_join(
    syph_diag_yr_race |>
      group_by(year) |>
      summarise(
        cases = sum(cases),
        denominator = sum(denominator)
      ) |>
      mutate(syph_diag_rate = 100000 * cases / denominator),
    by = "year"
  ) |>
  mutate(
    syph_diag_rate_preg = diag_rate * 100000,
    positivity_rate = pos_rate * 100000,
    mean = P_mean * 100000,
    lower = P_lower * 100000,
    upper = P_upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses among women of reproductive age")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses among women of reproductive age", group = 1)) +
  geom_point(aes(y = syph_diag_rate_preg, color = "Diagnoses among pregnant women")) +
  geom_line(aes(y = syph_diag_rate_preg, color = "Diagnoses among pregnant women", group = 1)) +
  geom_point(aes(y = positivity_rate, color = "Birth certificate positivity")) +
  geom_line(aes(y = positivity_rate, color = "Birth certificate positivity", group = 1)) +
  geom_point(aes(y = mean, color = "Estimated prevalence")) +
  geom_line(aes(y = mean, color = "Estimated prevalence", group = 1)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated prevalence"),
                 linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Diagnoses among women of reproductive age" = "#107c10",
      "Diagnoses among pregnant women" = diagnoses_color,
      "Birth certificate positivity" = positives_color,
      "Estimated prevalence" = prevalence_color_main
    ),
    breaks = c(
      "Estimated prevalence",
      "Birth certificate positivity",
      "Diagnoses among pregnant women",
      "Diagnoses among women of reproductive age"
    )
  ) +
  labs(x = "Year", y = "Cases per 100,000 women", color = NULL) +
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
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

pregnant_compare

ggsave(
  plot = pregnant_compare,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "2_diagnoses_compare_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

# 3 Infer stillbirths ---------------------------------------------------------------------------------------------

stillbirth_plot = prevalence_yr |>
  left_join(cs_cause_yr, by = "year") |>
  mutate(across(P:CS_prop, ~ .x * 100000)) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = P, color = "Live births (Estimated)")) +
  geom_point(aes(y = CS_prop, color = "Stillbirths")) +
  geom_line(aes(y = P, color = "Live births (Estimated)", group = 1)) +
  geom_line(aes(y = CS_prop, color = "Stillbirths", group = 1)) +
  geom_linerange(aes(ymin = P_lower, ymax = P_upper, color = "Live births (Estimated)"),
                 linewidth = 0.5) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Year", y = "Syphilis prevalence per 100,000 women", color = NULL) +
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
      "Live births (Estimated)" = prevalence_color_main,
      "Stillbirths" = "#5c005c"
    ),
    labels = c(
      "Live births (Estimated)" = "Estimated prevalence among live births",
      "Stillbirths" = "Prevalence based on reported stillbirths attributable to congenital syphilis among stillbirths"
    )
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

stillbirth_plot

ggsave(
  plot = stillbirth_plot,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3_stillbirths_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

# 4 Sensitivity analysis 1 ----------------------------------------------------------------------------------------

fit_yr_race_SA1 = import(here("output", "stan_output", "2_SA1.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

P_main_yr_race = fit_yr_race |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  ungroup()

P_sa_yr_race = fit_yr_race_SA1 |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  ungroup()

P_yr_race_manual = P_main_yr_race |>
  bind_rows(P_sa_yr_race, .id = "type") |>
  left_join(positive_rate, by = c("year", "race")) |>
  mutate(across(c(mean:positivity_rate), ~ . * 100000)) |>
  mutate(type = ifelse(type == 1, "Main analysis", "Sensitivity analysis")) |>
  ggplot(aes(x = year, color = type)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.5,
                 position = position_dodge2(width = 0.9)) +
  geom_point(aes(y = positivity_rate, color = "Birth certificate positivity")) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  scale_color_manual(
    values = c(
      "Main analysis" = prevalence_color_main,
      "Sensitivity analysis" = prevalence_color_SA,
      "Birth certificate positivity" = positives_color
    ),
    breaks = c("Main analysis", "Sensitivity analysis", "Birth certificate positivity")
  ) +
  labs(x = "Year", y = "Estimated syphilis prevalence per 100,000 live births", color = NULL) +
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
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

P_yr_race_manual

ggsave(
  plot = P_yr_race_manual,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "4_prevalence_yr_race_SA1.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)
