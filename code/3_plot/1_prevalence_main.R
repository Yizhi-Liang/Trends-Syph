

rm(list = ls())

pacman::p_load(tidyverse,
               ggplot2,
               here,
               rio,
               tidybayes)

theme_set(theme_tidybayes() + cowplot::panel_border())
source(here("code", "0_functions", "order_levels.R"))

# 0 load data -----------------------------------------------------------------------------------------------------

pregnancies_yr_race = import(here("data",
                                  "parameter",
                                  "4_pregnancies.csv")) |>
  order_levels(c("year", "race")) |>
  select(year, race, Pop_mean = births, fetal_deaths, total) |>
  group_by(year) |>
  mutate(wt = Pop_mean / sum(Pop_mean)) |>
  ungroup()

fit_yr_race = import(here("output", "stan_output", "1_main.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

# 1 prevalence by year and race -----------------------------------------------------------------------------------

## by year

P_yr = fit_yr_race |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  ungroup() |>
  left_join(pregnancies_yr_race |> select(year:Pop_mean),
            by = c("year", "race")) |>
  mutate(
    numerator = Pop_mean * mean,
    numerator_lower = Pop_mean * lower,
    numerator_upper = Pop_mean * upper
  ) |>
  group_by(year) |>
  summarise(
    Pop_mean = sum(Pop_mean),
    numerator = sum(numerator),
    numerator_lower = sum(numerator_lower),
    numerator_upper = sum(numerator_upper)
  ) |>
  mutate(
    P = 1e5 * numerator / Pop_mean,
    lower = 1e5 * numerator_lower / Pop_mean,
    upper = 1e5 * numerator_upper / Pop_mean
  ) |>
  select(year, P:upper)

P_yr |>
  export(here("output", "table", "1_prev_yr.csv"))

# by year and race
P_yr_race_manual = fit_yr_race |>
  group_by(year, race) |>
  summarise(
    mean = mean(P),
    lower = quantile(P, probs = 0.025),
    upper = quantile(P, probs = 0.975)
  ) |>
  mutate(across(c(mean, lower, upper), ~ .x * 1e5)) |>
  ungroup() |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = mean)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.5) +
  scale_y_continuous(limits = c(0, NA),
                     labels = scales::label_number(1)) +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  labs(x = "Year",
       y = "Estimated syphilis prevalence per 100,000 women",
       color = "Race/Ethnicity") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

P_yr_race_manual

ggsave(
  plot = P_yr_race_manual,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "1_prevalence_yr_race_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

# 2 compare with diagnoses among general women by year and race ---------------------------------------------------

positivities_yr_race = import(here("data",
                                   "parameter",
                                   "0_positivities.csv")) |>
  order_levels(c("year", "race"))

positive_rate = pregnancies_yr_race |>
  left_join(positivities_yr_race, by = c('year', 'race')) |>
  mutate(positivity_rate = Positivities / Pop_mean) |>
  select(year, race, positivity_rate)

syph_diag_yr_race <- import(here("data",
                                 "surveillance",
                                 "atlas_syphilis_yr_race.csv")) |>
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
    type = case_when(str_detect(type, "Primary and Secondary") ~ "PS",
                     TRUE ~ "Early"),
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

general_compare = syph_diag_yr_race |>
  group_by(year, race) |>
  summarise(cases = sum(cases),
            denominator = sum(denominator)) |>
  mutate(syph_diag_rate = ifelse(denominator == 0, 0, cases / denominator)) |>
  ungroup() |>
  left_join(positive_rate, by = c("year", "race")) |>
  select(year, race, syph_diag_rate, positivity_rate) |>
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
    syph_diag_rate = syph_diag_rate * 100000,
    positivity_rate = positivity_rate * 100000,
    mean = mean * 100000,
    lower = lower * 100000,
    upper = upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_point(aes(y = mean, color = "Estimated Prevalance")) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated Prevalance"),
                 linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726",
      "Estimated Prevalance" = "#bc3d29"
    )
  ) +
  facet_wrap(~ race, scales = "free") +
  theme_bw() +
  scale_y_log10(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women in log scale",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

general_compare

ggsave(
  plot = general_compare,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "2_general_compare_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

syph_diag_yr_race |>
  group_by(year, race) |>
  summarise(cases = sum(cases),
            denominator = sum(denominator)) |>
  mutate(syph_diag_rate = ifelse(denominator == 0, 0, cases / denominator)) |>
  ungroup() |>
  left_join(positive_rate, by = c("year", "race")) |>
  select(year, race, syph_diag_rate, positivity_rate) |>
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
    syph_diag_rate = syph_diag_rate * 100000,
    positivity_rate = positivity_rate * 100000,
    mean = mean * 100000,
    lower = lower * 100000,
    upper = upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_point(aes(y = mean, color = "Estimated Prevalance")) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated Prevalance"),
                 linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726",
      "Estimated Prevalance" = "#bc3d29"
    )
  ) +
  facet_wrap(~ race, scales = "free") +
  theme_bw() +
  scale_y_continuous(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  plot = last_plot(),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "2_general_compare_nolog_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

# 3 compare with diagnoses among pregnant women by year -----------------------------------------------------------

syph_pregnant <- import(here("data", "surveillance",
                             "syphilis_yr_pregnant.csv")) |>
  as_tibble() |>
  mutate(year = as.numeric(Year),
         cases = as.numeric(Cases)) |>
  select(year, cases) |>
  mutate(year = factor(year,
                       levels = c(2016, 2017, 2018,
                                  2019, 2020, 2021,
                                  2022))) |>
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

pregnant_compare = syph_pregnant |>
  left_join(positivies_pregnant_yr, by = 'year') |>
  left_join(P_pregnant_yr,
            by = "year") |>
  mutate(
    syph_diag_rate = diag_rate * 100000,
    positivity_rate = pos_rate * 100000,
    mean = P_mean * 100000,
    lower = P_lower * 100000,
    upper = P_upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses", group = 1)) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_line(aes(y = positivity_rate, color = "Positivities", group = 1)) +
  geom_point(aes(y = mean, color = "Estimated Prevalance")) +
  geom_line(aes(y = mean, color = "Estimated Prevalance", group = 1)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated Prevalance"),
                 linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726",
      "Estimated Prevalance" = "#bc3d29"
    )
  ) +
  theme_bw() +
  scale_y_log10(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women in log scale",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

pregnant_compare

ggsave(
  plot = pregnant_compare,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3_pregnant_compare_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

syph_pregnant |>
  left_join(positivies_pregnant_yr, by = 'year') |>
  left_join(P_pregnant_yr,
            by = "year") |>
  mutate(
    syph_diag_rate = diag_rate * 100000,
    positivity_rate = pos_rate * 100000,
    mean = P_mean * 100000,
    lower = P_lower * 100000,
    upper = P_upper * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses", group = 1)) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_line(aes(y = positivity_rate, color = "Positivities", group = 1)) +
  geom_point(aes(y = mean, color = "Estimated Prevalance")) +
  geom_line(aes(y = mean, color = "Estimated Prevalance", group = 1)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = "Estimated Prevalance"),
                 linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726",
      "Estimated Prevalance" = "#bc3d29"
    )
  ) +
  theme_bw() +
  scale_y_continuous(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  plot = last_plot(),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3_pregnant_compare_nolog_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

syph_pregnant |>
  left_join(positivies_pregnant_yr, by = 'year') |>
  mutate(
    syph_diag_rate = diag_rate * 100000,
    positivity_rate = pos_rate * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses", group = 1)) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_line(aes(y = positivity_rate, color = "Positivities", group = 1)) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726"
    )
  ) +
  theme_bw() +
  scale_y_log10(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women in log scale",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  plot = last_plot(),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3.1_pregnant_compare.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

syph_pregnant |>
  left_join(positivies_pregnant_yr, by = 'year') |>
  mutate(
    syph_diag_rate = diag_rate * 100000,
    positivity_rate = pos_rate * 100000
  ) |>
  ggplot(aes(x = year)) +
  geom_point(aes(y = syph_diag_rate, color = "Diagnoses")) +
  geom_line(aes(y = syph_diag_rate, color = "Diagnoses", group = 1)) +
  geom_point(aes(y = positivity_rate, color = "Positivities")) +
  geom_line(aes(y = positivity_rate, color = "Positivities", group = 1)) +
  scale_color_manual(
    values = c(
      "Diagnoses" = "#0073b5",
      "Positivities" = "#e28726"
    )
  ) +
  theme_bw() +
  scale_y_continuous(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Cases per 100,000 women",
       color = NULL) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  plot = last_plot(),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "3.1_pregnant_compare_nolog.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)


# 4 temporal change with 2014 -------------------------------------------------------------------------------------

P_time_trend_yr_race = fit_yr_race |>
  group_by(race, .chain, .iteration) |>
  mutate(P = P * 100000) |>
  mutate(P_diff = P - P[year == "2014"],
         P_rr = P / P[year == "2014"],
         P_rel_ratio = P_diff / P[year == "2014"]) |>
  ungroup() |>
  filter(year != "2014")

## absolute difference
time_trend_diff = P_time_trend_yr_race |>
  ggplot(aes(x = year)) +
  
  geom_point(
    data = P_time_trend_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_diff),
        lower = quantile(P_diff, probs = 0.025),
        upper = quantile(P_diff, probs = 0.975)
      ) |>
      ungroup() |>
      filter(year != "2014"),
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = P_time_trend_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_diff),
        lower = quantile(P_diff, probs = 0.025),
        upper = quantile(P_diff, probs = 0.975)
      ) |>
      ungroup() |>
      filter(year != "2014"),
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "#595959") +
  facet_wrap(~ race, ncol = 3) +
  scale_y_continuous(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Absolute difference per 100,000 women",
       fill = "Race/Ethnicity") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

time_trend_diff

ggsave(
  plot = time_trend_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "4.1_time_trend_diff_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

time_trend_ratio = P_time_trend_yr_race |>
  ggplot(aes(x = year)) +
  
  geom_point(
    data = P_time_trend_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_rr),
        lower = quantile(P_rr, probs = 0.05),
        upper = quantile(P_rr, probs = 0.95)
      ) |>
      ungroup() |>
      filter(year != "2014"),
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = P_time_trend_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_rr),
        lower = quantile(P_rr, probs = 0.025),
        upper = quantile(P_rr, probs = 0.975)
      ) |>
      ungroup() |>
      filter(year != "2014"),
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "#595959") +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  scale_y_log10(labels = scales::label_number(1)) +
  labs(x = "Year",
       y = "Prevalence ratio in log scale") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

time_trend_ratio

ggsave(
  plot = time_trend_ratio,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "4.2_time_trend_ratio_main.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

## plot together
cowplot::plot_grid(
  time_trend_diff +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  time_trend_ratio +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  rel_heights = c(1, 0.1)
)

ggsave(
  plot = last_plot(),
  filename = here("output",
                  "figure",
                  "est_prevalence",
                  "4_time_trend_main.pdf"),
  width = 16.5,
  height = 12,
  dpi = 1200
)

# 5 racial ratio with white ---------------------------------------------------------------------------------------

disparity_yr_race = fit_yr_race |>
  group_by(year, .chain, .iteration) |>
  mutate(P = P * 100000) |>
  mutate(
    P_diff = P - P[race == "White"],
    P_rr = P / P[race == "White"],
    P_relative_ratio = P_diff / P[race == "White"]
  ) |>
  ungroup() |>
  filter(race != "White")

## absolute difference
disparity_race_diff = disparity_yr_race |>
  ggplot(aes(x = year)) +
  
  geom_point(
    data = disparity_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_diff),
        lower = quantile(P_diff, probs = 0.025),
        upper = quantile(P_diff, probs = 0.975)
      ) |>
      ungroup() |>
      filter(race != "Non-H White"),
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = disparity_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_diff),
        lower = quantile(P_diff, probs = 0.025),
        upper = quantile(P_diff, probs = 0.975)
      ) |>
      ungroup() |>
      filter(race != "Non-H White"),
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "#595959") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  facet_wrap(~ race, ncol = 3) +
  labs(x = "Year",
       y = "Absolute difference per 100,000 women",
       fill = "Race/Ethnicity") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

disparity_race_diff

ggsave(
  plot = disparity_race_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "5.1_disparity_diff_main.pdf"
  ),
  width = 16.5,
  height = 12,
  dpi = 1200
)

## ratio
disparity_race_ratio = disparity_yr_race |>
  ggplot(aes(x = year)) +
  geom_point(
    data = disparity_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_rr),
        lower = quantile(P_rr, probs = 0.05),
        upper = quantile(P_rr, probs = 0.95)
      ) |>
      ungroup() |>
      filter(race != "Non-H White"),
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = disparity_yr_race |>
      group_by(year, race) |>
      summarise(
        mean = mean(P_rr),
        lower = quantile(P_rr, probs = 0.025),
        upper = quantile(P_rr, probs = 0.975)
      ) |>
      ungroup() |>
      filter(race != "Non-H White"),
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "#595959") +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  labs(x = "Year",
       y = "Prevalence ratio in log scale",
       fill = "Race/Ethnicity") +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "#595959",
                                    fill = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

disparity_race_ratio

ggsave(
  plot = disparity_race_ratio,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "5.2_disparity_ratio_main.pdf"
  ),
  width = 16.5,
  height = 12,
  dpi = 1200
)

## together
cowplot::plot_grid(
  disparity_race_diff +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  disparity_race_ratio +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  label_size = 12,
  rel_heights = c(1, 0.1)
)

ggsave(
  plot = last_plot(),
  filename = here("output",
                  "figure",
                  "est_prevalence",
                  "5_disparity_main.pdf"),
  width = 16,
  height = 12,
  dpi = 1200
)

# 6 index of disparity --------------------------------------------------------------------------------------------

P_mean_yr <- fit_yr_race |>
  group_by(year) |>
  summarise(P_mean = mean(P)) |>
  ungroup()

index_df <- fit_yr_race |>
  left_join(P_mean_yr, by = "year") |>
  select(year, race, P, .chain, .iteration, wt, P_mean) |>
  group_nest(year, .chain, .iteration) |>
  mutate(index = map_dbl(data,
                         ~ {
                           .x |>
                             mutate(P_diff = abs(P - P_mean),
                                    index = 100 * sum(P_diff * wt) / P_mean) |>
                             pull(index) |>
                             mean()
                         })) |>
  select(year, index)

disparity_index = index_df |>
  group_by(year) |>
  summarise(
    mean_index = mean(index),
    lower_index = quantile(index, probs = 0.025),
    upper_index = quantile(index, probs = 0.975)
  ) |>
  ungroup() |>
  ggplot(aes(x = year, y = mean_index)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_linerange(aes(ymin = lower_index, ymax = upper_index),
                 linewidth = 0.5) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(x = "Year",
       y = "Index of disparity") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

disparity_index

ggsave(
  plot = disparity_index,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "6_index_disparity_main.pdf"
  ),
  width = 10,
  height = 7.5,
  dpi = 1200
)

# 7 infer the stillbirths -----------------------------------------------------------------------------------------

pregnancies_yr = pregnancies_yr_race |> 
  group_by(year) |> 
  summarise(fetal_deaths = sum(fetal_deaths))
  
cs_cause_yr <- import(here("data",
                           "surveillance",
                           "CS_stillbirth.csv")) |>
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

stillbirth_plot = prevalence_yr |> 
  left_join(cs_cause_yr, by = "year") |> 
  mutate(across(
    P:CS_prop,
    ~ .x * 100000
  )) |> 
  pivot_longer(
    cols = c("P", "CS_prop"),
    names_to = "variable",
    values_to = "value"
  ) |> 
  mutate(type = ifelse(variable == "P", "Live births (Estimated)", "Stillbirths")) |>
  ggplot(aes(x = year, y = value, shape = type)) +

  geom_linerange(aes(ymin = P_lower, ymax = P_upper), 
                 linewidth = 0.5) +
  geom_point() +
  geom_line(aes(y = value, group = type, linetype = type)) +
  
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Syphilis prevalence per 100,000 women",
    linetype = NULL,
    shape = NULL
  ) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  scale_shape_manual(values = c(16, 4))

ggsave(
  plot = stillbirth_plot,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "7_stillbirths_main.pdf"
  ),
  width = 10,
  height = 7.5,
  dpi = 1200
)

prevalence_yr |> 
  left_join(cs_cause_yr, by = "year") |> 
  mutate(
    ratio = CS_prop / P,
    lower_ratio = CS_prop / P_upper,
    upper_ratio = CS_prop / P_lower,
    CS_prop_100k = CS_prop * 100000
  )


