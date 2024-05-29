
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

fit_yr_race_SA1 = import(here("output", "stan_output", "2_SA1.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

# 1 prevalence by year and race -----------------------------------------------------------------------------------

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

P_yr_race_manual =  P_main_yr_race |> 
  bind_rows(P_sa_yr_race, .id = "type") |> 
  mutate(type = ifelse(type == "1", "Main Analysis", "Sensitivity Analysis")) |>
  mutate(across(c(mean, lower, upper), ~ . * 100000)) |> 
  ggplot(aes(x = year)) +
  geom_point(aes(y = mean, shape = type), 
             position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.5,
                 position = position_dodge2(width = 0.5)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = "Year",
    y = "Estimated syphilis prevalence per 100,000 women",
    shape = NULL
  ) +
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
  guides(color = "none",
         shape = guide_legend(nrow = 1, byrow = TRUE))

P_yr_race_manual

ggsave(
  plot = P_yr_race_manual,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "8_prevalence_yr_race_SA1.pdf"
  ),
  width = 12,
  height = 9,
  dpi = 1200
)

