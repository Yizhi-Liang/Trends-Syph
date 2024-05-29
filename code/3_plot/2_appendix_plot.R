
rm(list = ls())

pacman::p_load(tidyverse, ggplot2, here, rio, tidybayes, cowplot)

theme_set(theme_half_open() + panel_border())
source(here("code", "0_functions", "order_levels.R"))

size_fig = 18

# color set
prevalence_color_main = "#a4262c"

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

# 1 temporal trend ------------------------------------------------------------------------------------------------

# theme
p_theme = theme(
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
)


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
  p_theme

time_trend_diff

ggsave(
  plot = time_trend_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "1.1_time_trend_diff_main.pdf"
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
  P_theme

time_trend_ratio

ggsave(
  plot = time_trend_ratio,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "1.2_time_trend_ratio_main.pdf"
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
                  "appendix",
                  "1_time_trend_main.pdf"),
  width = 20,
  height = 12,
  dpi = 1200
)


# 2 racial/ethnic disparity ---------------------------------------------------------------------------------------

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
  p_theme

disparity_race_diff

ggsave(
  plot = disparity_race_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "2.1_disparity_diff_main.pdf"
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
  p_theme

disparity_race_ratio

ggsave(
  plot = disparity_race_ratio,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "2.2_disparity_ratio_main.pdf"
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
                  "appendix",
                  "2_disparity_main.pdf"),
  width = 20,
  height = 12,
  dpi = 1200
)

# 3 Sensitivity analysis 2 ----------------------------------------------------------------------------------------

fit_yr_race_SA2 = import(here("output", "stan_output", "3_SA2.rds")) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  mutate(id = factor(
    case_when(
      P_beta == 1 ~ "Beta(1, 1)",
      P_beta == 2 ~ "Beta(1, 2)",
      P_beta == 3 ~ "Beta(1, 3)",
      P_beta == 4 ~ "Beta(1, 4)",
      P_beta == 5 ~ "Beta(1, 5)"
    ),
    levels = c(
      "Beta(1, 5)",
      "Beta(1, 4)",
      "Beta(1, 3)",
      "Beta(1, 2)",
      "Beta(1, 1)"
    )
  ))

P_yr_race_SA2 = fit_yr_race_SA2 |> 
  group_by(id, year, race) |> 
  summarise(
    P_mean = mean(P),
    lower = quantile(P, 0.025),
    upper = quantile(P, 0.975)
  ) |> 
  mutate(across(P_mean:upper, ~.x * 1e5)) |> 
  ungroup() |>
  ggplot(aes(x = year, color = id)) +
  geom_point(aes(y = P_mean, shape = id),
             position = position_dodge(width = 0.9)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.5,
                 position = position_dodge2(width = 0.9)) +
  scale_shape_manual(values = c(16, 17, 1, 2, 3)) +
  ggsci::scale_color_nejm() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ race, ncol = 3, scales = "free") +
  labs(
    x = "Year",
    y = "Estimated syphilis prevalence per 100,000 women",
    shape = NULL
  ) +
  p_theme +
  guides(color = "none",
         shape = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  plot = P_yr_race_SA2,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "3_prevalence_yr_race_SA2.pdf"
  ),
  width = 4*3.5,
  height = 3*3.5,
  dpi = 1200
)
