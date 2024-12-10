rm(list = ls())

pacman::p_load(tidyverse, ggplot2, here, rio, tidybayes, cowplot)

theme_set(theme_half_open())
source(here("code", "0_functions", "order_levels.R"))

size_fig <- 18

# color set
detected_color <- "#40a9ff"
prevalence_color_S1 <- "#a4262c"
prevalence_color_SA <- "#ff80ab"
prevalence_color_S2 <- "#ef6c00"
diagnoses_color <- "#ffb900"

# 0 load data -----------------------------------------------------------------------------------------------------

pregnancies_yr_race <- import(here("data", "parameter", "4_pregnancies.csv")) |>
  order_levels(c("year", "race")) |>
  select(year, race, Pop_mean = births, fetal_deaths, total) |>
  group_by(year) |>
  mutate(wt = Pop_mean / sum(Pop_mean)) |>
  ungroup()

fit_yr_race_S1 <- import(here("output", "stan_output", "1_S1.rds"), trust = TRUE) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

fit_yr_race_S2 <- import(here("output", "stan_output", "2_S2.rds"), trust = TRUE) |>
  as_tibble() |>
  order_levels(c("year", "race")) |>
  left_join(pregnancies_yr_race, by = c("year", "race"))

# 1 temporal trend ------------------------------------------------------------------------------------------------

# theme
p_theme <- theme(
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

########################################################################################
## data processing
P_time_func <- function(data) {
  result_df <- data |>
    group_by(race, .chain, .iteration) |>
    mutate(P = P * 100000) |>
    mutate(
      P_diff = P - P[year == "2016"],
      P_rr = P / P[year == "2016"]
    ) |>
    ungroup() |>
    filter(year != "2016")

  return(result_df)
}

P_time_trend_yr_race_S1 <- P_time_func(fit_yr_race_S1)
P_time_trend_yr_race_S2 <- P_time_func(fit_yr_race_S2)
P_time_trend_yr_race <- bind_rows(P_time_trend_yr_race_S1, P_time_trend_yr_race_S2, .id = "type") |>
  mutate(
    type = ifelse(
      type == 1,
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2"
    )
  )

########################################################################################
## absolute difference
P_time_trend_yr_race_diff <- P_time_trend_yr_race |>
  group_by(type, year, race) |>
  summarise(
    mean = median(P_diff),
    lower = quantile(P_diff, probs = 0.025),
    upper = quantile(P_diff, probs = 0.975)
  ) |>
  ungroup() |>
  filter(year != "2016")

time_trend_diff <- P_time_trend_yr_race_diff |>
  ggplot(aes(x = year, color = type)) +
  geom_point(
    data = P_time_trend_yr_race_diff,
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = P_time_trend_yr_race_diff,
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#595959"
  ) +
  facet_wrap(~race, ncol = 3) +
  scale_y_continuous(labels = scales::label_number(1)) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2
    ),
    labels = c(
      "Scenario 1",
      "Scenario 2"
    )
  ) +
  labs(
    x = "Year",
    y = "Absolute difference per 100,000 live births",
    color = "Scenario"
  ) +
  p_theme

time_trend_diff

ggsave(
  plot = time_trend_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "1.1_time_trend_diff.svg"
  ),
  width = 12,
  height = 9,
  dpi = 800
)

########################################################################################
## time ratio
P_time_trend_yr_race_ratio <- P_time_trend_yr_race |>
  group_by(type, year, race) |>
  summarise(
    mean = median(P_rr),
    lower = quantile(P_rr, probs = 0.025),
    upper = quantile(P_rr, probs = 0.975)
  ) |>
  ungroup() |>
  filter(year != "2016")

time_trend_ratio <- P_time_trend_yr_race_ratio |>
  ggplot(aes(x = year, color = type)) +
  geom_point(
    data = P_time_trend_yr_race_ratio,
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = P_time_trend_yr_race_ratio,
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "#595959"
  ) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2
    ),
    labels = c(
      "Scenario 1",
      "Scenario 2"
    )
  ) +
  facet_wrap(~race, ncol = 3, scales = "free") +
  scale_y_log10(labels = scales::label_number(1)) +
  labs(
    x = "Year",
    y = "Prevalence ratio in log scale",
    color = "Scenario"
  ) +
  p_theme

time_trend_ratio

ggsave(
  plot = time_trend_ratio,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "1.2_time_trend_ratio.svg"
  ),
  width = 12,
  height = 9,
  dpi = 800
)

## plot together

ggsave(
  plot = cowplot::plot_grid(
  time_trend_diff +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  time_trend_ratio +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  rel_heights = c(1, 0.1)
),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "1_time_trend.svg"
  ),
  width = 20,
  height = 12,
  dpi = 800
)

# 2 racial/ethnic disparity ---------------------------------------------------------------------------------------
disparity_yr_race_func <- function(data) {
  data |>
    group_by(year, .chain, .iteration) |>
    mutate(P = P * 100000) |>
    mutate(
      P_diff = P - P[race == "White"],
      P_rr = P / P[race == "White"]
    ) |>
    ungroup() |>
    filter(race != "White")
}

disparity_yr_race_S1 <- disparity_yr_race_func(fit_yr_race_S1)
disparity_yr_race_S2 <- disparity_yr_race_func(fit_yr_race_S2)
disparity_yr_race <- bind_rows(disparity_yr_race_S1, disparity_yr_race_S2, .id = "type") |>
  mutate(
    type = ifelse(
      type == 1,
      "Estimated prevalence in Scenario 1",
      "Estimated prevalence in Scenario 2"
    )
  )

## absolute difference
disparity_diff <- disparity_yr_race |>
  group_by(type, year, race) |>
  summarise(
    mean = median(P_diff),
    lower = quantile(P_diff, probs = 0.025),
    upper = quantile(P_diff, probs = 0.975)
  ) |>
  ungroup() |>
  filter(race != "White")

disparity_race_diff <- disparity_diff |>
  ggplot(aes(x = year, color = type)) +
  geom_point(
    data = disparity_diff,
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = disparity_diff,
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#595959"
  ) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2
    ),
    labels = c(
      "Scenario 1",
      "Scenario 2"
    )
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  facet_wrap(~race, ncol = 3) +
  labs(
    x = "Year",
    y = "Absolute difference per 100,000 live births",
    color = "Scenario"
  ) +
  p_theme

disparity_race_diff

ggsave(
  plot = disparity_race_diff,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "2.1_disparity_diff.svg"
  ),
  width = 16.5,
  height = 12,
  dpi = 800
)

## ratio
disparity_ratio <- disparity_yr_race |>
  group_by(type, year, race) |>
  summarise(
    mean = median(P_rr),
    lower = quantile(P_rr, probs = 0.05),
    upper = quantile(P_rr, probs = 0.95)
  ) |>
  ungroup() |>
  filter(race != "White")

disparity_race_ratio <- disparity_ratio |>
  ggplot(aes(x = year, color = type)) +
  geom_point(
    data = disparity_ratio,
    aes(y = mean),
    position = position_dodge(0.9)
  ) +
  geom_linerange(
    data = disparity_ratio,
    aes(ymin = lower, ymax = upper),
    position = position_dodge(0.9),
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "#595959"
  ) +
  scale_color_manual(
    values = c(
      "Estimated prevalence in Scenario 1" = prevalence_color_S1,
      "Estimated prevalence in Scenario 2" = prevalence_color_S2
    ),
    labels = c(
      "Scenario 1",
      "Scenario 2"
    )
  ) +
  facet_wrap(~race, ncol = 3, scales = "free") +
  labs(
    x = "Year",
    y = "Prevalence ratio in log scale",
    color = "Scenario"
  ) +
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
    "2.2_disparity_ratio.svg"
  ),
  width = 16.5,
  height = 12,
  dpi = 800
)

## together

ggsave(
  plot = cowplot::plot_grid(
  disparity_race_diff +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  disparity_race_ratio +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ),
  label_size = 12,
  rel_heights = c(1, 0.1)
),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "2_disparity.svg"
  ),
  width = 20,
  height = 12,
  dpi = 800
)

# 3 Sensitivity analysis with priors ----------------------------------------------------------------------------------------

fit_yr_race_S1_SA <- import(here("output", "stan_output", "1_S1_SA.rds")) |>
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

P_yr_race_S1_SA <- fit_yr_race_S1_SA |>
  group_by(id, year, race) |>
  summarise(
    P_mean = mean(P),
    lower = quantile(P, 0.025),
    upper = quantile(P, 0.975)
  ) |>
  mutate(across(P_mean:upper, ~ .x * 1e5)) |>
  ungroup() |>
  ggplot(aes(x = year, color = id)) +
  geom_point(aes(y = P_mean, shape = id),
    position = position_dodge(width = 0.9)
  ) +
  geom_linerange(aes(ymin = lower, ymax = upper),
    linewidth = 0.5,
    position = position_dodge2(width = 0.9)
  ) +
  scale_shape_manual(values = c(16, 17, 1, 2, 3)) +
  ggsci::scale_color_nejm() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~race, ncol = 3, scales = "free") +
  labs(
    x = "Year",
    y = "Estimated syphilis prevalence per 100,000 live births",
    shape = NULL,
    color = NULL
  ) +
  p_theme +
  guides(
    color = "none", shape = "none"
  )

fit_yr_race_S2_SA <- import(here("output", "stan_output", "2_S2_SA.rds")) |>
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

P_yr_race_S2_SA <- fit_yr_race_S2_SA |>
  group_by(id, year, race) |>
  summarise(
    P_mean = mean(P),
    lower = quantile(P, 0.025),
    upper = quantile(P, 0.975)
  ) |>
  mutate(across(P_mean:upper, ~ .x * 1e5)) |>
  ungroup() |>
  ggplot(aes(x = year, color = id)) +
  geom_point(aes(y = P_mean, shape = id),
    position = position_dodge(width = 0.9)
  ) +
  geom_linerange(aes(ymin = lower, ymax = upper),
    linewidth = 0.5,
    position = position_dodge2(width = 0.9)
  ) +
  scale_shape_manual(values = c(16, 17, 1, 2, 3)) +
  ggsci::scale_color_nejm() +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~race, ncol = 3, scales = "free") +
  labs(
    x = "Year",
    y = "Estimated syphilis prevalence per 100,000 live births",
    shape = "Weakly informative priors for prevalence",
  ) +
  p_theme +
  guides(
    color = "none",
    shape = guide_legend(nrow = 1, byrow = TRUE)
  )

ggsave(
  plot = cowplot::plot_grid(
    P_yr_race_S1_SA, P_yr_race_S2_SA,
    labels = c("Scenario 1", "Scenario 2"),
    ncol = 1
  ),
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "appendix",
    "3_prevalence_yr_race_SA.svg"
  ),
  width = 4 * 3.5,
  height = 6 * 3.5,
  dpi = 800
)
