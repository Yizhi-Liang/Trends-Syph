

rm(list = ls())

pacman::p_load(tidyverse,
               ggplot2,
               here,
               rio,
               tidybayes)

theme_set(theme_tidybayes() + cowplot::panel_border())
source(here("code", "0_functions", "order_levels.R"))

# 0 load data -----------------------------------------------------------------------------------------------------

fit_yr_race = import(here("output", "stan_output", "3_SA2.rds")) |>
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

# 1 prevalence by year and race -----------------------------------------------------------------------------------

prev_yr_SA2 = fit_yr_race |> 
  group_by(id, year, race) |>
  summarise(
    P_mean = mean(P),
    lower = quantile(P, 0.025),
    upper = quantile(P, 0.975)
  ) |> 
  mutate(across(P_mean:upper, ~.x * 1e5)) |> 
  pivot_wider(
    id_cols = c(year, race),
    names_from = id,
    values_from = c(P_mean, lower, upper)
  )

prev_yr_SA2 |> 
  export(here("output", "table", "1_prev_yr_SA2.csv"))

# prev_yr_SA2 |> 
#   export(here("output", "table", "1_prev_yr_SA2.xlsx"))

## plot

P_yr_race_SA2 = fit_yr_race |> 
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

ggsave(
  plot = P_yr_race_SA2,
  filename = here(
    "output",
    "figure",
    "est_prevalence",
    "9_prevalence_yr_race_SA2.pdf"
  ),
  width = 4*3.5,
  height = 3*3.5,
  dpi = 1200
)
