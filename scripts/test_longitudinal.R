devtools::load_all("/Users/billyquilty/Documents/Work/synthpop")
library(tidyverse)

set.seed(42)

# --- Generate a synthetic population with viral load trajectories ---
n_individuals <- 200
n_timepoints <- 14

pop <- tibble(
  id = 1:n_individuals,
  age = round(rnorm(n_individuals, mean = 45, sd = 15)),
  sex = sample(c("M", "F"), n_individuals, replace = TRUE),
  vaccinated = sample(c(0, 1), n_individuals, replace = TRUE, prob = c(0.4, 0.6))
)

generate_vl <- function(age, vaccinated) {
  peak_day <- rnorm(1, mean = 4, sd = 0.8)
  peak_vl <- rnorm(1, mean = 7.5 - 0.8 * vaccinated - 0.01 * age, sd = 0.5)
  decay_rate <- rnorm(1, mean = 0.4 + 0.15 * vaccinated, sd = 0.05)
  growth_rate <- rnorm(1, mean = 1.5, sd = 0.2)

  days <- 1:n_timepoints
  vl <- peak_vl - growth_rate * pmax(0, peak_day - days) -
    decay_rate * pmax(0, days - peak_day)
  vl <- pmax(0, vl + rnorm(n_timepoints, 0, 0.3))
  vl
}

vl_long <- pop |>
  rowwise() |>
  mutate(vl = list(generate_vl(age, vaccinated))) |>
  unnest_longer(vl, indices_to = "day")

cat("Long format:", nrow(vl_long), "rows x", ncol(vl_long), "cols\n")

vl_wide <- vl_long |>
  ungroup() |>
  pivot_wider(names_from = day, values_from = vl, names_prefix = "vl_day")

cat("Wide format:", nrow(vl_wide), "rows x", ncol(vl_wide), "cols\n\n")

# --- Test 1: Wide format ---
cat("=== Test 1: Wide format (one row per individual) ===\n")
wide_df <- as.data.frame(vl_wide |> dplyr::select(-id))
wide_df$sex <- factor(wide_df$sex)
wide_df$vaccinated <- factor(wide_df$vaccinated)

syn_wide <- syn(wide_df, method = "cart", seed = 2024)

cat("\nReal vs synthetic correlation between day 1 and day 7 VL:\n")
cat("  Real:", round(cor(wide_df$vl_day1, wide_df$vl_day7), 3), "\n")
cat("  Synthetic:", round(cor(syn_wide$syn$vl_day1, syn_wide$syn$vl_day7), 3), "\n")

cat("\nReal vs synthetic correlation between day 3 and day 5 VL:\n")
cat("  Real:", round(cor(wide_df$vl_day3, wide_df$vl_day5), 3), "\n")
cat("  Synthetic:", round(cor(syn_wide$syn$vl_day3, syn_wide$syn$vl_day5), 3), "\n")

real_cors <- sapply(2:n_timepoints, function(d) {
  cor(wide_df[[paste0("vl_day", d - 1)]], wide_df[[paste0("vl_day", d)]])
})
syn_cors <- sapply(2:n_timepoints, function(d) {
  cor(syn_wide$syn[[paste0("vl_day", d - 1)]], syn_wide$syn[[paste0("vl_day", d)]])
})

cat("\nLag-1 autocorrelation (consecutive days):\n")
cat("  Real  mean:", round(mean(real_cors), 3), " range:",
    round(min(real_cors), 3), "-", round(max(real_cors), 3), "\n")
cat("  Synth mean:", round(mean(syn_cors), 3), " range:",
    round(min(syn_cors), 3), "-", round(max(syn_cors), 3), "\n")

cat("\nMean peak VL (day 4) by vaccination status:\n")
cat("  Real  - unvacc:", round(mean(wide_df$vl_day4[wide_df$vaccinated == "0"]), 2),
    " vacc:", round(mean(wide_df$vl_day4[wide_df$vaccinated == "1"]), 2), "\n")
cat("  Synth - unvacc:", round(mean(syn_wide$syn$vl_day4[syn_wide$syn$vaccinated == "0"]), 2),
    " vacc:", round(mean(syn_wide$syn$vl_day4[syn_wide$syn$vaccinated == "1"]), 2), "\n")

ug_wide <- utility.gen(syn_wide, wide_df, print.flag = FALSE)
cat("\nUtility (wide): S_pMSE =", round(ug_wide$S_pMSE, 2),
    " SPECKS =", round(ug_wide$SPECKS, 3), "\n")

# --- Test 2: Long format (ignoring individual structure) ---
cat("\n=== Test 2: Long format (ignoring individual structure) ===\n")
long_df <- as.data.frame(vl_long |> dplyr::select(-id))
long_df$sex <- factor(long_df$sex)
long_df$vaccinated <- factor(long_df$vaccinated)

syn_long <- syn(long_df, method = "cart", seed = 2024)

ug_long <- utility.gen(syn_long, long_df, print.flag = FALSE)
cat("Utility (long): S_pMSE =", round(ug_long$S_pMSE, 2),
    " SPECKS =", round(ug_long$SPECKS, 3), "\n")

# --- Test 3: Baseline cross-sectional only ---
cat("\n=== Test 3: Baseline cross-sectional only ===\n")
baseline_df <- as.data.frame(pop |> dplyr::select(-id))
baseline_df$sex <- factor(baseline_df$sex)
baseline_df$vaccinated <- factor(baseline_df$vaccinated)

syn_baseline <- syn(baseline_df, method = "cart", seed = 2024)
ug_baseline <- utility.gen(syn_baseline, baseline_df, print.flag = FALSE)
cat("Utility (baseline): S_pMSE =", round(ug_baseline$S_pMSE, 2),
    " SPECKS =", round(ug_baseline$SPECKS, 3), "\n")

# --- Plot: individual trajectories ---
set.seed(1)
sample_ids_real <- sample(unique(vl_long$id), 20)
sample_real <- vl_long |>
  filter(id %in% sample_ids_real) |>
  mutate(source = "Real",
         sex = as.character(sex),
         vaccinated = as.numeric(vaccinated))

sample_syn_wide <- syn_wide$syn
sample_syn_wide$id <- 1:nrow(sample_syn_wide)
sample_syn_long <- sample_syn_wide |>
  mutate(across(starts_with("vl_day"), ~ as.numeric(as.character(.))),
         sex = as.character(sex),
         vaccinated = as.numeric(as.character(vaccinated))) |>
  pivot_longer(starts_with("vl_day"), names_to = "day", values_to = "vl",
               names_prefix = "vl_day") |>
  mutate(day = as.integer(day))

set.seed(1)
sample_ids_syn <- sample(unique(sample_syn_long$id), 20)
sample_syn <- sample_syn_long |>
  filter(id %in% sample_ids_syn) |>
  mutate(source = "Synthetic")

plot_data <- bind_rows(sample_real, sample_syn)

p <- ggplot(plot_data, aes(x = day, y = vl, group = interaction(id, source),
                           colour = source)) +
  geom_line(alpha = 0.5) +
  scale_colour_manual(values = c(Real = "#6680c0", Synthetic = "#fc8d62")) +
  labs(x = "Day", y = "Viral load (log10 copies/mL)",
       title = "Individual viral load trajectories: real vs synthetic (wide-format CART)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"))

ggsave("/Users/billyquilty/Documents/Work/synthpop/outputs/vl_trajectories_comparison.png", p,
       width = 10, height = 6, dpi = 150)
cat("\nPlot saved to outputs/vl_trajectories_comparison.png\n")

# Mean trajectory by vaccination status
mean_real <- vl_long |>
  mutate(vaccinated = factor(vaccinated)) |>
  group_by(day, vaccinated) |>
  summarise(mean_vl = mean(vl), .groups = "drop") |>
  mutate(source = "Real")

mean_syn <- sample_syn_long |>
  group_by(day, vaccinated) |>
  summarise(mean_vl = mean(vl), .groups = "drop") |>
  mutate(source = "Synthetic", vaccinated = factor(vaccinated))

mean_data <- bind_rows(mean_real, mean_syn)

p2 <- ggplot(mean_data, aes(x = day, y = mean_vl, colour = source,
                             linetype = vaccinated)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(Real = "#6680c0", Synthetic = "#fc8d62")) +
  labs(x = "Day", y = "Mean viral load (log10 copies/mL)",
       title = "Mean viral load trajectory by vaccination status",
       linetype = "Vaccinated") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"))

ggsave("/Users/billyquilty/Documents/Work/synthpop/outputs/vl_mean_trajectories.png", p2,
       width = 10, height = 6, dpi = 150)
cat("Plot saved to outputs/vl_mean_trajectories.png\n")
