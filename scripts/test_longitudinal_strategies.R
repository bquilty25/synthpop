devtools::load_all("/Users/billyquilty/Documents/Work/synthpop")
library(tidyverse)

set.seed(42)

# --- Simulate longitudinal data similar to a real cohort ---
# Scale up to make the problem realistic: many individuals, many timepoints
n_individuals <- 500
n_timepoints <- 30
n_clusters <- 10

pop <- tibble(
  id = 1:n_individuals,
  cluster = sample(1:n_clusters, n_individuals, replace = TRUE),
  age = round(rnorm(n_individuals, mean = 45, sd = 15)),
  sex = sample(c("M", "F"), n_individuals, replace = TRUE),
  vaccinated = sample(c(0, 1), n_individuals, replace = TRUE, prob = c(0.4, 0.6))
)

generate_vl <- function(age, vaccinated) {
  peak_day <- rnorm(1, mean = 5, sd = 1.2)
  peak_vl <- rnorm(1, mean = 7.5 - 0.8 * vaccinated - 0.01 * age, sd = 0.5)
  decay_rate <- rnorm(1, mean = 0.3 + 0.1 * vaccinated, sd = 0.05)
  growth_rate <- rnorm(1, mean = 1.2, sd = 0.2)
  days <- 1:n_timepoints
  vl <- peak_vl - growth_rate * pmax(0, peak_day - days) -
    decay_rate * pmax(0, days - peak_day)
  pmax(0, vl + rnorm(n_timepoints, 0, 0.3))
}

vl_long <- pop |>
  rowwise() |>
  mutate(vl = list(generate_vl(age, vaccinated))) |>
  unnest_longer(vl, indices_to = "day") |>
  ungroup()

cat("Data:", n_individuals, "individuals x", n_timepoints, "timepoints =",
    nrow(vl_long), "rows\n\n")

# Helper: compute lag-1 autocorrelation within individuals
lag1_autocor <- function(long_df, id_col = "id", day_col = "day",
                         value_col = "vl") {
  long_df |>
    arrange(.data[[id_col]], .data[[day_col]]) |>
    group_by(.data[[id_col]]) |>
    summarise(r = cor(.data[[value_col]][-n()],
                      .data[[value_col]][-1],
                      use = "complete.obs"),
              .groups = "drop") |>
    pull(r) |>
    mean(na.rm = TRUE)
}

# Helper: vaccination effect at peak
vax_effect <- function(wide_df, day = 5) {
  col <- paste0("vl_day", day)
  if (!col %in% names(wide_df)) return(c(unvacc = NA, vacc = NA))
  vacc_col <- if (is.factor(wide_df$vaccinated)) {
    as.numeric(as.character(wide_df$vaccinated))
  } else {
    wide_df$vaccinated
  }
  c(unvacc = mean(wide_df[[col]][vacc_col == 0], na.rm = TRUE),
    vacc = mean(wide_df[[col]][vacc_col == 1], na.rm = TRUE))
}

real_lag1 <- lag1_autocor(vl_long)
cat("Real data lag-1 autocorrelation:", round(real_lag1, 3), "\n")

real_wide <- vl_long |>
  pivot_wider(names_from = day, values_from = vl, names_prefix = "vl_day")

real_vax <- vax_effect(real_wide)
cat("Real vaccination effect at day 5:",
    round(real_vax["unvacc"], 2), "vs", round(real_vax["vacc"], 2), "\n\n")

results <- list()

# ============================================================
# Strategy 1: Wide format (one row per individual)
# This is the naive approach. Works if n_timepoints is modest.
# ============================================================
cat("=== Strategy 1: Wide format ===\n")
t1 <- system.time({
  wide_df <- as.data.frame(
    real_wide |> dplyr::select(-id, -cluster)
  )
  wide_df$sex <- factor(wide_df$sex)
  wide_df$vaccinated <- factor(wide_df$vaccinated)
  syn_wide <- tryCatch(
    syn(wide_df, method = "cart", seed = 2024),
    error = function(e) {
      cat("  FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
})

if (!is.null(syn_wide)) {
  syn_wide_long <- syn_wide$syn |>
    mutate(id = row_number(),
           across(starts_with("vl_day"), ~ as.numeric(as.character(.))),
           vaccinated = as.numeric(as.character(vaccinated))) |>
    pivot_longer(starts_with("vl_day"), names_to = "day", values_to = "vl",
                 names_prefix = "vl_day") |>
    mutate(day = as.integer(day))
  syn_lag1 <- lag1_autocor(syn_wide_long)
  syn_vax <- vax_effect(syn_wide$syn)
  ug <- utility.gen(syn_wide, wide_df, print.flag = FALSE)
  results$wide <- list(
    time = t1["elapsed"], lag1 = syn_lag1,
    vax = syn_vax, S_pMSE = ug$S_pMSE, SPECKS = ug$SPECKS
  )
  cat("  Time:", round(t1["elapsed"], 1), "s\n")
  cat("  Lag-1 autocorrelation:", round(syn_lag1, 3),
      "(real:", round(real_lag1, 3), ")\n")
  cat("  Vax effect:", round(syn_vax["unvacc"], 2), "vs",
      round(syn_vax["vacc"], 2), "\n")
  cat("  S_pMSE:", round(ug$S_pMSE, 2), " SPECKS:", round(ug$SPECKS, 3), "\n")
} else {
  results$wide <- list(time = t1["elapsed"], failed = TRUE)
}

# ============================================================
# Strategy 2: Long format (ignoring individual structure)
# Fast but destroys within-individual correlation.
# ============================================================
cat("\n=== Strategy 2: Long format (no individual structure) ===\n")
t2 <- system.time({
  long_df <- as.data.frame(
    vl_long |> dplyr::select(-id, -cluster)
  )
  long_df$sex <- factor(long_df$sex)
  long_df$vaccinated <- factor(long_df$vaccinated)
  syn_long <- syn(long_df, method = "cart", seed = 2024)
})

syn_long_data <- syn_long$syn |>
  mutate(id = rep(1:(n() %/% n_timepoints), each = n_timepoints)[1:n()])
syn_lag1_long <- lag1_autocor(syn_long_data)
ug_long <- utility.gen(syn_long, long_df, print.flag = FALSE)

results$long <- list(
  time = t2["elapsed"], lag1 = syn_lag1_long,
  S_pMSE = ug_long$S_pMSE, SPECKS = ug_long$SPECKS
)
cat("  Time:", round(t2["elapsed"], 1), "s\n")
cat("  Lag-1 autocorrelation:", round(syn_lag1_long, 3),
    "(real:", round(real_lag1, 3), ")\n")
cat("  S_pMSE:", round(ug_long$S_pMSE, 2),
    " SPECKS:", round(ug_long$SPECKS, 3), "\n")

# ============================================================
# Strategy 3: Two-stage synthesis
# Stage 1: synthesise baseline characteristics
# Stage 2: synthesise trajectory summary parameters per individual,
#           then reconstruct timepoints from those parameters
# This preserves individual-level structure without wide explosion.
# ============================================================
cat("\n=== Strategy 3: Two-stage (baseline + trajectory summaries) ===\n")
t3 <- system.time({
  # Summarise each individual's trajectory into a few parameters
  traj_summaries <- vl_long |>
    group_by(id) |>
    summarise(
      vl_peak = max(vl),
      vl_peak_day = day[which.max(vl)],
      vl_mean = mean(vl),
      vl_auc = sum(vl),
      vl_last = vl[n()],
      .groups = "drop"
    )

  baseline_with_traj <- pop |>
    dplyr::select(-cluster) |>
    left_join(traj_summaries, by = "id") |>
    dplyr::select(-id)

  baseline_df <- as.data.frame(baseline_with_traj)
  baseline_df$sex <- factor(baseline_df$sex)
  baseline_df$vaccinated <- factor(baseline_df$vaccinated)

  syn_baseline <- syn(baseline_df, method = "cart", seed = 2024)

  # Reconstruct trajectories from synthetic summaries using a simple
  # parametric model (rise-then-decay)
  reconstruct_trajectory <- function(peak_vl, peak_day, n_days) {
    peak_day <- max(1, min(peak_day, n_days))
    growth <- peak_vl / peak_day
    decay <- peak_vl / (n_days - peak_day + 1)
    days <- 1:n_days
    vl <- ifelse(days <= peak_day,
                 growth * days,
                 peak_vl - decay * (days - peak_day))
    pmax(0, vl + rnorm(n_days, 0, 0.3))
  }

  syn_base <- syn_baseline$syn
  syn_base$id <- 1:nrow(syn_base)

  syn_traj_long <- syn_base |>
    rowwise() |>
    mutate(
      vl = list(reconstruct_trajectory(
        as.numeric(as.character(vl_peak)),
        as.integer(as.character(vl_peak_day)),
        n_timepoints
      ))
    ) |>
    unnest_longer(vl, indices_to = "day") |>
    ungroup()
})

syn_lag1_twostage <- lag1_autocor(syn_traj_long)

# Vaccination effect from the synthesised baseline summaries
syn_vax_ts <- c(
  unvacc = mean(
    as.numeric(as.character(
      syn_base$vl_peak[syn_base$vaccinated == levels(syn_base$vaccinated)[1]]
    )),
    na.rm = TRUE
  ),
  vacc = mean(
    as.numeric(as.character(
      syn_base$vl_peak[syn_base$vaccinated == levels(syn_base$vaccinated)[2]]
    )),
    na.rm = TRUE
  )
)

results$twostage <- list(
  time = t3["elapsed"], lag1 = syn_lag1_twostage, vax_peak = syn_vax_ts
)
cat("  Time:", round(t3["elapsed"], 1), "s\n")
cat("  Lag-1 autocorrelation:", round(syn_lag1_twostage, 3),
    "(real:", round(real_lag1, 3), ")\n")
cat("  Synthetic peak VL by vax:",
    round(syn_vax_ts["unvacc"], 2), "vs",
    round(syn_vax_ts["vacc"], 2), "\n")
cat("  Real peak VL by vax:",
    round(real_vax["unvacc"], 2), "vs",
    round(real_vax["vacc"], 2), "\n")

# ============================================================
# Strategy 4: Synthesise per cluster, then stack
# Each cluster is synthesised independently in wide format.
# Keeps the problem small per call but loses between-cluster variation.
# ============================================================
cat("\n=== Strategy 4: Per-cluster wide synthesis ===\n")
t4 <- system.time({
  cluster_syns <- list()
  for (cl in sort(unique(pop$cluster))) {
    cl_ids <- pop$id[pop$cluster == cl]
    cl_wide <- real_wide |>
      filter(id %in% cl_ids) |>
      dplyr::select(-id, -cluster) |>
      as.data.frame()
    cl_wide$sex <- factor(cl_wide$sex)
    cl_wide$vaccinated <- factor(cl_wide$vaccinated)
    cl_syn <- tryCatch(
      syn(cl_wide, method = "cart", seed = 2024 + cl),
      error = function(e) {
        cat("  Cluster", cl, "FAILED:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(cl_syn)) {
      cluster_syns[[as.character(cl)]] <- cl_syn$syn
    }
  }
  syn_cluster_wide <- bind_rows(cluster_syns, .id = "cluster")
})

syn_cluster_long <- syn_cluster_wide |>
  mutate(id = row_number(),
         across(starts_with("vl_day"), ~ as.numeric(as.character(.))),
         vaccinated = as.numeric(as.character(vaccinated))) |>
  pivot_longer(starts_with("vl_day"), names_to = "day", values_to = "vl",
               names_prefix = "vl_day") |>
  mutate(day = as.integer(day))

syn_lag1_cluster <- lag1_autocor(syn_cluster_long)

results$per_cluster <- list(
  time = t4["elapsed"], lag1 = syn_lag1_cluster,
  n_clusters_ok = length(cluster_syns)
)
cat("  Time:", round(t4["elapsed"], 1), "s\n")
cat("  Clusters synthesised:", length(cluster_syns), "/", n_clusters, "\n")
cat("  Lag-1 autocorrelation:", round(syn_lag1_cluster, 3),
    "(real:", round(real_lag1, 3), ")\n")

# ============================================================
# Strategy 5: Wide format with parametric method
# Parametric fits linear models, which are faster and more stable
# than CART for many columns.
# ============================================================
cat("\n=== Strategy 5: Wide format, parametric method ===\n")
t5 <- system.time({
  wide_df_p <- as.data.frame(
    real_wide |> dplyr::select(-id, -cluster)
  )
  wide_df_p$sex <- factor(wide_df_p$sex)
  wide_df_p$vaccinated <- factor(wide_df_p$vaccinated)
  syn_wide_p <- tryCatch(
    syn(wide_df_p, method = "parametric", seed = 2024),
    error = function(e) {
      cat("  FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
})

if (!is.null(syn_wide_p)) {
  syn_wp_long <- syn_wide_p$syn |>
    mutate(id = row_number(),
           vaccinated = as.numeric(as.character(vaccinated))) |>
    pivot_longer(starts_with("vl_day"), names_to = "day", values_to = "vl",
                 names_prefix = "vl_day") |>
    mutate(day = as.integer(day))
  syn_lag1_p <- lag1_autocor(syn_wp_long)
  syn_vax_p <- vax_effect(syn_wide_p$syn)
  ug_p <- utility.gen(syn_wide_p, wide_df_p, print.flag = FALSE)
  results$wide_param <- list(
    time = t5["elapsed"], lag1 = syn_lag1_p,
    vax = syn_vax_p, S_pMSE = ug_p$S_pMSE, SPECKS = ug_p$SPECKS
  )
  cat("  Time:", round(t5["elapsed"], 1), "s\n")
  cat("  Lag-1 autocorrelation:", round(syn_lag1_p, 3),
      "(real:", round(real_lag1, 3), ")\n")
  cat("  Vax effect:", round(syn_vax_p["unvacc"], 2), "vs",
      round(syn_vax_p["vacc"], 2), "\n")
  cat("  S_pMSE:", round(ug_p$S_pMSE, 2),
      " SPECKS:", round(ug_p$SPECKS, 3), "\n")
} else {
  results$wide_param <- list(time = t5["elapsed"], failed = TRUE)
}

# ============================================================
# Summary
# ============================================================
cat("\n\n========== SUMMARY ==========\n")
cat(sprintf("%-30s %8s %8s\n", "Strategy", "Time(s)", "Lag-1 r"))
cat(sprintf("%-30s %8s %8s\n", "Real data", "", round(real_lag1, 3)))
for (nm in names(results)) {
  r <- results[[nm]]
  lag <- if (!is.null(r$lag1)) round(r$lag1, 3) else "FAILED"
  cat(sprintf("%-30s %8.1f %8s\n", nm, r$time, lag))
}
cat("\nKey question: which strategy best preserves within-individual\n")
cat("correlation (lag-1 r) AND between-group effects (vax effect)?\n")
