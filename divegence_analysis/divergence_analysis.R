# === Required Libraries ===
library(dplyr)
library(tidyverse)
library(broom)
library(forcats)
library(car)

# === PART 1: SYNONYMOUS GLM ===

# Load data
df_raw_syn <- read.csv("path/to/synonymous_data.csv")

# Compute jumps and prepare dataframe
df_syn <- df_raw_syn %>%
  filter(patient_id != "N7") %>%
  arrange(patient_id, time) %>%
  group_by(patient_id) %>%
  mutate(jump = if_else(row_number() == 1, 0, divergence - lag(divergence))) %>%
  ungroup()

# Plot observed synonymous divergence
ggplot(df_syn, aes(x = time, y = divergence, color = patient_id)) +
  geom_point(size = 2) +
  labs(x = "Days", y = "Observed synonymous divergence", color = "Patient ID") +
  theme_minimal(base_size = 16)

# Fit GLM with binary jump threshold
rate_thresh <- 1e-5
df2 <- df_raw_syn %>%
  filter(patient_id != "N7") %>%
  arrange(patient_id, time) %>%
  group_by(patient_id) %>%
  mutate(dt = time - lag(time),
         delta_div = divergence - lag(divergence),
         rate = delta_div / dt,
         jump = replace_na(rate > rate_thresh, FALSE)) %>%
  ungroup()

fit <- glm(divergence ~ time + jump, family = gaussian(), data = df2)
coef_tab <- broom::tidy(fit, conf.int = TRUE)

mu_row <- coef_tab %>% filter(term == "time")
mu_hat <- mu_row$estimate
mu_low <- mu_row$conf.low
mu_high <- mu_row$conf.high

sci <- function(x, digits = 3) formatC(x, format = "e", digits = digits)
cat("Mutation-rate estimate:", sci(mu_hat), "95% CI:", sci(mu_low), "-", sci(mu_high), "\n")
syn_rate <- mu_hat
syn_CI <- c(mu_low, mu_high)

# Evaluate thresholds
grid <- tibble(rate_thresh = 10^seq(-6, -4, length = 9))
rate_grid <- grid %>% mutate(mu_hat = map_dbl(rate_thresh, function(th) {
  df_tmp <- df_raw_syn %>%
    filter(patient_id != "N7") %>%
    arrange(patient_id, time) %>%
    group_by(patient_id) %>%
    mutate(dt = time - lag(time),
           delta_div = divergence - lag(divergence),
           rate = delta_div / dt,
           jump = replace_na(rate > th, FALSE)) %>%
    ungroup() %>%
    filter(!is.na(dt))
  coef(glm(divergence ~ time + jump, data = df_tmp, family = gaussian()))["time"]
}))

print(rate_grid)

# Bootstrap for confidence interval
set.seed(123)
n_boot <- 1000
bootstrap_results2 <- data.frame(mutation_rate = numeric(n_boot))

for (i in 1:n_boot) {
  boot_df <- df2 %>% sample_frac(replace = TRUE)
  glm_fit <- glm(divergence ~ time + jump, data = boot_df, family = gaussian())
  bootstrap_results2$mutation_rate[i] <- coef(glm_fit)["time"]
}

bootstrap_median_syn_rate2 <- median(bootstrap_results2$mutation_rate)
bootstrap_se2 <- sd(bootstrap_results2$mutation_rate)
ci <- quantile(bootstrap_results2$mutation_rate, probs = c(0.025, 0.975))
cat("Bootstrapped Mutation-rate median estimate:", sci(bootstrap_median_syn_rate2), "95% CI:", ci, "\n")

# GLM with time limit
max_days <- 180
filtered_df <- df_syn %>% filter(time <= max_days)
model_filtered <- glm(divergence ~ time + jump, data = filtered_df, family = gaussian())
summary(model_filtered)

# === PART 2: NON-SYNONYMOUS REGRESSION ===
data_dir <- "path/to/non_syn_data/"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

all_res <- map_dfr(files, function(f) {
  patient_id <- sub("_.*$", "", tools::file_path_sans_ext(basename(f)))
  if (patient_id == "N7") return(NULL)

  df <- read_csv(f, show_col_types = FALSE) %>%
    rename(time = timepoint, divergence = final_freq, patient_id = patient) %>%
    arrange(time) %>%
    mutate(jump = divergence - lag(divergence))

  if (nrow(df) <= 2) return(NULL)

  mod <- glm(divergence ~ time, data = df, family = gaussian())
  sm <- summary(mod)
  est <- coef(sm)["time", "Estimate"]
  se <- coef(sm)["time", "Std. Error"]
  pval <- coef(sm)["time", "Pr(>|t|)"]
  se_diff <- sqrt(se^2 + bootstrap_se2^2)
  z <- (est - syn_rate) / se_diff
  p_gt_ds <- pnorm(z, lower.tail = FALSE)

  tibble(
    patient = patient_id,
    rate = est,
    p.value = pval,
    p_gt_ds = p_gt_ds,
    lower_CI = est - 1.96 * se,
    upper_CI = est + 1.96 * se
  )
}) %>%
  mutate(
    reg_sig = p.value < 0.05,
    star = case_when(
      p_gt_ds < 0.005 ~ "***",
      p_gt_ds < 0.01  ~ "**",
      p_gt_ds < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    dn_ds = rate / syn_rate
  )

# Plot non-syn divergence vs syn rate
breaks <- pretty(range(c(all_res$lower_CI, all_res$upper_CI, syn_CI)), n = 5)

ggplot(all_res, aes(x = rate, y = fct_reorder(patient, rate))) +
  annotate("rect", xmin = syn_CI[1], xmax = syn_CI[2], ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.15) +
  geom_vline(xintercept = syn_rate, color = "steelblue", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +
  geom_point(aes(shape = reg_sig), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
  geom_text(aes(label = star), nudge_y = 0.25, size = 6) +
  scale_x_continuous(
    name = "divN (mut/site/day)",
    breaks = breaks,
    labels = function(x) sprintf("%.1e", x),
    sec.axis = sec_axis(~ . / syn_rate, name = "divN/divS", breaks = breaks / syn_rate, labels = function(x) sprintf("%.1f", x))
  ) +
  labs(y = "Patient") +
  theme_minimal(base_size = 16)

# === PART 3: Detect large divergence jumps ===
all_patients_data <- map_dfr(files, function(f) {
  patient_id <- sub("_.*$", "", tools::file_path_sans_ext(basename(f)))
  if (patient_id == "N7") return(NULL)
  df <- read_csv(f, show_col_types = FALSE) %>%
    rename(time = timepoint, divergence = final_freq, patient_id = patient) %>%
    arrange(time) %>%
    mutate(jump = divergence - lag(divergence), patient_id = patient_id)
  df
}) %>%
  arrange(patient_id, time) %>%
  group_by(patient_id) %>%
  mutate(delta_divergence = divergence - lag(divergence),
         delta_time = time - lag(time),
         observed_rate = delta_divergence / delta_time) %>%
  ungroup() %>%
  mutate(big_jump = observed_rate > (2 * syn_rate))

plot_df <- all_patients_data %>%
  group_by(patient_id) %>%
  mutate(max_day = max(time)) %>%
  ungroup() %>%
  mutate(patient_id = fct_reorder(patient_id, max_day, .desc = FALSE))

ggplot(plot_df, aes(x = time, y = divergence, color = big_jump)) +
  geom_point(size = 2) +
  geom_line(aes(group = patient_id), color = "gray50") +
  facet_wrap(~ patient_id, scales = "free_x") +
  scale_color_manual(values = c("black", "red"), labels = c("Normal", "Elevated")) +
  labs(x = "Days", y = "Non-synonymous divergence", color = "Rate Status") +
  theme_minimal(base_size = 14)