library(dplyr)      # data manipulation and dataframe processing
library(ggplot2)    # plotting
library(patchwork)  # combine plots into one figure

### set paths ###

# directory containing pairwise FST files
fst_dir <- file.path("bash_results", "fst_results")

# directory containing permutation CSV files
permutation_dir <- file.path("bash_results", "permutations")

# directory for R outputs
output_dir <- "r_results"

dir.create(output_dir, showWarnings = FALSE)

### read pairwise FST results ###

fst_files <- list.files(
  path = fst_dir,
  pattern = "\\.windowed\\.weir\\.fst$",
  full.names = TRUE
)

# read each FST file and obtain population names from the filename
fst_list <- lapply(fst_files, function(file) {
  fst_data <- read.table(
    file,
    header = TRUE
  )
  # remove VCFtools suffix from filename
  pair_name <- sub(
    "\\.windowed\\.weir\\.fst$",
    "",
    basename(file)
  )
  # split filename into population names
  populations <- strsplit(pair_name, "\\.")[[1]]
  fst_data %>%
    mutate(
      Pop1 = populations[1],
      Pop2 = populations[2]
    )
})

# combine all pairwise FST results 
obs_fst_combined <- bind_rows(fst_list) %>%
  select(
    Pop1,
    Pop2,
    WEIGHTED_FST,
    MEAN_FST
  )

# calculate genome-wide average FST
obs_summary_fst <- obs_fst_combined %>%
  group_by(Pop1, Pop2) %>%
  summarise(
    Weighted_FST = mean(WEIGHTED_FST, na.rm = TRUE),
    Mean_FST = mean(MEAN_FST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pair = paste(Pop1, Pop2, sep = ".")
  )

print(obs_summary_fst)

# Save summarized observed FST values
write.csv(
  obs_summary_fst,
  file.path(output_dir, "pairwise_fst_summary.csv"),
  row.names = FALSE
)

### read permutation results ###

permutation_files <- list.files(
  path = permutation_dir,
  pattern = "\\.average_fst_permutation\\.csv$",
  full.names = TRUE
)

fst_permutation <- lapply(permutation_files, function(file) {
  data <- read.csv(file)
  # extract population-pair name from filename
  data$pair <- sub(
    "\\.average_fst_permutation\\.csv$",
    "",
    basename(file)
  )
  data
}) %>%
  bind_rows()

### plot observed and permuted weighted FST ###

plot_weighted <- ggplot(
  fst_permutation,
  aes(x = average_weighted_fst)
) +
  geom_histogram(
    binwidth = 0.001,
    fill = "skyblue",
    color = "black"
  ) +
  geom_vline(
    data = obs_summary_fst,
    aes(xintercept = Weighted_FST),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_wrap(
    ~pair,
    scales = "free_y"
  ) +
  labs(
    x = expression("Permuted Weighted " * italic(F)[ST]),
    y = "Count"
  ) +
  theme_minimal()

### plot observed and permuted mean FST ###

plot_mean <- ggplot(
  fst_permutation,
  aes(x = average_mean_fst)
) +
  geom_histogram(
    binwidth = 0.001,
    fill = "lightgreen",
    color = "black"
  ) +
  geom_vline(
    data = obs_summary_fst,
    aes(xintercept = Mean_FST),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_wrap(
    ~pair,
    scales = "free_y"
  ) +
  labs(
    x = expression("Permuted Mean " * italic(F)[ST]),
    y = "Count"
  ) +
  theme_minimal()

### combine permutation plots ###

combined_plot <- plot_mean / plot_weighted +
  plot_annotation(tag_levels = "a")

combined_plot

# save plots
ggsave(
  file.path(output_dir, "fst_permutations.png"),
  plot = combined_plot,
  width = 8.5,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(output_dir, "fst_weighted_permutations.png"),
  plot = plot_weighted,
  width = 8.5,
  height = 10,
  dpi = 300
)

### calculate empirical P-values ###

# weighted FST
p_weighted <- fst_permutation %>%
  left_join(
    obs_summary_fst,
    by = "pair"
  ) %>%
  group_by(pair) %>%
  summarise(
    p_weighted = round(
      (sum(average_weighted_fst >= Weighted_FST, na.rm = TRUE) + 1) /
        (sum(!is.na(average_weighted_fst)) + 1),
      4
    ),
    .groups = "drop"
  )

# mean FST
p_mean <- fst_permutation %>%
  left_join(
    obs_summary_fst,
    by = "pair"
  ) %>%
  group_by(pair) %>%
  summarise(
    p_mean = round(
      (sum(average_mean_fst >= Mean_FST, na.rm = TRUE) + 1) /
        (sum(!is.na(average_mean_fst)) + 1),
      4
    ),
    .groups = "drop"
  )

# combine and save P-values

fst_pvalues <- left_join(
  p_weighted,
  p_mean,
  by = "pair"
)

print(fst_pvalues)

write.csv(
  fst_pvalues,
  file.path(output_dir, "fst_permutation_pvalues.csv"),
  row.names = FALSE
)
