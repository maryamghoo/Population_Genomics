library(dplyr)        # data manipulation and dataframe processing
library(ggplot2)      # plotting
library(patchwork)    # combine plots into one figure
library(multcompView) # generate significance letters (a, b, c, etc.) from pairwise tests

# set paths
input_dir  <- "r_data"
output_dir <- "r_results"

# creates r_results folder if it does not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# minimum number of SNPs required per window for pi analysis
min_variants <- 10

# genetic clusters in desired plotting order
clusters <- c(
  "Merritt Island",
  "SL",
  "CBC",
  "SBC",
  "Archbold",
  "Ocala"
)

# file prefixes
file_names <- c(
  "Merritt Island" = "merritt",
  "SL"             = "sl",
  "CBC"            = "cbc",
  "SBC"            = "sbc",
  "Archbold"       = "archbold",
  "Ocala"          = "ocala"
)

# colors
cluster_colors <- c(
  "Merritt Island" = "#d62728",
  "SL"             = "royalblue1",
  "CBC"            = "navy",
  "SBC"            = "chartreuse3",
  "Archbold"       = "mediumorchid4",
  "Ocala"          = "darkorange"
)


########## NUCLEOTIDE DIVERSITY (Pi) ##########

# read and combine pi files
pi_data <- bind_rows(
  lapply(clusters, function(cluster) {
    file <- file.path(
      input_dir,
      paste0(file_names[cluster], ".windowed.pi")
    )
    dat <- read.table(file, header = TRUE)
    dat$cluster <- cluster
    dat
  })
)

pi_data$cluster <- factor(pi_data$cluster, levels = clusters)

# filter windows
pi_filtered <- pi_data %>%
  filter(
    N_VARIANTS >= min_variants,
    !is.na(PI)
  )

# summary statistics
pi_summary <- pi_filtered %>%
  group_by(cluster) %>%
  summarise(
    mean_PI   = mean(PI, na.rm = TRUE),
    median_PI = median(PI, na.rm = TRUE),
    sd_PI     = sd(PI, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

print(pi_summary)

write.table(
  pi_summary,
  file.path(output_dir, "pi_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# overall significance test
pi_kruskal <- kruskal.test(PI ~ cluster, data = pi_filtered)

print(pi_kruskal)

capture.output(
  pi_kruskal,
  file = file.path(output_dir, "pi_kruskal_wallis.txt")
)

# pairwise comparisons
pi_pairwise <- pairwise.wilcox.test(
  pi_filtered$PI,
  pi_filtered$cluster,
  p.adjust.method = "bonferroni",
  exact = FALSE
)

print(pi_pairwise)

write.table(
  pi_pairwise$p.value,
  file.path(output_dir, "pi_pairwise_wilcoxon.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# generate significance letters automatically
get_significance_letters <- function(pairwise_result, groups) {
  pmat <- pairwise_result$p.value
  
  # create full symmetric matrix
  full_mat <- matrix(
    1,
    nrow = length(groups),
    ncol = length(groups),
    dimnames = list(groups, groups)
  )
  
  for (i in rownames(pmat)) {
    for (j in colnames(pmat)) {
      
      if (!is.na(pmat[i, j])) {
        full_mat[i, j] <- pmat[i, j]
        full_mat[j, i] <- pmat[i, j]
      }
    }
  }
  
  diag(full_mat) <- 1
  
  letters <- multcompLetters(
    full_mat,
    threshold = 0.05
  )$Letters
  
  data.frame(
    cluster = names(letters),
    label = letters,
    stringsAsFactors = FALSE
  )
}

pi_letters <- get_significance_letters(
  pi_pairwise,
  clusters
)

print(pi_letters)

write.table(
  pi_letters,
  file.path(output_dir, "pi_significance_letters.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

########## HETEROZYGOSITY ##########

# read and combine heterozygosity files
het_data <- bind_rows(
  lapply(clusters, function(cluster) {
    file <- file.path(
      input_dir,
      paste0(file_names[cluster], ".het")
    )
    dat <- read.table(file, header = TRUE)
    dat$cluster <- cluster
    dat
  })
)

het_data$cluster <- factor(
  het_data$cluster,
  levels = clusters
)

# calculate observed and expected heterozygosity
het_data <- het_data %>%
  mutate(
    # Observed heterozygosity
    Ho = (N_SITES - O.HOM.) / N_SITES,
    
    # Expected heterozygosity
    He = (N_SITES - E.HOM.) / N_SITES
  )

# calculate sample size for each genetic cluster
sample_sizes <- het_data %>%
  group_by(cluster) %>%
  summarise(
    N = n(),
    .groups = "drop"
  )

print(sample_sizes)

# correct expected heterozygosity for sampling bias
# He_unbiased = (2N / (2N - 1)) * He
het_data <- het_data %>%
  left_join(sample_sizes, by = "cluster") %>%
  mutate(
    He_unbiased = (2 * N) / (2 * N - 1) * He
  )

# summary statistics
het_summary <- het_data %>%
  group_by(cluster) %>%
  summarise(
    sample_size = first(N),
    
    mean_Ho = mean(Ho, na.rm = TRUE),
    sd_Ho = sd(Ho, na.rm = TRUE),
    
    mean_He = mean(He, na.rm = TRUE),
    sd_He = sd(He, na.rm = TRUE),
    
    mean_He_unbiased = mean(He_unbiased, na.rm = TRUE),
    sd_He_unbiased = sd(He_unbiased, na.rm = TRUE),
    
    .groups = "drop"
  )

print(het_summary)

# save summary
write.table(
  het_summary,
  file.path(output_dir, "heterozygosity_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# save individual-level heterozygosity
write.table(
  het_data,
  file.path(output_dir, "heterozygosity_individuals.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# significance test for observed heterozygosity
het_kruskal <- kruskal.test(
  Ho ~ cluster,
  data = het_data
)

print(het_kruskal)

capture.output(
  het_kruskal,
  file = file.path(output_dir, "heterozygosity_kruskal_wallis.txt")
)

# pairwise comparisons
het_pairwise <- pairwise.wilcox.test(
  het_data$Ho,
  het_data$cluster,
  p.adjust.method = "bonferroni",
  exact = FALSE
)

print(het_pairwise)

write.table(
  het_pairwise$p.value,
  file.path(output_dir, "heterozygosity_pairwise_wilcoxon.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# generate significance letters automatically
het_letters <- get_significance_letters(
  het_pairwise,
  clusters
)

print(het_letters)

write.table(
  het_letters,
  file.path(output_dir, "heterozygosity_significance_letters.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

########## PREPARE SIGNIFICANCE LETTER POSITIONS ##########

# put letters slightly above the highest observed value
pi_range <- range(pi_filtered$PI, na.rm = TRUE)

pi_letter_y <- pi_range[2] +
  0.05 * diff(pi_range)

pi_letters$y <- pi_letter_y


het_range <- range(het_data$Ho, na.rm = TRUE)

het_letter_y <- het_range[2] +
  0.05 * diff(het_range)

het_letters$y <- het_letter_y


# make sure order is consistent
pi_letters$cluster <- factor(
  pi_letters$cluster,
  levels = clusters
)

het_letters$cluster <- factor(
  het_letters$cluster,
  levels = clusters
)

########## PLOTS ##########

base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0
    ),
    
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    
    legend.position = "none"
  )

# plot a: nucleotide diversity
p_pi <- ggplot(
  pi_filtered,
  aes(x = cluster, y = PI, fill = cluster)
) +
  
  geom_boxplot() +
  
  geom_text(
    data = pi_letters,
    aes(
      x = cluster,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  
  labs(
    title = "(a) Nucleotide Diversity Across Genetic Clusters",
    x = "",
    y = expression(
      "Nucleotide Diversity (" * italic(pi) * ")"
    )
  ) +
  
  base_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# plot b: Observed heterozygosity

p_het <- ggplot(
  het_data,
  aes(x = cluster, y = Ho, fill = cluster)
) +
  
  geom_boxplot() +
  
  geom_text(
    data = het_letters,
    aes(
      x = cluster,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  
  labs(
    title = "(b) Observed Heterozygosity Across Genetic Clusters",
    x = "Genetic Cluster",
    y = expression(
      "Observed Heterozygosity (" * italic(H)[o] * ")"
    )
  ) +
  
  base_theme

# combine plots

combined_plot <- p_pi / p_het
combined_plot

# save
ggsave(
  file.path(output_dir, "pi_heterozygosity.pdf"),
  combined_plot,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  file.path(output_dir, "pi_heterozygosity.png"),
  combined_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(output_dir, "pi_heterozygosity.svg"),
  combined_plot,
  width = 10,
  height = 8,
  units = "in"
)
