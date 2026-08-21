library(dplyr)        # data manipulation and dataframe processing
library(ggplot2)      # plotting
library(patchwork)    # combine plots into one figure
library(multcompView) # generate significance letters from pairwise tests

# set paths
input_dir  <- "r_data"
output_dir <- "r_results"

# creates r_results folder if it does not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# callable genome size in Mb
callable_genome_mb <- 948.29

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


########## Helper Function ##########

# generate significance letters automatically
get_significance_letters <- function(pairwise_result, groups) {
  pmat <- pairwise_result$p.value
  
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


########## FIS ##########

# read and combine vcftools --het files
het_data <- bind_rows(
  lapply(clusters, function(cluster) {
    file <- file.path(
      input_dir,
      paste0(file_names[cluster], ".het")
    )
    
    dat <- read.table(
      file,
      header = TRUE
    )
    
    dat$cluster <- cluster
    
    dat
  })
)

# rename individual ID column
het_data <- het_data %>%
  rename(IID = INDV)

# set cluster order
het_data$cluster <- factor(
  het_data$cluster,
  levels = clusters
)

# summary statistics
fis_summary <- het_data %>%
  group_by(cluster) %>%
  summarise(
    mean_F = mean(F, na.rm = TRUE),
    median_F = median(F, na.rm = TRUE),
    min_F = min(F, na.rm = TRUE),
    max_F = max(F, na.rm = TRUE),
    sd_F = sd(F, na.rm = TRUE),
    n = sum(!is.na(F)),
    .groups = "drop"
  )

print(fis_summary)

write.table(
  fis_summary,
  file.path(output_dir, "fis_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# overall significance test
fis_kruskal <- kruskal.test(
  F ~ cluster,
  data = het_data
)

print(fis_kruskal)

capture.output(
  fis_kruskal,
  file = file.path(output_dir, "fis_kruskal_wallis.txt")
)

# pairwise comparisons
fis_pairwise <- pairwise.wilcox.test(
  het_data$F,
  het_data$cluster,
  p.adjust.method = "bonferroni",
  exact = FALSE
)

print(fis_pairwise)

write.table(
  fis_pairwise$p.value,
  file.path(output_dir, "fis_pairwise_wilcoxon.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# generate significance letters
fis_letters <- get_significance_letters(
  fis_pairwise,
  clusters
)

print(fis_letters)

write.table(
  fis_letters,
  file.path(output_dir, "fis_significance_letters.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# test whether F differs from zero within each cluster
fis_zero_test <- het_data %>%
  group_by(cluster) %>%
  summarise(
    p_value = wilcox.test(F, mu = 0)$p.value,
    median_F = median(F, na.rm = TRUE),
    n = sum(!is.na(F)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "bonferroni"
    )
  )

print(fis_zero_test)

write.table(
  fis_zero_test,
  file.path(output_dir, "fis_zero_test.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


########## FROH ##########

# read individual-level PLINK ROH summary
roh_indiv <- read.csv(
  file.path(
    input_dir,
    "roh.hom.indiv.csv"
  )
)

# set cluster order
roh_indiv$cluster <- factor(
  roh_indiv$cluster,
  levels = clusters
)

# convert total ROH length from kb to Mb
# and calculate FROH
roh_indiv <- roh_indiv %>%
  mutate(
    MB = KB / 1000,
    FROH = MB / callable_genome_mb
  )

# summary statistics
froh_summary <- roh_indiv %>%
  group_by(cluster) %>%
  summarise(
    mean_FROH = mean(FROH, na.rm = TRUE),
    median_FROH = median(FROH, na.rm = TRUE),
    min_FROH = min(FROH, na.rm = TRUE),
    max_FROH = max(FROH, na.rm = TRUE),
    sd_FROH = sd(FROH, na.rm = TRUE),
    n = sum(!is.na(FROH)),
    se_FROH = sd_FROH / sqrt(n),
    .groups = "drop"
  )

print(froh_summary)

write.table(
  froh_summary,
  file.path(output_dir, "froh_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# overall significance test
froh_kruskal <- kruskal.test(
  FROH ~ cluster,
  data = roh_indiv
)

print(froh_kruskal)

capture.output(
  froh_kruskal,
  file = file.path(output_dir, "froh_kruskal_wallis.txt")
)

# pairwise comparisons
froh_pairwise <- pairwise.wilcox.test(
  roh_indiv$FROH,
  roh_indiv$cluster,
  p.adjust.method = "bonferroni",
  exact = FALSE
)

print(froh_pairwise)

write.table(
  froh_pairwise$p.value,
  file.path(output_dir, "froh_pairwise_wilcoxon.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# generate significance letters
froh_letters <- get_significance_letters(
  froh_pairwise,
  clusters
)

print(froh_letters)

write.table(
  froh_letters,
  file.path(output_dir, "froh_significance_letters.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


########## Prepare significance letters position ##########

# FIS letters
fis_range <- range(
  het_data$F,
  na.rm = TRUE
)

fis_letter_y <- fis_range[2] +
  0.05 * diff(fis_range)

fis_letters$y <- fis_letter_y

fis_letters$cluster <- factor(
  fis_letters$cluster,
  levels = clusters
)

# FROH letters
froh_range <- range(
  roh_indiv$FROH,
  na.rm = TRUE
)

froh_letter_y <- froh_range[2] +
  0.05 * diff(froh_range)

froh_letters$y <- froh_letter_y

froh_letters$cluster <- factor(
  froh_letters$cluster,
  levels = clusters
)


########## FIS and FROH plot ##########

# common plot theme
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


# plot a: FIS
p_fis <- ggplot(
  het_data,
  aes(
    x = cluster,
    y = F,
    fill = cluster
  )
) +
  
  geom_boxplot() +
  
  geom_text(
    data = fis_letters,
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
    expand = expansion(
      mult = c(0.02, 0.12)
    )
  ) +
  
  labs(
    title = expression(
      "(a) Inbreeding Coefficient (" *
        italic(F)[IS] *
        ") Across Genetic Clusters"
    ),
    x = "",
    y = expression(
      "Inbreeding Coefficient (" *
        italic(F)[IS] *
        ")"
    )
  ) +
  
  base_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# plot b: FROH
p_froh <- ggplot(
  roh_indiv,
  aes(
    x = cluster,
    y = FROH,
    fill = cluster
  )
) +
  
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.2,
    size = 1.5,
    alpha = 0.6
  ) +
  
  geom_text(
    data = froh_letters,
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
    expand = expansion(
      mult = c(0.02, 0.12)
    )
  ) +
  
  labs(
    title = expression(
      "(b) Inbreeding Coefficient (" *
        italic(F)[ROH] *
        ") Across Genetic Clusters"
    ),
    x = "Genetic Cluster",
    y = expression(
      "Inbreeding Coefficient (" *
        italic(F)[ROH] *
        ")"
    )
  ) +
  
  base_theme


# combine plots
combined_plot <- p_fis / p_froh
combined_plot

# save plots
ggsave(
  file.path(
    output_dir,
    "fis_froh.pdf"
  ),
  combined_plot,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  file.path(
    output_dir,
    "fis_froh.png"
  ),
  combined_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "fis_froh.svg"
  ),
  combined_plot,
  width = 10,
  height = 8,
  units = "in"
)


########## ROH summary statistics ##########

# read segment-level PLINK ROH data
roh_data <- read.csv(
  file.path(
    input_dir,
    "roh.hom.csv"
  )
)

# set cluster order
roh_data$cluster <- factor(
  roh_data$cluster,
  levels = clusters
)

# convert ROH length from kb to Mb
roh_data <- roh_data %>%
  mutate(
    MB = KB / 1000
  )

# summarize ROH segments by cluster
roh_summary <- roh_data %>%
  group_by(cluster) %>%
  summarise(
    total_roh_length_kb = sum(KB, na.rm = TRUE),
    total_roh_segments = n(),
    average_roh_length_kb = mean(KB, na.rm = TRUE),
    .groups = "drop"
  )

print(roh_summary)

write.table(
  roh_summary,
  file.path(output_dir, "roh_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


########## ROH length distribution ##########

# assign ROH segments to length categories
roh_length_data <- roh_data %>%
  filter(
    MB >= 1,
    MB <= 8
  ) %>%
  
  mutate(
    length_class = cut(
      MB,
      breaks = c(
        1,
        2,
        3,
        4,
        5,
        8
      ),
      labels = c(
        "1–2 Mb",
        "2–3 Mb",
        "3–4 Mb",
        "4–5 Mb",
        "5–7 Mb"
      ),
      include.lowest = TRUE
    )
  )

# calculate proportion of ROH segments in each length category
roh_length_summary <- roh_length_data %>%
  count(
    cluster,
    length_class
  ) %>%
  
  group_by(cluster) %>%
  
  mutate(
    proportion = n / sum(n)
  ) %>%
  
  ungroup()

print(roh_length_summary)

write.table(
  roh_length_summary,
  file.path(
    output_dir,
    "roh_length_distribution.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# plot ROH length distribution
roh_length_distribution <- ggplot(
  roh_length_summary,
  aes(
    x = length_class,
    y = proportion,
    fill = cluster
  )
) +
  
  geom_col(
    position = "dodge",
    color = "black"
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  labs(
    x = "ROH Length",
    y = "Proportion of ROH Segments per Genetic Cluster",
    fill = "Genetic Cluster"
  ) +
  
  theme_light() +
  
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    axis.line = element_line(
      color = "black"
    ),
    
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

roh_length_distribution

# save ROH length distribution
ggsave(
  file.path(
    output_dir,
    "roh_length_distribution.png"
  ),
  roh_length_distribution,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "roh_length_distribution.pdf"
  ),
  roh_length_distribution,
  width = 10,
  height = 6
)


########## Distribution of ROH across chromosomes ##########

# read ROH data containing chromosome information
roh_chr <- read.csv(
  file.path(
    input_dir,
    "roh.hom_chr.csv"
  )
)

# chromosome order
chromosomes <- c(
  "1", "1A", "2", "3", "4", "4A",
  "5", "6", "7", "8", "9", "10",
  "11", "12", "13", "14", "15",
  "17", "18", "19", "20", "21",
  "22", "23", "27"
)

# set factor order
roh_chr$cluster <- factor(
  roh_chr$cluster,
  levels = clusters
)

roh_chr$CHR_NUM <- factor(
  roh_chr$CHR_NUM,
  levels = chromosomes
)

# count ROH segments by chromosome and cluster
roh_chr_summary <- roh_chr %>%
  group_by(
    cluster,
    CHR_NUM
  ) %>%
  
  summarise(
    n_segments = n(),
    .groups = "drop"
  )

print(roh_chr_summary)
    
write.table(
  roh_chr_summary,
  file.path(
    output_dir,
    "roh_chromosome_distribution.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# plot
p_roh_chr <- ggplot(
  roh_chr_summary,
  aes(
    x = CHR_NUM,
    y = n_segments,
    fill = cluster
  )
) +
  
  geom_col(
    position = "stack",
    color = "black",
    alpha = 0.7
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  labs(
    title = "Distribution of ROH Segments Across the Genome",
    x = "Chromosome Number",
    y = "Number of ROH Segments",
    fill = "Genetic Cluster"
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    ),
    
    axis.title = element_text(size = 14),
    
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(size = 12)
  )

p_roh_chr

# save
ggsave(
  file.path(
    output_dir,
    "roh_seg_distribution.png"
  ),
  p_roh_chr,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "roh_seg_distribution.pdf"
  ),
  p_roh_chr,
  width = 10,
  height = 6
)

########## Prepare combined inbreeding data ##########

# keep required columns
fis_subset <- het_data %>%
  select(
    IID,
    F
  )

froh_subset <- roh_indiv %>%
  select(
    IID,
    FROH,
    cluster,
    location
  )

# merge FIS and FROH
inbreeding_data <- full_join(
  fis_subset,
  froh_subset,
  by = "IID"
)

# save combined data
write.csv(
  inbreeding_data,
  file.path(
    output_dir,
    "inbreeding_data.csv"
  ),
  row.names = FALSE
)


########## FIS and FROH correlation ##########

# overall Spearman correlation
overall_corr <- cor.test(
  inbreeding_data$F,
  inbreeding_data$FROH,
  method = "spearman",
  exact = FALSE
)

print(overall_corr)

capture.output(
  overall_corr,
  file = file.path(
    output_dir,
    "fis_froh_correlation.txt"
  )
)

# Spearman correlation within each cluster
cluster_corr <- inbreeding_data %>%
  group_by(cluster) %>%
  summarise(
    rho = cor(
      F,
      FROH,
      method = "spearman",
      use = "complete.obs"
    ),
    
    p_value = cor.test(
      F,
      FROH,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    n = sum(
      complete.cases(
        F,
        FROH
      )
    ),
    
    .groups = "drop"
  )

print(cluster_corr)

write.table(
  cluster_corr,
  file.path(
    output_dir,
    "fis_froh_correlation_by_cluster.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


########## Correlation plot ##########

# prepare labels for cluster-specific correlations
cluster_stats <- cluster_corr %>%
  mutate(
    label = paste0(
      "ρ = ",
      round(rho, 2),
      ", p = ",
      format.pval(
        p_value,
        digits = 2,
        eps = 0.001
      )
    )
  )

# overall correlation label
overall_label <- paste0(
  "Overall: ρ = ",
  round(
    unname(overall_corr$estimate),
    2
  ),
  ", p = ",
  format.pval(
    overall_corr$p.value,
    digits = 2,
    eps = 0.001
  )
)

# order labels
cluster_stats$cluster <- factor(
  cluster_stats$cluster,
  levels = clusters
)

cluster_stats <- cluster_stats %>%
  arrange(cluster)

# correlation plot
p_corr <- ggplot(
  inbreeding_data,
  aes(
    x = F,
    y = FROH,
    color = cluster
  )
) +
  
  geom_point(
    alpha = 0.6,
    size = 2
  ) +
  
  # cluster-specific trend lines
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1.2
  ) +
  
  # overall trend line
  geom_smooth(
    aes(group = 1),
    method = "lm",
    color = "black",
    se = FALSE,
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  scale_color_manual(
    values = cluster_colors
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = overall_label,
    hjust = -0.1,
    vjust = 2,
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text(
    data = cluster_stats,
    aes(
      x = -Inf,
      y = Inf,
      label = label,
      color = cluster
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = seq(
      4,
      by = 1.2,
      length.out = nrow(cluster_stats)
    ),
    size = 4,
    show.legend = FALSE
  ) +
  
  labs(
    x = expression(
      italic(F)[IS]
    ),
    y = expression(
      italic(F)[ROH]
    ),
    color = "Genetic Cluster"
  ) +
  
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )

p_corr

# save correlation plot
ggsave(
  file.path(
    output_dir,
    "fis_froh_correlation.pdf"
  ),
  p_corr,
  width = 8,
  height = 6,
  units = "in"
)

ggsave(
  file.path(
    output_dir,
    "fis_froh_correlation.png"
  ),
  p_corr,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
