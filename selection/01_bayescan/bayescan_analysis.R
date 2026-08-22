library(ggplot2)  # plotting
library(dplyr)    # data manipulation
library(qqman)    # Manhattan plot

# set paths
data_dir <- "r_data"
bash_dir <- "bash_results"
output_dir <- "r_results"

# creates r_results folder if it does not exist
dir.create(output_dir, showWarnings = FALSE)

##### load data #####

# load prepared BayeScan results
df <- read.table(
  file.path(bash_dir, "final_bayescan_fst.txt"),
  header = TRUE,
  sep = "",
  strip.white = TRUE,
  fill = TRUE
)

# load chromosome label mapping
chr_map <- read.table(
  file.path(data_dir, "chromosomes.txt"),
  header = FALSE,
  fill = TRUE
)

colnames(chr_map) <- c("Chr_Label", "Chr_Num", "CHR")

# merge chromosome information with BayeScan results
df <- left_join(df, chr_map, by = "CHR")

# keep only SNPs mapped to chromosomes
# removes scaffolds without chromosome mapping
df_chr <- df %>%
  filter(!is.na(Chr_Num))

# number of SNPs per chromosome
snp_counts <- df_chr %>%
  count(Chr_Num, name = "Num_SNPs") %>%
  arrange(as.numeric(gsub("[^0-9]", "", Chr_Num)))

snp_counts


##### FST distribution #####

# distribution of FST for all SNPs
fst_distribution <- ggplot(df, aes(x = fst)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_density(
    color = "darkblue",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = mean(df$fst, na.rm = TRUE),
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of FST Values",
    x = expression(italic(F)[ST]),
    y = "Density"
  ) +
  theme_minimal()

fst_distribution

# save plot
ggsave(
  file.path(output_dir, "fst_distribution.png"),
  plot = fst_distribution,
  width = 8,
  height = 6,
  dpi = 300
)


# distribution of FST for chromosome-mapped SNPs
fst_chr_distribution <- ggplot(df_chr, aes(x = fst)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_density(
    color = "darkblue",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = mean(df_chr$fst, na.rm = TRUE),
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of Chromosome FST Values",
    x = expression(italic(F)[ST]),
    y = "Density"
  ) +
  theme_minimal()

fst_chr_distribution

# save plot
ggsave(
  file.path(output_dir, "fst_chr_distribution.png"),
  plot = fst_chr_distribution,
  width = 8,
  height = 6,
  dpi = 300
)


##### identify outliers using FDR #####

# FDR < 0.05
fdr05_chr_outliers <- df_chr %>%
  filter(qval < 0.05) %>%
  select(ID, CHR, Chr_Label, POS, REF, ALT)

fdr05_outliers <- df %>%
  filter(qval < 0.05) %>%
  select(ID, CHR, Chr_Label, POS, REF, ALT)

# number of outliers
nrow(fdr05_chr_outliers)
nrow(fdr05_outliers)

# save FDR outliers
write.table(
  fdr05_chr_outliers,
  file.path(output_dir, "fdr05_chr_outliers.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

write.table(
  fdr05_outliers,
  file.path(output_dir, "fdr05_outliers.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)


##### Manhattan plot #####

# prepare dataframe for qqman
manhattan_df <- data.frame(
  CHR = as.numeric(df_chr$Chr_Num),
  BP = df_chr$POS,
  P = df_chr$fst,
  SNP = df_chr$ID
)

# function for Manhattan plot
plot_manhattan <- function() {
  
  manhattan(
    manhattan_df,
    highlight = fdr05_chr_outliers$ID,
    highlight.col = "green",
    chr = "CHR",
    bp = "BP",
    p = "P",
    logp = FALSE,
    ylim = c(0, 0.6),
    col = c("blue4", "orange3"),
    ylab = expression(italic(F)[ST]),
    xlab = "Chromosome",
    cex.axis = 1.5,
    cex.lab = 1.8
  )
}


##### save Manhattan plot #####

# PNG
png(
  file.path(output_dir, "bayescan_manhattan.png"),
  width = 14,
  height = 7,
  units = "in",
  res = 300
)

par(mar = c(5.1, 7.5, 4.1, 2.1))
plot_manhattan()
dev.off()


# PDF
pdf(
  file.path(output_dir, "bayescan_manhattan.pdf"),
  width = 12,
  height = 6
)

par(mar = c(5.1, 7.5, 4.1, 2.1))
plot_manhattan()
dev.off()
