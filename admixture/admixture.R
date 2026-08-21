library(ggplot2)  # plotting
library(dplyr)    # data manipulation and dataframe processing

# set paths
data_dir <- "r_data"
output_dir <- "r_results"

# creates r_results folder if it does not exist
dir.create(output_dir, showWarnings = FALSE)

# load CV error data
kvalue <- read.csv(file.path(data_dir, "kvalue.csv"))

# CV error plot
cv_error <- ggplot(kvalue, aes(K, CVerror)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_discrete(limits = kvalue$K) +
  labs(
    x = expression(italic(k)),
    y = "CV error"
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# save CV error plot
ggsave(
  filename = file.path(output_dir, "cv_error.png"),
  plot = cv_error,
  width = 8,
  height = 6,
  dpi = 300
)

#! from cv error plot, we found k = 5 has the least error

# change K value as needed
k_value <- 5

input_data <- read.csv(
  file.path(data_dir, paste0("adinput.", k_value, ".Q.csv")),
  header = TRUE,
  stringsAsFactors = FALSE
)

##### prepare admixture dataframe #####

# Extract ancestry columns based on K value
# For K = 5, this means columns 4 to 8
popMatrix <- as.matrix(input_data[, 4:(3 + k_value)])

# Number of individuals
N <- nrow(input_data)

# Create empty dataframe with K rows per individual
admix_df <- data.frame(matrix(0, ncol = 5, nrow = k_value * N))

names(admix_df)[1] <- "Region"
names(admix_df)[2] <- "Location"
names(admix_df)[3] <- "SampleID"
names(admix_df)[4] <- "POPGROUP"
names(admix_df)[5] <- "Fraction"

# Extract ancestry proportion columns
new_input <- input_data[, 4:(3 + k_value)]

# Fill Fraction column
# This turns K ancestry columns per sample into K rows per sample
admix_df <- admix_df %>%
  mutate(Fraction = matrix(t(new_input), ncol = 1))

# Fill metadata columns by repeating each individual K times
admix_df[, 1] <- rep(input_data[, 1], each = k_value)
admix_df[, 2] <- rep(input_data[, 2], each = k_value)
admix_df[, 3] <- rep(input_data[, 3], each = k_value)

# Assign ancestry group IDs: 1, 2, ..., K for each individual
admix_df[, 4] <- rep(1:k_value, nrow(input_data))

# save processed dataframe
write.csv(
  admix_df,
  file.path(output_dir, paste0("admix_k", k_value, ".csv")),
  row.names = FALSE
)

##### plotting #####

# define plotting order
location_order <- c(
  "Shiloh", "CNS", "Happy Creek", "Schwartz",  "RR", "TEL4",
  "SL", "Fox", "BL", "DCE", "Viera", "Cruickshank", "Malabar", 
  "Valkaria", "Micco", "Archbold", "Ocala")

admix_df$Location <- factor(
  admix_df$Location,
  levels = location_order
)

# define colors
region_colors <- c(
  "#d62728",
  "mediumorchid4",
  "chartreuse3",
  "royalblue1",
  "navy"
)

# admixture barplot
admix_plot <- ggplot(
  admix_df,
  aes(
    factor(SampleID),
    Fraction,
    fill = factor(POPGROUP)
  )
) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(
    ~Location,
    switch = "x",
    scales = "free",
    space = "free"
  ) +
  theme_minimal() +
  labs(
    x = "Individuals",
    title = expression(italic(k) ~ "= 5"),
    y = "Ancestry"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(
    values = region_colors,
    guide = FALSE
  ) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 13),
    panel.grid = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 15,
      b = 25,
      l = 15
    )
  )

admix_plot

##### save outputs #####

# PDF
ggsave(file.path(output_dir, paste0("admix_k", k_value, ".pdf")),
  plot = admix_plot,
  width = 18,
  height = 8,
  dpi = 600
)

# PNG
ggsave(
  file.path(output_dir, paste0("admix_k", k_value, ".png")),
  plot = admix_plot,
  width = 18,
  height = 8,
  dpi = 600
)

# SVG
ggsave(
  file.path(output_dir, paste0("admix_k", k_value, ".svg")),
  plot = admix_plot,
  width = 18,
  height = 8,
  dpi = 600
)
