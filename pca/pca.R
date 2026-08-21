# load libraries
library(ggplot2)      # plotting
library(cowplot)      # combine multiple plots into one figure
library(RColorBrewer) # color palettes for consistent colors

# set paths
data_dir <- "r_data"
output_dir <- "r_results"

# creates r_results folder if it does not exist
dir.create(output_dir, showWarnings = FALSE)

##### helper functions #####

# rotate PCA axes
rotate_pca <- function(df, angle_deg) {
  angle_rad <- angle_deg * pi / 180
  
  df$PC1_rot <- df$PC01 * cos(angle_rad) - df$PC02 * sin(angle_rad)
  df$PC2_rot <- df$PC01 * sin(angle_rad) + df$PC02 * cos(angle_rad)
  
  return(df)
}

# create PCA plot
create_pca_plot <- function(df, color_var, shape_var = NULL,
                            pc_variance, title, colors) {
  
  p <- ggplot(df, aes(x = PC1_rot, y = PC2_rot, colour = .data[[color_var]]))
  
  if (!is.null(shape_var)) {
    p <- p + aes(shape = .data[[shape_var]])
  }
  
  p +
    geom_point(size = 3) +
    stat_ellipse(aes(color = .data[[color_var]]), type = "t", size = 1) +
    scale_color_manual(
      values = colors,
      breaks = if (!is.null(names(colors))) names(colors) else waiver()
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.position = "right"
    ) +
    labs(
      title = title,
      x = paste0("PC1 (", round(pc_variance["PC01"], 1), "%)"),
      y = paste0("PC2 (", round(pc_variance["PC02"], 1), "%)")
    )
}

# scree plot
create_scree_plot <- function(eigenvalues) {
  ggplot(eigenvalues, aes(x = V1, y = V2)) +
    geom_bar(stat = "identity") +
    labs(title = "PCA eigenvalues") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(), 
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = -14)),
      plot.title.position = "panel",
      panel.grid = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0))
}

# load data
load_pca_data <- function(prefix) {
  eigenvalues <- read.csv(
    file.path(data_dir, paste0(prefix, ".eigenvalues.csv")),
    header = FALSE
  )
  
  eigenvectors <- read.csv(
    file.path(data_dir, paste0(prefix, ".eigenvectors.csv"))
  )
  
  list(values = eigenvalues, vectors = eigenvectors)
}

##### 1. all regions PCA #####

all_data <- load_pca_data("all_regions")

all_data$vectors <- rotate_pca(all_data$vectors, angle_deg = -180)

# % of variance explained by PC01 and PC02 obtained from the PCA eigenvalues
pc_variance_all <- c(PC01 = 10.7, PC02 = 7.2)

region_colors <- c(
  "Merritt Island" = "#d62728",
  "NBC" = "royalblue1",
  "CBC" = "navy",
  "SBC" = "chartreuse3",
  "Archbold" = "mediumorchid4",
  "Ocala" = "darkorange"
)

pca_all <- create_pca_plot(
  df = all_data$vectors,
  color_var = "Region",
  pc_variance = pc_variance_all,
  title = "PCA of individuals across sampling regions",
  colors = region_colors
)

scree_all <- create_scree_plot(all_data$values)

combined_all <- ggdraw() +
  draw_plot(pca_all) +
  draw_plot(scree_all, x = 0.04, y = 0.6, width = 0.3, height = 0.3) # adjust x, y based on scree plot position

##### 2. Brevard County PCA #####

bc_data <- load_pca_data("mainland_bc")

bc_data$vectors <- rotate_pca(bc_data$vectors, angle_deg = -90)

# % of variance explained by PC01 and PC02 obtained from the PCA eigenvalues
pc_variance_bc <- c(PC01 = 5.2, PC02 = 4.4)
 
# set the number of colors from Set1 palette based on the number of locations
bc_colors <- colorRampPalette(brewer.pal(8, "Set1"))(15)

bc_data$vectors$Region <- factor(
  bc_data$vectors$Region,
  levels = c("NBC", "CBC", "SBC")
)

bc_data$vectors[["Region: Location"]] <- paste(
  bc_data$vectors$Region,
  bc_data$vectors$Location,
  sep = ": "
)

bc_data$vectors[["Region: Location"]] <- factor(
  bc_data$vectors[["Region: Location"]],
  levels = c(
    "NBC: Buck Lake", "NBC: South Lake", "NBC: Fox Lake", "NBC: Dicerandra",
    "CBC: Viera","CBC: Cruickshank", "SBC: Malabar", "SBC: Valkaria", "SBC: Micco"
  )
)

pca_bc <- create_pca_plot(
  df = bc_data$vectors,
  color_var = "Region: Location",
  shape_var = "Region",
  pc_variance = pc_variance_bc,
  title = "PCA of individuals across mainland Brevard County",
  colors = bc_colors
)

##### 3. Central and South Brevard County PCA #####

cbc_sbc_data <- load_pca_data("cbc_sbc")

cbc_sbc_data$vectors <- rotate_pca(cbc_sbc_data$vectors, angle_deg = 90)

# % of variance explained by PC01 and PC02 obtained from the PCA eigenvalues
pc_variance_cbc_sbc <- c(PC01 = 5.0, PC02 = 3.6)

cbc_sbc_data$vectors[["Region: Location"]] <- paste(
  cbc_sbc_data$vectors$Region,
  cbc_sbc_data$vectors$Location,
  sep = ": "
)

cbc_sbc_data$vectors[["Region: Location"]] <- factor(
  cbc_sbc_data$vectors[["Region: Location"]],
  levels = c(
    "CBC: Viera","CBC: Cruickshank", "SBC: Malabar", "SBC: Valkaria", "SBC: Micco"
  )
)

pca_cbc_sbc <- create_pca_plot(
  df = cbc_sbc_data$vectors,
  color_var = "Region: Location",
  shape_var = "Region",
  pc_variance = pc_variance_cbc_sbc,
  title = "PCA of individuals across CBC and SBC",
  colors = bc_colors
)

##### combine all plots #####

final_plot <- plot_grid(
  combined_all,
  pca_bc,
  pca_cbc_sbc,
  ncol = 1,
  labels = c("a", "b", "c"),
  label_size = 18
)

##### save outputs #####

# all regions
ggsave(file.path(output_dir, "pca_all.pdf"), combined_all, width = 10, height = 6)
ggsave(file.path(output_dir, "pca_all.png"), combined_all, width = 10, height = 6, dpi = 600)
ggsave(file.path(output_dir, "pca_all.svg"), combined_all, width = 10, height = 6)

# Mainland BC
ggsave(file.path(output_dir, "pca_bc.pdf"), pca_bc, width = 10, height = 6)
ggsave(file.path(output_dir, "pca_bc.png"), pca_bc, width = 10, height = 6, dpi = 600)
ggsave(file.path(output_dir, "pca_bc.svg"), pca_bc, width = 10, height = 6)

# CBC, SBC
ggsave(file.path(output_dir, "pca_cbc_sbc.pdf"), pca_cbc_sbc, width = 10, height = 6)
ggsave(file.path(output_dir, "pca_cbc_sbc.png"), pca_cbc_sbc, width = 10, height = 6, dpi = 600)
ggsave(file.path(output_dir, "pca_cbc_sbc.svg"), pca_cbc_sbc, width = 10, height = 6)

# combined
ggsave(file.path(output_dir, "pca_combined.pdf"), final_plot, width = 10, height = 12)
ggsave(file.path(output_dir, "pca_combined.png"), final_plot, width = 10, height = 12, dpi = 600)
ggsave(file.path(output_dir, "pca_combined.svg"), final_plot, width = 10, height = 12)
