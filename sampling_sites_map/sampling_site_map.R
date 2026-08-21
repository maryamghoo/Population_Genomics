# load libraries
library(raster)      # create spatial extents and crop geographic data
library(rworldmap)   # obtain high-resolution world/country map boundaries
library(dplyr)       # data manipulation and dataframe processing
library(ggplot2)     # plotting
library(sf)          # handle, transform, and project spatial vector data
library(ggspatial)   # add spatial features such as scale bars and north arrows
library(grid)        # control graphical objects (grobs), sizes, and positioning

# set paths

data_dir <- "r_data"
output_dir <- "r_results"

# create output directory if it does not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##### 1. define map boundaries

# Geographic extent of the main study area
main_boundary <- extent(
  -82, -80.5,
  27, 29.68)

# Geographic extent of Florida for the inset map
florida_boundary <- extent(
  -87.6349, -80.031,
  24.3963, 31)


##### 2. load base map

# get high-resolution geographic boundaries
world_map <- getMap(resolution = "high")


##### 3. create Florida inset map

# crop the world map to Florida
florida_map_data <- crop(
  world_map,
  y = florida_boundary) %>%
  fortify()

# the red rectangle indicates the extent of the main study area
inset_map <- ggplot() +
  geom_polygon(
    data = florida_map_data,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_rect(
    aes(
      xmin = -82,
      xmax = -80.5,
      ymin = 27,
      ymax = 29.5
    ),
    fill = NA,
    colour = "red",
    linewidth = 0.75
  ) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = NA,
      colour = NA
    )
  )


##### 4. prepare main map

# crop the base map to the main study area
main_sf <- st_as_sf(
  crop(world_map, y = main_boundary)
)

# project the study area to UTM Zone 17N (EPSG:32617).
# UTM uses meters, allowing the scale bar to represent real distances.
main_sf_utm <- st_transform(
  main_sf,
  crs = 32617)


##### 5. read and project sampling locoations

# read sampling-site coordinates
data2 <- read.csv(file.path(data_dir, "coordinates.csv")) %>%
  mutate(
    longitude = as.numeric(LatX),
    latitude = as.numeric(LatY))

# convert sampling locations from longitude/latitude (WGS84) to UTM Zone 17N so they match the main map projection
data2_utm <- data2 %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326) %>%
  st_transform(
    crs = 32617)


##### 6. convert spatial data to X/Y coordinates

main_coords <- st_coordinates(main_sf_utm)

main_map_data_utm <- data.frame(
  X = main_coords[, "X"],
  Y = main_coords[, "Y"],
  group = main_coords[, "L1"])

# extract UTM coordinates for sampling locations
point_coords <- st_coordinates(data2_utm)

data2_plot <- data.frame(
  X = point_coords[, "X"],
  Y = point_coords[, "Y"])


##### 7. convert map display limits to UTM

# define the desired plotting extent in longitude/latitude
limit_box <- st_as_sfc(
  st_bbox(
    c(
      xmin = -82.1,
      xmax = -80.1,
      ymin = 26.85,
      ymax = 29.75
    ),
    crs = st_crs(4326))
)

# transform the plotting extent to UTM Zone 17N
limit_box_utm <- st_transform(
  limit_box,
  crs = 32617
)

# extract projected limits for use in coord_fixed()
utm_bbox <- st_bbox(limit_box_utm)


##### 8. create main map

main_map <- ggplot() +
  # main study-area polygon
  geom_polygon(
    data = main_map_data_utm,
    aes(
      x = X,
      y = Y,
      group = group
    ),
    fill = "grey",
    colour = "black",
    linewidth = 0.5
  ) +
  # Sampling locations
  geom_point(
    data = data2_plot,
    aes(
      x = X,
      y = Y
    ),
    colour = "black",
    size = 3
  ) +
  # automatic metric scale bar
  annotation_scale(
    location = "br",
    width_hint = 0.35,
    style = "ticks",
    unit_category = "metric",
    plot_unit = "m",
    text_cex = 0.8,
    line_width = 0.7,
    pad_x = unit(3, "cm"),
    pad_y = unit(1.1, "cm")
  ) +
  # north arrow
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm")
  ) +
  # set plotting extent and desired map appearance
  coord_fixed(
    ratio = 0.8,
    xlim = c(
      as.numeric(utm_bbox["xmin"]),
      as.numeric(utm_bbox["xmax"])
    ),
    ylim = c(
      as.numeric(utm_bbox["ymin"]),
      as.numeric(utm_bbox["ymax"])
    ),
    expand = FALSE
  ) +
  # remove axes, labels, and gridlines
  theme_void() +
  # set water/background color and remove plot margins
  theme(
    panel.background = element_rect(
      fill = "lightsteelblue2"
    ),
    plot.margin = unit(
      c(0, 0, 0, 0),
      "cm"
    )
  )


##### 9. add Florida inset

# convert the Florida inset to a graphical object (grob)
inset_grob <- ggplotGrob(inset_map)

# overlay the inset on the upper-right portion of the main map.
gg <- main_map +
  annotation_custom(
    grob = inset_grob,
    xmin = as.numeric(utm_bbox["xmin"]) +
      0.65 * (
        as.numeric(utm_bbox["xmax"]) -
          as.numeric(utm_bbox["xmin"])
      ),
    xmax = as.numeric(utm_bbox["xmin"]) +
      1.00 * (
        as.numeric(utm_bbox["xmax"]) -
          as.numeric(utm_bbox["xmin"])
      ),
    ymin = as.numeric(utm_bbox["ymin"]) +
      0.67 * (
        as.numeric(utm_bbox["ymax"]) -
          as.numeric(utm_bbox["ymin"])
      ),
    ymax = as.numeric(utm_bbox["ymin"]) +
      0.92 * (
        as.numeric(utm_bbox["ymax"]) -
          as.numeric(utm_bbox["ymin"])
      )
  )


##### 10. save final map as svg for additional labeling and adding polygins manually

ggsave(
  filename = file.path(output_dir, "raw_sampling_map.svg"),
  plot = gg,
  width = 8,
  height = 7.3,
  units = "in"
)
