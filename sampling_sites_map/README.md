# Sampling Location Map

## Overview
This workflow generates a geographic map showing the sampling locations used in the study. The map focuses on the main study area in Florida and includes a Florida inset map showing the geographic extent of the study region.
The workflow uses geographic boundary data together with sampling-site coordinates and projects the main map to UTM Zone 17N for accurate spatial representation and distance scaling.

The workflow includes:

- Loading high-resolution geographic boundary data
- Defining the main study area and Florida inset boundaries
- Creating a Florida inset map
- Importing sampling-site coordinates
- Converting geographic coordinates to spatial data
- Projecting the study area and sampling locations to UTM Zone 17N
- Plotting sampling locations on the study-area map
- Adding a metric scale bar and north arrow
- Adding a Florida inset showing the location of the main study area
- Exporting the final map in SVG format for additional manual editing

See:
`sampling_site_map.R`

### Directory Structure
```
project/
│
├── r_data/
│   └── coordinates.csv
│
├── r_results/
│   └── raw_sampling_map.svg
│
├── sampling_site_map.R
└── README.md
```

## 1. R Sampling Location Map Workflow

### Purpose

The R workflow generates a geographic map of the study area and sampling locations in Florida.

### R Packages
The script requires:

- `raster`
- `rworldmap`
- `dplyr`
- `ggplot2`
- `sf`
- `ggspatial`
- `grid`

### Expected Input
The sampling-site coordinate table can be found in the `r_data/` folder:
`coordinates.csv`
The current script expects the following columns:
- `LatX` — longitude
- `LatY` — latitude
Coordinates are expected to be provided in longitude/latitude using the WGS84 coordinate reference system (EPSG:4326).

Example:
```
LatX,LatY
-81.35,28.62
-81.12,28.91
-80.95,27.84
```

### Output
`raw_sampling_map.svg`
saved in: `r_results/`. SVG format is useful for opening the map in vector graphic editors such as Inkscape, allowing manual adjustment, labeling, and addition of polygons or other graphical elements without loss of image quality.

### Usage

- Place the sampling coordinate file in:
```
r_data/coordinates.csv
```
- Make sure the coordinate table contains the expected `LatX` and `LatY` columns.
- Update the following in `sampling_site_map.R` if needed:
`data_dir <- "r_data"`
`output_dir <- "r_results"`
- Adjust the main study-area extent if needed:
```
main_boundary <- extent(
  -82, -80.5,
  27, 29.68
)
```
- Adjust the Florida inset extent if needed:
```
florida_boundary <- extent(
  -87.6349, -80.031,
  24.3963, 31
)
```
- Adjust the final map display limits if needed:
```
xmin = -82.1
xmax = -80.1
ymin = 26.85
ymax = 29.75
```
- Adjust sampling-point size, map colors, scale-bar position, north-arrow position, inset-map position, and figure dimensions as needed.
