library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)

##########
#load data
##########

#ice sheet
url <- "https://services.arcgis.com/XSeYKQzfXnEgju9o/arcgis/rest/services/Last_Glacial_Maximum/FeatureServer/0/query?where=1%3D1&outFields=*&f=geojson"
lgm <- st_read(url)

#base map
world <- ne_countries(
       scale = "medium",
       returnclass = "sf"
   )
world <- st_transform(world, st_crs(lgm))

###################
##coastline at LGM
###################

bathymetryFile <- "https://raw.githubusercontent.com/jorgeassis/rGIS/master/Data/BathymetryDepthMean.tif"

bathymetry <- rast(bathymetryFile)

lgm_raster <- ifel(
  bathymetry >= -120,
  1,
  NA
)

lgm_landmass <- as.polygons(
  lgm_raster,
  dissolve = TRUE,
  na.rm = TRUE
)

lgm_landmass <- st_as_sf(lgm_landmass)

lgm_landmass <- st_transform(
  lgm_landmass,
  4326
)

###########################
#boreal forest delimitation
###########################

boreal_forest<-st_read("/media/ssd/picoides_genetic_delimitation/Data/Raw/boreal/NABoreal.shp")

boreal_forest <- st_transform(
  boreal_forest,
  4326
)

########
#plot
#######

#Plot LGM

ggplot() +
  geom_sf(
    data = world,
    fill = "grey90",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_sf(
    data = lgm_landmass,
    fill='#E3C98A',
    color = "#E3C98A",
    linewidth = 0.3
  ) +
  geom_sf(
    data = lgm,
    fill = "lightblue",
    color = "lightblue",
    alpha = 0.4
  ) +
  coord_sf(
    xlim = c(-180, -40),
    ylim = c(30, 85),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )

#Plot boreal forest delimitation
ggplot() +
  geom_sf(
    data = world,
    fill = "grey90",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_sf(
    data = boreal_forest,
    fill='#7FB77E',
    color = "#7FB77E",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(-180, -40),
    ylim = c(30, 85),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )

