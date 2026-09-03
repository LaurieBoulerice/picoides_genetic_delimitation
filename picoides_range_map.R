library(sf)
library(ggplot2)
library(ggpattern)


##########
#Dorsalis
##########

#these shapefiles were created in arcgis following birds of the world description of subspecies range
dorsalis_subsp<-st_read('/media/ssd/picoides_genetic_delimitation/Data/Raw/picoides_range/Dorsalis_subsp/P_dorsalis_subsp.shp')

d_subspecies<-c('bacatus','fasciatus','dorsalis')
dorsalis_subsp$d_subspecies<-d_subspecies
dorsalis_subsp<-dorsalis_subsp[,-(1:17)]

bacatus_sf <- dorsalis_subsp[dorsalis_subsp$d_subspecies == "bacatus", ]
fasciatus_sf <- dorsalis_subsp[dorsalis_subsp$d_subspecies == "fasciatus", ]
dorsalis_sf <- dorsalis_subsp[dorsalis_subsp$d_subspecies == "dorsalis", ]


mapview(bacatus_sf, col.regions = "#b2df8a", alpha = 0.6, layer.name = "Bacatus") +
  mapview(fasciatus_sf, col.regions = "#377eb8", alpha = 0.4, layer.name = "Fasciatus") +
  mapview(dorsalis_sf, col.regions = "#fb9a99", alpha = 0.4, layer.name = "Dorsalis")+
  mapview(p_dorsalis_sf)

###STATIC MAP###
#bounding box for xy limits
all_species_dorsalis <- rbind(
  bacatus_sf,
  fasciatus_sf,
  dorsalis_sf
)
bbox_dorsalis <- st_bbox(all_species_dorsalis)
dorsalis_bbox<-st_bbox(c(xmin = -163.08228, xmax = -52.63629, ymin = 34.07672  , ymax = 69.69269 )) #adjusted bbox

dorsalis_map<-ggplot() +
  geom_sf(data = world, fill = "white", color = "grey60") +
  geom_sf(data = all_species_dorsalis, aes(fill = d_subspecies), color = NA, alpha = 0.4) +
  coord_sf(
    xlim = c(-163.08228, -52.63629),
    ylim = c(34.07672, 69.69269),
    expand = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "bacatus" = "#b2df8a",  
      "fasciatus" = "#377eb8", 
      "dorsalis" = "#fb9a99"
    ),
    name = "Subspecies"
  ) 

#####BANG SUBSPECIES####

bang_americanus<-st_read('/media/ssd/picoides_genetic_delimitation/Data/Raw/picoides_range/Dorsalis_subsp/americanus/bang_americanus.shp')
bang_bacatus<-st_read('/media/ssd/picoides_genetic_delimitation/Data/Raw/picoides_range/Dorsalis_subsp/bacatus/bang_bacatus.shp')
bang_dorsalis<-st_read('/media/ssd/picoides_genetic_delimitation/Data/Raw/picoides_range/Dorsalis_subsp/dorsalis/bang_dorsalis.shp')
bang_labradorius<-st_read('/media/ssd/picoides_genetic_delimitation/Data/Raw/picoides_range/Dorsalis_subsp/labradorius/bang_labradorius.shp')

bang_subs <- bind_rows(
  bang_americanus %>% mutate(subspecies = "americanus"),
  bang_bacatus    %>% mutate(subspecies = "bacatus"),
  bang_dorsalis   %>% mutate(subspecies = "dorsalis"),
  bang_labradorius%>% mutate(subspecies = "labradorius")
)

all_species_dorsalis$d_subspecies <- factor(
  all_species_dorsalis$d_subspecies,
  levels = c("americanus", "bacatus", "dorsalis", "labradorius", "fasciatus")
)

bang_subs$subspecies <- factor(
  bang_subs$subspecies,
  levels = c("americanus", "bacatus", "dorsalis", "labradorius")
)

# Define stripe colors for Bang subspecies
stripe_colors <- c(
  "americanus" = "#1b9e77",
  "bacatus"    = "#d95f02",
  "dorsalis"   = "#7570b3",
  "labradorius"= "#e7298a"
)

# Plot
ggplot() +
  # Background map
  geom_sf(data = world, fill = "white", color = "grey60") +
  
  # Real subspecies polygons (filled colors)
  geom_sf(
    data = all_species_dorsalis,
    aes(fill = d_subspecies),
    color = NA,
    alpha = 0.4
  ) +
  
  # Bang subspecies polygons with stripes
  geom_sf_pattern(
    data = bang_subs,
    aes(pattern = subspecies, pattern_fill = subspecies),
    color = "black",
    fill = NA,                   # transparent base
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.03,
    size = 0.5
  ) +
  
  # Coordinate limits
  coord_sf(
    xlim = c(-163.08228, -52.63629),
    ylim = c(34.07672, 69.69269),
    expand = FALSE
  ) +
  
  # Fill scale for real subspecies (kept original colors)
  scale_fill_manual(
    name = "Subspecies",
    values = c(
      "bacatus"   = "#33a02c",  # deep green
      "fasciatus" = "#6a3d9a",  # purple
      "dorsalis"  = "#a6cee3"   # light blue
    )
  ) +
  
  # Pattern fill scale for Bang subspecies (striped legend)
  scale_pattern_fill_manual(
    name = "Bang's subspecies",
    values = stripe_colors
  ) +
  
  # Pattern type scale (needed for pattern legend)
  scale_pattern_manual(
    values = c(
      "americanus" = "stripe",
      "bacatus"    = "stripe",
      "dorsalis"   = "stripe",
      "labradorius"= "stripe"
    )
  ) 

