#### Script Info ####
# Originally developed by: M Schrandt (April 3, 2020)
# Revisions by: M. Schram (August 7, 2026)

# Purpose: Assign the Tampa Bay sampling sites (gear 20) to the Tampa Bay Estuary Program Bay Segments
#          and remove sites outside the bay proper bay segments

#### Set-Up ####
library(tidyverse)
library(sf)          # classes and functions for vector data
library(spData)      # load geographic data
library(tmap)        # for static and interactive maps
library(here)

# Specify input catch file from the script "1_Query_FIM_Database"
phys <- 
  # Important Physical Master Data
  readRDS(
    here(
      "data/FIM_TB_PhysicalMaster.rds"
    )
  ) %>%
  # Filter data 
  filter(
    # Note: FIM_TB_PhysicalMaster.rds is already filtered on TB samples
    Bay == "TB",
    Type == "M",
    (Project_1 == "AM" | Project_2 == "AM" | Project_3 == "AM"),
    Gear %in% c("019", "020")
  )

FIMsites <- 
  readRDS(
    here(
      "data/FIM_TB_BiologyCounts.rds"
    )
  ) %>%
  filter(
    Reference %in% phys$Reference
  ) %>%
  left_join(
    select(
      phys,
      Reference,
      Project_1:Grid
    ),
    by = "Reference"
  ) %>%
  relocate(
    Project_1:Grid,
    .after = Month
  ) %>%
  select(
    Reference,
    Longitude,
    Latitude,
    Zone, 
    Grid
  ) %>%
  arrange(
    Reference
  ) %>%
  #remove duplicate reference numbers (since each species has its own row in the catch data)
  distinct()

# Map FIM sites and assign to TBEP management zones
# Read in the shape file of TBEP segments
shape <- 
  st_read(
    here(
      "data/TBEP_Segments_WGS84.shp"
    )
  ) %>%
  #clean up name of bay segment column
  mutate(
    bay_segment = case_when(
      BAY_SEGMEN == 1 ~ "OTB",
      BAY_SEGMEN == 2 ~ "HB",
      BAY_SEGMEN == 3 ~ "MTB",
      BAY_SEGMEN == 4 ~ "LTB",
      BAY_SEGMEN == 5 ~ "Boca Ceiga",
      BAY_SEGMEN == 6 ~ "Terra Ceia",
      BAY_SEGMEN == 7 ~ "Manatee River"
    )
  ) %>%
  #remove extra columns
  select(
    -FID_1, 
    -Count_, 
    -BAY_SEGMEN
  )

# Check out how it maps
plot(shape)

# Read in the shape file of Florida
shape.FL <- 
  st_read(
    here(
      "data/Florida.shp"
    )
  ) %>%
  select(geometry)

#Convert FIM site list into sf object
Sites_sf <- FIMsites %>%
  st_as_sf(
    coords = c("Longitude", "Latitude"),
    agr = "constant",
    crs = 4326,                # WGS84 projection, to match the shapefile projection
    stringsAsFactors = FALSE,
    remove = TRUE
  )

#Find points/sites within polygons of bay segments
#in 2022 I ran into an error of edge ## has duplicate vertex with edge ## and the following line is suggested fix
sf_use_s2(FALSE)

# perform the spatial join
site_in_seg <- 
  st_join(
    Sites_sf, 
    shape, 
    join = st_within
  )

#Get visual of sites
plot(Sites_sf)

# Look at how the assignments line up with the bay segments
# visualizing the map takes a little while, especially on VPN
map <- 
  tm_shape(
    shape
  ) +
  tm_polygons(
    "bay_segment",
    pallette = "div",
    id = "bay_segment"
  ) +
  #add in the sites, colored by what segment they've been assigned to
  tm_shape(
    site_in_seg
  ) +
  tm_symbols(
    size = 0.5, 
    col = "bay_segment",
    alpha = 1
  ) +
  #add in land on top so that we can really see the bay
  tm_shape(
    shape.FL
  ) +
  tm_fill(
    "gray"
  )

map

#Clean up some of the assignments; remove ones outside four main bay segments
TBEPSeg_assign <- 
  site_in_seg %>%
  #fix some of the assignments due to land changes
  #make sure grid 363 falls into LTB, and 203 & 220 into MTB, along with grid 204 with Latitude < 27.7937
  mutate(
    bay_segment = case_when(
      Grid == 363 ~ "LTB",
      Grid %in% c(203, 220) ~ "MTB",
      (Grid == 204 & unlist(map(.$geometry,1)) < 27.7937) ~ "MTB",
      TRUE ~ bay_segment
    )
  ) %>%
  #remove sites from outside boundary of bay proper (ones not to be included in the nekton index)
  filter(
    !bay_segment %in% c("Boca Ceiga", "Terra Ceia", "Manatee River")
    & !is.na(bay_segment)
  )

# save(TBEPSeg_assign, file = '../tbeptools/data/fimstations.RData', compress = 'xz')

#Preview all sites on map
map2 <- 
  tm_shape(
    shape
  ) +
  tm_polygons(
    "bay_segment", 
    pallette = "div", 
    id = "bay_segment"
  ) +
  tm_layout(
    main.title = "Tampa Bay Small Seine Samples Assigned to\nTBEP Segments (bay proper only",
    main.title.position = "center"
  ) +
  #add in the sites, colored by what segment they've been assigned to
  tm_shape(
    TBEPSeg_assign
  ) +
  tm_symbols(
    size = 0.5, 
    col = "bay_segment", 
    alpha = 1
  ) +
  #add in land on top so that we can really see the bay
  tm_shape(
    shape.FL
  ) +
  tm_fill(
    "gray"
  )

show(map2) #need show function in order for it to show up in viewer when sourced

#### Merge the TBEP Segment Assignment with the Catch Data ####
Catch_TBEPsegs <- 
  readRDS(
    here(
      "data/FIM_TB_BiologyCounts.rds"
    )
  ) %>%
  filter(
    Reference %in% phys$Reference
  ) %>%
  left_join(
    select(
      phys,
      Reference,
      Project_1:Grid
    ),
    by = "Reference"
  ) %>%
  right_join(
    TBEPSeg_assign, 
    by = c("Reference", "Zone", "Grid")
  ) %>%
  select(
    -geometry
  )

# write the TBEP Segment assignments to a .rds
saveRDS(
  Catch_TBEPsegs,
  file = here("data/TampaBay_NektonIndexTBEPSegAssign.rds")
)

