#### Script Info ####
# M Schrandt
# April 10, 2020

# Purpose: Map individual years of the TBNI

#### Set-Up ####
library(tidyverse)
library(sf)          # classes and functions for vector data
library(spData)      # load geographic data
library(tmap)        # for static and interactive maps
library(tmaptools)   # to crop one shape to another shape for bounding box of map

# Specify file name of the TBNI scores
input <- "TBNI_Scores_1998-2018.csv" #this comes from the script called "3_ProcessingCatchData"

# Specify the year you want to map
yr = 2018

# Read in score data and prep for mapping
scores1 <- read.csv(input, header = TRUE, stringsAsFactors = FALSE) %>%
  select(Reference, Year, bay_segment, Longitude, Latitude, TBNI_Score) %>%
  filter(Year == yr) %>%
  arrange(Reference)

AveBayScores <- scores1 %>%
  group_by(Year, bay_segment) %>%
  summarize(SegScore = round(mean(TBNI_Score))) %>%
  mutate(MapCol = case_when(SegScore <= 32 ~ "#FF0000",
                            SegScore >= 46 ~ "#008000",
                            TRUE ~ "#FFFF00"))

scores <- scores1 %>%
  left_join(AveBayScores, by = c("Year", "bay_segment"))

# Read in the shape file of TBEP segments
shape_TB <- read_sf(dsn = ".", layer = "TBEP_Segments_WGS84") %>%
  #clean up name of bay segment column
  mutate(bay_segment = case_when(.$BAY_SEGMEN == 1 ~ "OTB",
                              .$BAY_SEGMEN == 2 ~ "HB",
                              .$BAY_SEGMEN == 3 ~ "MTB",
                              .$BAY_SEGMEN == 4 ~ "LTB",
                              .$BAY_SEGMEN == 5 ~ "Boca Ceiga",
                              .$BAY_SEGMEN == 6 ~ "Terra Ceia",
                              .$BAY_SEGMEN == 7 ~ "Manatee River")) %>%
  #remove extra columns
  select(-FID_1, -Count_, -BAY_SEGMEN) %>%
  #filter for just the main bay segments
  filter(bay_segment %in% c("OTB", "HB", "MTB", "LTB")) %>%
  left_join(AveBayScores) %>%
  st_make_valid()

# Read in the shape file of Florida
shape_FL <- read_sf(dsn = ".", layer = "Florida_land 40k") %>%
  select(geometry) %>%
  st_make_valid()

#Convert TBNI scores and long/lat to sf object so we can map them
scores_sf <- scores %>%
  st_as_sf(
    coords = c("Longitude", "Latitude"),
    agr = "constant",
    crs = 4326,                # WGS84 projection, to match the shapefile projection
    stringsAsFactors = FALSE,
    remove = TRUE
  ) %>%
  mutate(MapCol = case_when(TBNI_Score <= 32 ~ "#FF0000",
                            TBNI_Score >= 46 ~ "#008000",
                            TRUE ~ "#FFFF00"))

#Need to specify our bounding box to be slightly larger (10%) than the FIM samples
bbox_new <- st_bbox(scores_sf)
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
bbox_new[1] <- bbox_new[1] - (0.25 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.1 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.1 * yrange) # ymax - top

bbox_new <- bbox_new %>%
  st_as_sfc() #make it a sf polygon

#Let's crop the other shape files to the size similar to the FIM points
shape_TB_crop <- tmaptools::crop_shape(x = shape_TB, y = bbox_new, polygon = TRUE)
shape_FL_crop <- tmaptools::crop_shape(x = shape_FL, y = bbox_new, polygon = TRUE)


#### Make Map(s) ####
Yr_map <- tm_shape(shape_TB_crop, bbox = bbox_new) +
  tm_polygons(col = "MapCol",
              alpha = 3/10,
              border.col = "darkgray",
              style = "fixed",
              breaks = c(32, 46, 100)) +
    tm_layout(legend.outside = FALSE) +
  tm_shape(shape_FL_crop, bbox = bbox_new) +
    tm_fill("lightgray") +
    tm_borders() +
  tm_shape(scores_sf, bbox = bbox_new) +
    tm_symbols(size = 0.3,
               alpha = 4.5/10,
               col = "MapCol") +
  tm_layout(frame = FALSE) +
  tm_credits("OTB", size = 1, col = "black", fontface = "bold", position = c(0.35, 0.74)) +
  tm_credits("HB", size = 1, col = "black", fontface = "bold", position = c(0.7, 0.62)) +
  tm_credits("MTB", size = 1, col = "black", fontface = "bold", position = c(0.47, 0.46)) +
  tm_credits("LTB", size = 1, col = "black", fontface = "bold", position = c(0.26, 0.24)) +
  tm_layout(title = paste0(max(shape_TB_crop$Year)), title.position = c("center", "top"),
            title.fontface = 2, title.size = 1.5) +
  tm_scale_bar(breaks = c(0,5,10,15), text.size = 1) +
  tm_compass(type = "arrow", position = c("left", "bottom"))
Yr_map



