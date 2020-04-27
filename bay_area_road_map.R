########################################################################################################################
# Bay Area Road Map                                                                                                    #
#                                                                                                                      #
# Author: Dan Lependorf                                                                                                #
# Date: 2020-04-24                                                                                                     #
#                                                                                                                      #
# This script pulls a bunch of publicly available mapping data from the US Census Bureau's TIGER/Line Geodatabases and #
# generates a road map of the San Francisco Bay Area, with the landmass, water features, roads, and road coloration    #
# mapped to city boundaries.                                                                                           #
########################################################################################################################

library(tidyverse)
library(sf)
library(sp)
library(rgeos)
library(rgdal)
library(magick)

# Note that this script requires this lovely MapColoring package by Philipp Hunziker that isn't on CRAN. It can be
# downloaded from Github directly.
# remotes::install_github("hunzikp/MapColoring")
library(MapColoring)

########################################################################################################################
# Set parameters for the final image. Feel free to change any of these, just know that if you adjust one, you may have #
# to adjust some of the others to compensate. Expanding the dimensions, for example, may make the roads look too thin, #
# so go ahead and experiment with what looks best.                                                                     #
########################################################################################################################

# The subfolder in this directory to save the Census shapefiles in.
data_subfolder <- "./data/"

# For hydrography and roads, the data is organized in files per county. Download data for these county FIPS codes.
county_fips_codes <- c("Alameda" = "06001",
                       "Contra Costa" = "06013",
                       "Marin" = "06041",
                       "Napa" = "06055",
                       "San Francisco" = "06075",
                       "San Mateo" = "06081",
                       "Santa Clara" = "06085",
                       "Santa Cruz" = "06087",
                       "Solano" = "06095",
                       "Sonoma" = "06097")

# For city boundaries, the data is organized in a single file per state. Download data for these state FIPS codes.
state_fips_codes <- c("California" = "06")

# The width and height of the image, respectively.
image_dimensions <- c(9600, 12000)

# The latitude and longitude of the bounding box, bottom to top, left to right.
image_lat <- c(37.2, 38.3)
image_lon <- c(-122.7565, -121.6435)

# The thicknesses of the plotted roads. This requires two numbers, one for the thinnest style of road and one for the
# thickest. All other road styles will be scaled to match.
road_thickness <- c(0.5, 5)

########################################################################################################################
# Download and unzip everything from TIGER/Line.                                                                       #
########################################################################################################################

# Grab and unzip hydrography data for each of the counties listed above.
water_filenames <- paste0("tl_2019_", county_fips_codes, "_areawater.zip")
water_temp_files <- map_chr(water_filenames, tempfile)

if (dir.exists(file.path(data_subfolder, "water"))==FALSE) {
    dir.create(file.path(data_subfolder, "water"), recursive=TRUE)
}

walk2(.x=water_filenames, .y=water_temp_files,
      ~download.file(url=paste0("https://www2.census.gov/geo/tiger/TIGER2019/AREAWATER/", .x), destfile=.y))
walk(.x=water_temp_files,
     ~unzip(zipfile=.x, exdir=file.path(data_subfolder, "water")))

# Grab all municipal city boundaries for all of the states listed above.
cities_filenames <- paste0("tl_2019_", state_fips_codes, "_place.zip")
cities_temp_files <- map_chr(cities_filenames, tempfile)

if (dir.exists(file.path(data_subfolder, "cities"))==FALSE) {
    dir.create(file.path(data_subfolder, "cities"), recursive=TRUE)
}

walk2(.x=cities_filenames, .y=cities_temp_files,
      ~download.file(url=paste0("https://www2.census.gov/geo/tiger/TIGER2019/PLACE/", .x), destfile=.y))
walk(.x=cities_temp_files,
     ~unzip(zipfile=.x, exdir=file.path(data_subfolder, "cities")))

# Finally, download all road data for each county.
roads_filenames <- paste0("tl_2019_", county_fips_codes, "_roads.zip")
roads_temp_files <- map_chr(roads_filenames, tempfile)

if (dir.exists(file.path(data_subfolder, "roads"))==FALSE) {
    dir.create(file.path(data_subfolder, "roads"), recursive=TRUE)
}

walk2(.x=roads_filenames, .y=roads_temp_files,
      ~download.file(url=paste0("https://www2.census.gov/geo/tiger/TIGER2019/ROADS/", .x), destfile=.y))
walk(.x=roads_temp_files,
     ~unzip(zipfile=.x, exdir=file.path(data_subfolder, "roads")))


########################################################################################################################
# Load in 2019 area hydrography shapefiles from TIGER/Line and prep all water features.                                #
########################################################################################################################

# Load in all ".shp" shapefiles in the water subfolder and bind everything together.
water <- dir(file.path(data_subfolder, "water"), full.names=TRUE) %>%
    str_subset(".shp$") %>%
    # NOTE: If this is being run with an R version below 4.0.0, you will need to set stringsAsFactors=FALSE like so.
    # This applies to all st_read() calls in this script.
    # map(~st_read(.x, stringsAsFactors=FALSE)) %>%
    map(st_read) %>%
    do.call("rbind", .)

# The TIGER/Line shapefiles only cover water around the coast. Because of the shape of the Bay Area, the bottom-left
# corner is going to have a fairly large piece of empty ocean which isn't covered by these files. I'll have to draw the
# polygon myself.
water_ocean_fill <- c(-123, 37.87,
                      -122.65, 37.87,
                      -122.55, 37.8,
                      -122.55, 37.5,
                      -122.46, 37.45,
                      -122.46, 37.2,
                      -122.35, 37.13,
                      -122.35, 37,
                      -123, 37,
                      -123, 37.87) %>%
    matrix(ncol=2, byrow=TRUE) %>%
    list() %>%
    st_polygon() %>%
    st_sfc() %>%
    st_set_crs(4269)

########################################################################################################################
# Load in municipal boundary shapefiles and prep the cities data for plotting.                                         #
########################################################################################################################

cities <- dir(file.path(data_subfolder, "cities"), full.names=TRUE) %>%
    str_subset(".shp$") %>%
    st_read() %>%
    # The full dataset contains shapefiles for all Census-designated places in California, including a bunch of stuff I
    # don't want to plot, like unincorporated communities and military bases. I only want to keep actual incorporated
    # cities, so I'll need to filter on CLASSFP=="C1"
    filter(CLASSFP=="C1")

# This is where things get a little complicated and a little un-tidy. The MapColoring package I'm using is great, since
# it solves the four-color contiguous polygon map problem for me (see the "Four color theorem" Wikipedia page for more
# details), but it's old enough that it can't use the newer sf package for spatial data. This means I'll have to use the
# older rgdal package for this specific section. This readOGR() function is from rgdal and loads the shapefile as a
# SpatialPolygonsDataFrame, rather than as an sf object.
cities_sp <- dir(file.path(data_subfolder, "cities"), full.names=TRUE) %>%
    str_subset(".shp$") %>%
    readOGR()

# And this is where it unfortunately gets a little bit hackier. Aside from massive performance benefits, sf also has the
# upside of storing the spatial data in a dataframe-like object which has methods for all of the usual dplyr dataframe
# manipulation functions. This made the C1 filtering really easy. The sp version of this, no such luck. Since I'm just
# using this version of the city data to feed into the MapColoring package, I don't actually need to modify this
# properly and make sure it's actually plottable: I can just manually cut out the non-C1 data in the relevant S4 slots
# of this object, and MapColoring will be none the wiser.
cities_sp@polygons[which(cities_sp@data$CLASSFP!="C1")] <- NULL
cities_sp@data <- cities_sp@data[which(cities_sp@data$CLASSFP=="C1"),]

# Now, feed cities_sp into MapColoring::getColoring() to get city color assignments, ensuring that cities that border
# each other will never have the same color.
set.seed(17)
cities_color_assignments <- cities_sp %>%
    getColoring() %>%
    tibble(NAME=cities_sp@data$NAME,
           color_keys=.) %>%
    # I can get away with only four colors, except there's one single five that apparently needs to be assigned: Cudahy,
    # which is in southeast LA and is completely off of the map. Chuck it.
    filter(color_keys != 5) %>%
    # Mathematically, this is good, but visually it's not quite there. There's way more color #1 than anything else, and
    # the map ends up looking kind of one note. I'm going to randomly divide color assignment 1 into two equal subsets
    # to get a little more variety here. It's still not completely equal, but it's a lot closer.
    mutate(rng=sample(c(1, 5), size=n(), replace=TRUE),
           color_keys=ifelse(color_keys==1, rng, color_keys)) %>%
    # Also, there are a few cities that I'm going to adjust manually, because either they're the same color as cities
    # that are very close but not technically touching, or just because there's too much of that color in that region.
    mutate(color_keys=case_when(NAME=="Colma" ~ 5,
                                NAME=="Santa Clara" ~ 1,
                                NAME=="Berkeley" ~ 1,
                                NAME=="Concord" ~ 5,
                                NAME=="Pittsburg" ~ 4,
                                NAME=="Dublin" ~ 4,
                                NAME=="Woodside" ~ 4,
                                NAME=="Los Altos Hills" ~ 1,
                                NAME=="Emeryville" ~ 3,
                                TRUE ~ color_keys)) %>%
    select(-rng) %>%
    # Now that we have four colors of roughly equal proportion (and another fifth one in a small quantity), assign
    # colors from a custom color palette.
    mutate(color=case_when(color_keys==1 ~ "#4380BA", #blue
                           color_keys==2 ~ "#E47949", #orange
                           color_keys==3 ~ "#CD4169", #pink
                           color_keys==4 ~ "#E4A449", #yellow
                           color_keys==5 ~ "#319A72")) #green

# I'm also going to need these assignments as a named vector for ggplot2 later on.
cities_color_vector <- cities_color_assignments %>%
    select(color_keys, color) %>%
    distinct() %>%
    arrange(color_keys) %>%
    deframe()

# And finally, join these color assignments back to the main cities dataset.
cities <- left_join(cities, cities_color_assignments, by="NAME")

########################################################################################################################
# Read in all roads data and prepare it for plotting.                                                                  #
########################################################################################################################

roads <- dir(file.path(data_subfolder, "roads"), full.names=TRUE) %>%
    str_subset(".shp$") %>%
    map(st_read) %>%
    do.call("rbind", .) %>%
    # Much like filtering the city data, I only want to keep four road types: highways, primary roads, secondary roads,
    # and ramps. All other road types are things I don't want to plot, like unpaved trails and fire access roads.
    filter(MTFCC %in% c("S1100", "S1200", "S1400", "S1630")) %>%
    # Set thicknesses so the highways are thicker than the primary roads, etc. Note that these thicknesses aren't the
    # actual plot size values, since ggplot2 will scale them all anyway. This is just to make sure that they're all scaled
    # relative to each other.
    mutate(line_thickness=case_when(MTFCC=="S1100" ~ 5,
                                    MTFCC=="S1200" ~ 2,
                                    MTFCC %in% c("S1400", "S1630") ~ 1))

########################################################################################################################
# Plot! I have a complicated visual layering problem to solve: I want the intersection of the road and city layers,    #
# which will color roads that fall within city boundaries with the assigned color of that city. I can't do this while  #
# generating the main plot, since I want to isolate the road/city intersection and the main plot has water and         #
# background land layers to worry about, so I'm going to compute the road/city intersection first, and then overlay    #
# that on top of the rest of the map.                                                                                  #
########################################################################################################################

# First, generate the base layer, containing the land, water, and a full gray road layer to cover all roads that fall
# outside of city boundaries. Unfortunately, I don't know of a better way of doing this directly, so I'm going to have
# to generate each layer separately, save them to temp files, and then read them back in with magick to perform the
# layering.
base_layer_file <- tempfile()

png(filename=base_layer_file,
    width=image_dimensions[[1]],
    height=image_dimensions[[2]],
    units="px",
    type="cairo")
ggplot() +
    geom_sf(data=water, fill="#0D192E", color="#0D192E") +
    geom_sf(data=water_ocean_fill, fill="#0D192E", color="#0D192E") +
    geom_sf(data=roads, aes(size=line_thickness), color="#B3B1A1") +
    scale_size(range=road_thickness) +
    coord_sf(xlim=image_lon, ylim=image_lat, expand=FALSE) +
    theme_void() +
    theme(panel.background=element_rect(fill="#5C5C5C"),
          legend.position="none")
dev.off()

base_layer <- image_read(base_layer_file)

# Now generate the road and city layers.
road_layer_file <- tempfile()
city_layer_file <- tempfile()

png(filename=road_layer_file,
    width=image_dimensions[[1]],
    height=image_dimensions[[2]],
    units="px",
    type="cairo",
    bg="transparent")
ggplot() +
    geom_sf(data=roads, aes(size=line_thickness), color="#B3B1A1") +
    scale_size(range=road_thickness) +
    coord_sf(xlim=image_lon, ylim=image_lat, expand=FALSE) +
    theme_void() +
    theme(legend.position="none")
dev.off()

png(filename=city_layer_file,
    width=image_dimensions[[1]],
    height=image_dimensions[[2]],
    units="px",
    type="cairo",
    bg="transparent")
ggplot() +
    geom_sf(data=cities, aes(fill=as.character(color_keys), color=as.character(color_keys))) +
    scale_fill_manual(values=cities_color_vector) +
    scale_color_manual(values=cities_color_vector) +
    coord_sf(xlim=image_lon, ylim=image_lat, expand=FALSE) +
    theme_void() +
    theme(legend.position="none")
dev.off()

road_layer <- image_read(road_layer_file)
city_layer <- image_read(city_layer_file)

# This magick::image_flatten() function with the "In" operator generates the visual intersection of these two layers,
# pixel by pixel.
combined_road_city_layers <- image_flatten(c(road_layer, city_layer), operator="In")

# Layer everything together and save.
final_image <- base_layer %>%
    c(combined_road_city_layers) %>%
    image_flatten(operator="Over")

image_write(final_image, "./bay_area_road_map.png", format="png")

