################################################################################

# Calculate foraging ranges based on Grecian et al 2012 Methods Steps 1-3

################################################################################

# README #

# Grecian, W. J., Witt, M. J., Attrill, M. J., Bearhop, S., Godley,
# B. J., Grémillet, D., Hamer, K. C., & Votier, S. C. (2012).
# A novel projection technique to identify important at-sea areas for seabird
# conservation: An example using Northern gannets breeding in the North East Atlantic.
# Biological Conservation, 156, 43–52. https://doi.org/10.1016/j.biocon.2011.12.010

# land polygons can be downloaded from:
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/

################################################################################

# load libraries
library(dplyr)
library(sf)
library(stars)

# source functions
source("Functions/create_bbox.R")
source("Functions/create_colonies_sf.R")
source("Functions/create_landseamask.R")
source("Functions/select_colony_year.R")
source("Functions/process_posterior_run69.R")

# define foraging range grecian et al 2012
foraging_range<-function(P){
  return(0.344*P^(0.5)+40.062)
}

# year of interest
year<-2020

# read in colonies data
colonies<-read.csv("Data/posterior_run69.csv")

# rename columns for use in other functions
colonies<-process_posterior_run69(colonies)

# filter for year 2020 (1900= year 1)
year<-121

# original crs
og_crs<-4326

# desired crs
des_crs<-3035

# select colonies surrounding British Isles
ukire<-c("Ailsa Craig", "Barra Head", "Bass Rock", "Bempton Cliff", "Copinsay",
         "Fair Isle", "Flannan Isles", "Foula", "Grassholm","Great Saltee",
         "Hermaness", "Holy Isle Arran","Isle of May","Les Etacs", "Lundy", "Noss", "Ortac", "Rockall", "Scar Rocks",
         "Shiant Islands", "St Kilda", "St. Margarets", "Sula Sgeir", "Sule Skerry",
         "Sule Stack", "Westray", "Troup Head", "Bull Rock", "Clare Island",
         "Great Saltee", "Ireland's Eye", "Lambay", "Little Skellig", "Rouzic")


colonies_uk<- dplyr::filter(colonies, Colony %in% ukire)

# select specific colony year
colonies.year<-select_colony_year(colonies=colonies_uk, year=year)

# select necessary columns
colonies.year<-data.frame(ID=colonies.year$ID,median=colonies.year$median, latitude=colonies.year$latitude, longitude=colonies.year$longitude)

# create colonies .shp file
colonies.sf<-create_colonies_sf(colonies.year, og_crs, des_crs)

# define max foraging range
forage<-350000

# define resolution in units of des_crs
res<-5000

# create bbox arround colonies
bbox<-create_bbox(colonies_sf=colonies.sf, og_crs=og_crs, des_crs=des_crs, max_forage=forage)

# read in land polygon shapefile
land_poly<-st_read("")

land_polyt<-st_transform(land_poly, crs=3035)

# create landseamask from polygons shapefile
land_rast<-create_landseamask(bbox, des_crs=des_crs, res=res, land_poly=land_poly)

# set land to 0 and sea to 1
land_rast[land_rast==0]<-1
land_rast[land_rast==-999]<-0

# make data frame of colony populations
populations<-data.frame(colonies.sf$median)

# calculate foraging ranges
foraging_ranges<-foraging_range(populations)

# create foraging range buffers
buffer_sf<-st_buffer(colonies.sf, dist=foraging_ranges$colonies.sf.median*1000)

# calculate buffer area
buffer_area<-st_area(buffer_sf)

# convert buffers to stars
# need to create multi-layer stars object for each buffer
n<-dim(buffer_sf)[1]

# rasterise first buffer at 5km res
buffer_stars<-st_rasterize(buffer_sf[1,2,,],dy=5000, dx=5000)

# warp to bbox of all buffers
buffer_stars<-st_warp(buffer_stars,dest=land_rast)

# remove land areas within buffer
buffer_stars<-buffer_stars*land_rast

# calculate density of usage in buffer
buffer_stars<-buffer_stars*buffer_sf$median[1]/sum(buffer_stars[[1]], na.rm=T)

# rasterise buffers for all colonies
for (i in 2:n){

  # rasterise
  buffer_stars_temp<-st_rasterize(buffer_sf[i,2,,],dy=5000, dx=5000)

  # warp to bbox of all buffers
  buffer_stars_temp<-st_warp(buffer_stars_temp,dest=land_rast)

  # remove land areas within buffer
  buffer_stars_temp<-buffer_stars_temp*land_rast

  # calculate density of usage in buffer
  buffer_stars_temp<-buffer_stars_temp*buffer_sf$median[i]/sum(buffer_stars_temp[[1]], na.rm=T)

  # combine into one raster object, each layer is the buffer from one colony
  buffer_stars<-c(buffer_stars, buffer_stars_temp)
}

### save foraging range density surfaces to file
for (i in 1:n){

  write_stars(buffer_stars[i,,],paste("Data/", buffer_sf$ID[i],".nc",sep=""))
}




