################################################################################

# Estimate exposure to planned wind farms in 2020
# Using model-derived HRs, projected distributions and foraging ranges

################################################################################

# README #

# wind farm polygon data can be openly accessed from EMODnet Human Activities at (accessed on 2023-11-24):
# https://ows.emodnet-humanactivities.eu/geonetwork/srv/eng/catalog.search#/metadata/8201070b-4b0b-4d54-8910-abcea5dce57f

# code to produce projected distributions available from:
# D3. Critchley, E. J., Grecian, W. J., Kane, A., Jessopp, M. J., & Quinn, J. L. (2018).
# Marine protected areas show low overlap with projected distributions of seabird
# populations in Britain and Ireland. Biological Conservation, 224, 309–317.
# https://doi.org/10.1016/j.biocon.2018.06.007

################################################################################

# load libraries
library(sf)
library(stars)
library(dplyr)
library(gtools)

# source functions
source("Functions/select_colony_year.R")
source("Functions/process_posterior_run69.R")

################################################################################

# Initialise colonies data

################################################################################

# read in colonies data
colonies<-read.csv("Data/posterior_run69.csv")

# rename columns e.t.c for use in other functions
colonies<-process_posterior_run69(colonies)

# select year
year<-121

# set resolution used to estimate HRs
res<-5000

ukire<-c("Ailsa Craig", "Barra Head", "Bass Rock", "Bempton Cliff", "Copinsay",
         "Fair Isle", "Flannan Isles", "Foula", "Grassholm","Great Saltee",
         "Hermaness", "Holy Isle Arran","Isle of May","Les Etacs", "Lundy", "Noss", "Ortac", "Rockall", "Scar Rocks",
         "Shiant Islands", "St Kilda", "St. Margarets", "Sula Sgeir", "Sule Skerry",
         "Sule Stack", "Westray", "Troup Head", "Bull Rock", "Clare Island",
         "Great Saltee", "Ireland's Eye", "Lambay", "Little Skellig", "Rouzic")

# filter for British Isles colonies
colonies_uk<- dplyr::filter(colonies, Colony %in% ukire)

# select specific colony year
colonies.year<-select_colony_year(colonies=colonies_uk, year=year)

# arrange colonies in order
colonies.year_ordered<-arrange(colonies.year, colonies.year$ID)

################################################################################

# Reading in model-derived HRs (generated from "estimate_HR" see README.md for more details)

################################################################################

# ENTER FILE LOCATIONS
filenames_2020 <-list.files("", pattern="*home_range", full.names=TRUE)
filenames_2020<-mixedsort(filenames_2020)

ud<-read_stars(filenames_2020)


################################################################################

# Read in projected distributions

################################################################################

# ENTER FILE LOCATIONS
# code available at D3, see README above
filenames_proj <-list.files("",pattern="*distribution", full.names=TRUE)

filenames_proj<-mixedsort(filenames_proj)

proj_dist<-read_stars(filenames_proj)

################################################################################

# Reading in foraging ranges

################################################################################

# ENTER FILE LOCATIONS
# foraging ranges produced using "calculate_foraging_ranges.R" file
filenames_range<-list.files("", pattern="", full.names=TRUE)

filenames_range<-mixedsort(filenames_range)

forage_r<-read_stars(filenames_range)

################################################################################

# Reading in and sorting windfarm data

################################################################################

# select year for 'present day'
present<-2020

# read in wind farm data, ENTER FILE LOCATION
eu<-read_sf("Data/")

# check statuses
unique(eu$STATUS)

# make valid
eu<-st_make_valid(eu)

# transform to desired crs
eu<-st_transform(eu,crs=3035)

# filter for wind farms surrounding British Isles
eu<-eu %>%filter(COUNTRY %in% c("United Kingdom","Ireland","France","Belgium","Netherlands"))

# filter for planned, existing or in construction windfarms
eu_planned<-filter(eu, eu$STATUS %in% c("Planned", "Approved") |
                     (eu$STATUS %in% c("Production","Test site") & eu$YEAR>present))

# rasterise wind farms to HR raster size
eu_all.stars<-st_rasterize(eu_planned, dx=res,dy=res,template=st_as_stars(st_bbox(ud),
                                                                      dx=res, dy=res))
# set all windfarm presence to 1
all_wind<-eu_all.stars[1,,]
all_wind[all_wind>=1]<-1

# rasterise for projected distributions
all_wind_proj<-st_rasterize(eu_planned, dx=res, dy=res,
                            template=st_as_stars(st_bbox(proj_dist), dx=res,dy=res))
all_wind_proj<-all_wind_proj[1,,]
all_wind_proj[all_wind_proj>=1]<-1

# rasterise for range distributions
all_wind_range<-st_rasterize(eu_planned, dx=res, dy=res,
                            template=st_as_stars(st_bbox(forage_r), dx=res,dy=res))
all_wind_range<-all_wind_range[1,,]
all_wind_range[all_wind_range>=1]<-1

################################################################################

# Calculate Wind farm exposure

################################################################################

# create objects to store wind farm exposure for each colony
wind_ud_all<-ud*0
wind_proj_all<-proj_dist*0
wind_range_all<-forage_r*0

sum_effect_ud_all<-rep(0,(length(ud)-1))
sum_effect_proj_all<-rep(0,length(proj_dist))
sum_effect_range_all<-rep(0,length(forage_r))

# for each colony calculate exposure
for (k in 1:(length(ud))){
  wind_proj_all[[k]]<-proj_dist[[k]]*all_wind_proj[[1]]

  wind_range_all[[k]]<-forage_r[[k]]*all_wind_range[[1]]

  wind_ud_all[[k]]<-ud[[k]]*all_wind[[1]]

  sum_effect_proj_all[k]<-sum(wind_proj_all[[k]], na.rm=T)

  sum_effect_range_all[k]<-sum(wind_range_all[[k]], na.rm=T)

  sum_effect_ud_all[k]<-sum(wind_ud_all[[k]], na.rm=T)

}

# make dataframe of wind farm exposure
wind_exposure.df<-data.frame(method=c(rep("HRs",26), rep("Projected distributions",26),
                                      rep("Foraging range",26)),
                             Colony=rep(colonies.year_ordered,3),
                             absolute_exposure=c(wind_ud_all*colonies.year_ordered$median, wind_proj_all, wind_range_all),
                             percentage=scales::percent(round(c(sum_effect_ud_all,
                                                                (sum_effect_proj_all/colonies.year_ordered$median),
                                                                (sum_effect_range_all/colonies.year_ordered$median)),3)))

# filter for grassholm, great saltee and ailsa craig
wind_exposure_selected.df<-wind_exposure.df%>%filter(Colony%in% c("Ailsa Craig", "Grassholm", "Great Saltee"))

