################################################################################

# Estimate home ranges
# For each run: only change the following as needed:
# the year, the chosen colonies/regions and the saving location

################################################################################

# README #

# land polygons can be downloaded from:
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/

################################################################################

# load libraries
library(sf)
library(dplyr)
library(tidyverse)
library(stars)

# source functions
source("Functions/calculate_commute2.R")
source("Functions/calculate_foraging.R")
source("Functions/calculate_min_dist.R")
source("Functions/create_bbox.R")
source("Functions/create_colonies_sf.R")
source("Functions/create_colonies_stars_array.R")
source("Functions/create_landseamask.R")
source("Functions/create_p_layer.R")
source("Functions/move_colonies2.R")
source("Functions/normalise_foraging_commuting.R")
source("Functions/select_colony_year.R")
source("Functions/process_posterior_run69.R")

# read in colonies data
colonies<-read.csv("Data/posterior_run69.csv")

# rename columns e.t.c for use in other functions
colonies<-process_posterior_run69(colonies)

# desired crs
des_crs<-3035

# original crs
og_crs<-4326

# resolution of grid in units of des_crs
res<-5000

# CHOOSE YEAR (1 = year 1900, i.e. 112= year 2011)
year<-112

# define max foraging range for cropping grid
forage<-350000

# SELECT COLONIES
# British Isles colonies
col_regions<-c("Celtic Seas", "Minches & Western Scotland",
               "Irish Sea", "Scottish Continental Shelf","Northern North Sea",
               "South West Atlantic")

# if selecting specific colonies
#select_cols<-c("Grassholm","Great Saltee")

# select colonies in desired regions
colonies_reg<- dplyr::filter(colonies, regional_seas_split %in% col_regions)
#colonies_reg<- dplyr::filter(colonies, Colony %in% select_cols)

# select specific colony year
colonies.year<-select_colony_year(colonies=colonies_reg, year=year)

# calculate real year
colonies.year$colonisation_year<-colonies.year$ColonisationYear-1899
colonies.year$extinction_year<-colonies.year$ExtinctionYear-1899

# filter out uncolonised colonies
colonies.year<-filter(colonies.year, colonisation_year<=year)

# filter out extinct colonies
colonies.year<-filter(colonies.year, is.na(extinction_year)|extinction_year>year)

# arrange colonies in order of ID
colonies.year<-arrange(colonies.year, ID)

# create colonies .shp file
colonies.sf<-create_colonies_sf(colonies.year, og_crs, des_crs)

# create bbox around colonies account for forage range
bbox<-create_bbox(colonies_sf=colonies.sf, max_forage=forage)

# read in land polygon shapefile from file location
land_poly<-sf::st_read("")

# create landseamask from polygons shapefile
# MAKE THIS MORE FLEXIBLE - just being able to specify permitted and non-permitted space
land_rast<-create_landseamask(bbox, des_crs=des_crs, res=res, land_poly=land_poly)

# visualise landseamask
plot(land_rast)

# put colony populations into gridded space
list<-create_colonies_stars_array(colonies_sf=colonies.sf, landseamask=land_rast, res=res)

# extract function results (a array and central colony locations)
a<-list[[1]]
col_pos<-list[[2]]

# check if colonies on land and if so move
updated_positions<-move_colonies2(landseamask=land_rast, a=a,col_pos=col_pos)

# extract updated colony positions
a_update<-updated_positions[[1]]
col_pos_update<-updated_positions[[2]]

# true parameters (from fitting to simulated data, see 'fit_model_to_simulated_data.R')
true_par<-c(0.6, 15, 1, -0.9, -0.02)
# estimated parameters (from fitting to simulated data, see 'fit_model_to_simulated_data.R')
est_par<-c(0.51, 11.5, 3.48,-0.83, -0.028)
# estimated parameters from fitting to tracking data (see 'fit_model_to_tracking_data.R')
tracking_par<-c(0.5, 1.25, 1.87, -0.92, -0.0354)

# choose parameters to estimate HRs with
par<-tracking_par

# assign carrying capacity of sea
cap<-par[2]

# create permeability layer
p<-create_p_layer(mask=land_rast, cap=cap)

# calculate distances avoiding land
dist<-calculate_min_dist(a=a_update, mask=land_rast, col_pos=col_pos_update, avoid_land=T)
dist[is.na(dist)]<-10000

# DEFINE SAVING LOCATIONS FOR progress images and final outputs
save_im<-paste("Data/")
save_tif<-paste("Data/")

# run model
for (count in 1:1){

  start<-Sys.time()

  # calculate foraging surface
  a<-calculate_foraging(mask=land_rast, a=a_update, p=p, dist=dist, colonies_sf=colonies.sf,
                        col_pos=col_pos_update, tmax=10000, m=par[3], pspill=0.6,
                        a2=par[4], a3=par[5], t_store_length=50,
                        saveloc=save_im, savefreq=100,
                        plot_results=T, thresh_plot=0.0001, log_file_path="",
                        record_resids_path="")

  # save foraging surface
  for (k in 1:length(a)){
    write_stars(a, layer=k, paste(save_tif,"foraging_usage", colonies.sf$ID[k],".nc", sep=""))
  }

  # calculate commuting surface
  commute<-calculate_commute2(a=a, dist=dist, mask=land_rast, col_pos=col_pos_update,thresh_com=1*10^(-8),
                              save_loc=save_im,
                              save_freq=10, tmax=200, colonies_sf=colonies.sf,
                              plot_results=T, thresh_plot_a=0.0001)

  # extract commuting
  commuting<-commute[[1]]

  # extract final reversed positions
  reverse_a<-commute[[2]]

  # save commuting surface
  for (k in 1:length(a)){
    write_stars(commuting, layer=k, paste(save_tif,"commuting_usage",colonies.sf$ID[k],".nc", sep=""))
  }

  # NORMALISED surfaces
  norm<-normalise_foraging_commuting(a=a, commuting=commuting)
  anorm<-norm[[1]]
  comnorm<-norm[[2]]

  # weighting
  w<-par[1]

  # ESTIMATED HR
  ud<-(1-w)*anorm+w*comnorm

  # save HR
  for (k in 1:length(a)){
    write_stars(ud, layer=k, paste(save_tif,"home_range", colonies.sf$ID[k],".nc", sep=""))
  }

  # save colony information
  st_write(colonies.sf,paste(save_tif, "colonies_sf_.shp",sep=""))

  end<-Sys.time()
}


