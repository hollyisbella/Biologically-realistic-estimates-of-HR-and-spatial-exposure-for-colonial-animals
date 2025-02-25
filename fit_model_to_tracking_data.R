################################################################################

# Model fitting to tracking data

################################################################################

# README #

# land polygons can be downloaded from:
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/

# Northern gannet GPS tracking data used in this file openly available from the
# Seabird Tracking Database: https://www.seabirdtracking.org/
# Datasets: Grassholm (ID: 731,732) and Great Saltee (ID:723)

# load libraries
library(stars)
library(sf)
library(dplyr)

################################################################################

# Define function - run_simulation
# run_simulation outputs likelihood used for model fitting

################################################################################

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
source("Functions/pre_optim.R")
source("Functions/process_posterior_run69.R")

run_simulation<-function(par, a_init, mask, colonies_sf, col_pos, tmax, pspill,
                         tel_points, tel_points_grid, num_tel){


  # load libraries
  library(stars)
  library(sf)
  library(dplyr)

  # assign parameters
  w=par[1]
  cap=par[2]
  m=par[3]
  a2=par[4]
  a3=par[5]

  # create permeability layer
  p<-create_p_layer(mask=mask, cap=cap) # create permeability layer

  # calculate distances to colonies avoiding land
  dist<-calculate_min_dist(a=a_init, mask=mask, col_pos=col_pos, avoid_land=T)
  dist[is.na(dist)]<-10000

  a_dist<-calculate_foraging(mask=mask, a=a_init, p=p, dist=dist, colonies_sf=colonies.sf,
                             col_pos=col_pos, tmax=tmax, m=m, pspill=pspill,
                             a2=a2, a3=a3, t_store_length=50,
                             saveloc="images/", savefreq=100,
                             plot_results=F, thresh_plot=0.0001, log_file_path="",
                             record_resids_path="")

  # calculate commuting surface
  commute<-calculate_commute2(a=a_dist, dist=dist, mask=mask, col_pos=col_pos, thresh_com=1*10^(-60),
                              save_loc="Images/",
                              save_freq=20, tmax=200, colonies_sf=colonies_sf,
                              sd_rand=0, plot_results=F, thresh_plot_a=0.0001)

  commuting<-commute[[1]]

  # NORMALISED surfaces
  norm<-normalise_foraging_commuting(a=a_dist, commuting=commuting)
  anorm<-norm[[1]]
  comnorm<-norm[[2]]

  # ESTIMATED HOME RANGE
  ud<-(1-w)*anorm+w*comnorm

  # add small probability to all space to prevent infinite LL calculation
  for (k in 1:length(ud)){
    ud[[k]]<-ud[[k]]+(min(ud[[k]][ud[[k]]>0], na.rm=T))/1000
  }

  LL<-c(rep(0,length(a_init)))

  # calculate log likelihood for each colony with data
  for (k in 1:length(a_init)){

    LL[k]<-sum(dpois(x=tel_points_grid[[k]],lambda=num_tel[k]*ud[[k]], log=TRUE))

  }

  # sum -ve LLs
  sum_LL<-sum(-LL)

  # save parameter values tried to text file
  # SPECIFY LOCATION
  sink("optimise/", append=T)

  print(paste(w, cap, m, a2, a3, sum_LL))

  sink()

  return(sum_LL)
}

################################################################################

# Run pre-optim

################################################################################

# read in colonies data
colonies<-read.csv("Data/posterior_run69.csv")

# rename columns e.t.c for use in other functions
colonies<-process_posterior_run69(colonies)

# read in land polygon shapefile from file location
land_poly<-sf::st_read("")

des_crs<-3035
og_crs<-4326
res<-5000
year<-112
forage<-400000
thinning<-20

# colonies with tracking data
colony_names<-c("Grassholm", "Great Saltee")

colonies_track<- dplyr::filter(colonies, Colony %in% colony_names)

# run functions required before optimisation
return_list<-pre_optim(land_poly=land_poly, colonies=colonies_track, year=year, og_crs=og_crs, des_crs=des_crs, max_forage=forage, res=res)

# extract pre-optim objects
a<-return_list[[1]]
mask<-return_list[[2]]
col_pos<-return_list[[3]]
colonies.sf<-return_list[[4]]

################################################################################

# Read in tracking data and process

################################################################################


# read in tracking data
grassholm_raw.sf<-read_sf("")
great_saltee_raw.sf<-read_sf("")

### Remove any locations within 1km of colony centre

# create 1km buffer around  colony
buffer_1km<-sf::st_buffer(colonies.sf, dist=1000)

# intersection of buffer with tracking points
intersect_col_buffer_grassholm<-st_intersection(buffer_1km[1,,],grassholm_raw.sf)
intersect_col_buffer_great_saltee<-st_intersection(buffer_1km[2,,],great_saltee_raw.sf)

'%!in%'<-Negate('%in%')

# filter out points 1km from colony centre
grassholm.sf<-grassholm_raw.sf %>% filter(geometry %!in% intersect_col_buffer_grassholm$geometry)
great_saltee.sf<-great_saltee_raw.sf %>% filter(geometry %!in% intersect_col_buffer_great_saltee$geometry)

# specify thinning level, so interval = 40 mins
thinning = 20

# sub-sample tracking data
thinned_values_grass<-seq(from=1, to=length(grassholm.sf$geometry), by=thinning)
thinned_values_saltee<-seq(from=1, to=length(great_saltee.sf$geometry), by=thinning)

thinned_grass<-grassholm.sf[thinned_values_grass,]
thinned_saltee<-great_saltee.sf[thinned_values_saltee,]

# create list of tracking data
tel_points<-list(thinned_grass, thinned_saltee)

# create stars rasters of tracking points
# create count column
tel_points[[1]]$count<-1
tel_points[[2]]$count<-1
tel_points_grid1<-st_rasterize(tel_points[[1]][,"count"],template=st_as_stars(st_bbox(mask),
                                                                     dx=res, dy=res), options=c("MERGE_ALG=ADD"))

tel_points_grid2<-st_rasterize(tel_points[[2]][,"count"],template=st_as_stars(st_bbox(mask),
                                                                      dx=res, dy=res), options=c("MERGE_ALG=ADD"))

tel_points_grid<-c(tel_points_grid1, tel_points_grid2)

# visualise gridded tracking data
image(tel_points_grid[[1]])
image(tel_points_grid[[2]])

# store total location points
num_tel<-c(rep(0, length(a)))

# assign no. location points to each colony
for (k in 1:length(a)){
  num_tel[k]<-length(tel_points[[k]]$geometry)
}

################################################################################

# Run model fitting

################################################################################

# IN OPTIM
tmax<-10000

pspill<-0.6

# create text file to store tried parameters in model fitting
#writeLines(c(""),"Optimise/log1.txt")

# Model fitting with adaptive search in parallel
# change number of cores(num_cores) as needed
library(PopED)
opt_ARS52<-optim_ARS(par=c(0.8, 4, 4, -0.9,-0.02), fn=run_simulation,
                     lower=c(0,0.01,0,-1,-0.1),upper=c(1,30,10,-0.82,0),
                     iter_adapt=40,
                     trace_iter=1, max_run=140, parallel=T, num_cores=11,
                     a_init=a, mask=mask,
                     colonies_sf=colonies.sf, col_pos=col_pos,
                     tmax=tmax, pspill=pspill, tel_points=tel_points,
                     tel_points_grid=tel_points_grid, num_tel=num_tel)


