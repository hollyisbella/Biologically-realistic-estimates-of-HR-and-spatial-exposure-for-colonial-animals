################################################################################

### Model fitting to simulated data

################################################################################

# README #
# land polygons can be downloaded from:
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/

# load libraries
library(sf)
library(dplyr)
library(PopED)

################################################################################

# Define functions - generate_simulated_data and run_simulation
# generate_simulated_data: generates simulated data for use in model fitting
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

# FUNCTION to generate simulated data
generate_simulated_data<-function(colonies, year, og_crs, des_crs, max_forage, res){

  # select specific colony year
  colonies.year<-select_colony_year(colonies=colonies_wales, year=year)

  # create colonies .shp file
  colonies.sf<-create_colonies_sf(colonies.year, og_crs, des_crs)

  # calculate bbox function
  bbox<-create_bbox(colonies_sf=colonies.sf, max_forage=forage)

  # read in land polygon shapefile
  land_poly<-st_read("")

  # create landseamask from polygons shapefile
  land_rast<-create_landseamask(bbox=bbox, des_crs=des_crs, res=res, land_poly=land_poly)

  # create colonies raster and get colony center cell locations
  list<-create_colonies_stars_array(colonies_sf=colonies.sf, landseamask=land_rast, res=res)

  a<-list[[1]]
  col_pos<-list[[2]]

  # specify carrying capacity
  cap<-15

  # create permeability layer
  p<-create_p_layer(mask=land_rast, cap=cap)

  # calculate distances avoiding land
  dist<-calculate_min_dist(a=a, mask=land_rast, col_pos=col_pos, avoid_land=T)
  dist[is.na(dist)]<-10000

  # DEFINE SAVING LOCATIONS FOR progress images and final outputs
  save_im<-paste("Data/")

  a_dist<-calculate_foraging(mask=land_rast, a=a, p=p, dist=dist, colonies_sf=colonies.sf,
                             col_pos=col_pos, tmax=10000, m=1, pspill=0.6,
                             a2=-0.9, a3=-0.02, t_store_length=50,
                             saveloc=save_im, savefreq=100,
                             plot_results=T, thresh_plot=0.0001, log_file_path="",
                             record_resids_path="")

  # COMMUTING SURFACE
  commute<-calculate_commute2(a=a_dist, dist=dist, mask=land_rast,
                             col_pos=col_pos,thresh_com=1*10^(-8),
                             save_loc=save_im,
                             save_freq=10, tmax=500, colonies_sf=colonies.sf,
                            plot_results=F, thresh_plot_a=0.0001)

  commuting<-commute[[1]]

  # NORMALISATION
  norm<-normalise_foraging_commuting(a=a_dist, commuting=commuting)
  anorm<-norm[[1]]
  comnorm<-norm[[2]]

  # SIMULATED HOME RANGES
  w<-0.6
  ud<-(1-w)*anorm+w*comnorm

  ud1<-ud[1,,]
  ud2<-ud[2,,]

  # get size of grids
  x_size1<-unname(dim(ud1))[1]
  y_size1<-unname(dim(ud1))[2]
  x_size2<-unname(dim(ud2))[1]
  y_size2<-unname(dim(ud2))[2]

  # generate simulated data
  a1_data<-matrix(rpois(x_size1*y_size1, colonies_sf$median[[1]]*ud1[[1]]), nrow=x_size1, ncol=y_size1)
  a2_data<-matrix(rpois(x_size2*y_size2, colonies_sf$median[[2]]*ud2[[1]]), nrow=x_size2, ncol=y_size2)

  # return values
  return_list<-list(a_dist, a, land_rast, col_pos, colonies.sf, a1_data, a2_data)

  return(return_list)
}


################################################################################

run_simulation<-function(par, a_init, mask, colonies_sf, col_pos, tmax, pspill,
                         a1_data, a2_data){

  # load libaries
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
  p<-create_p_layer(mask=mask, cap=cap)


  # calculate distances avoiding land
  dist<-calculate_min_dist(a=a_init, mask=mask, col_pos=col_pos, avoid_land=T)
  dist[is.na(dist)]<-10000

  # calculate foraging surface
  # check image path if printing
  a_dist<-calculate_foraging(mask=mask, a=a_init, p=p, dist=dist, colonies_sf=colonies.sf,
                            col_pos=col_pos, tmax=tmax, m=m, pspill=pspill,
                            a2=a2, a3=a3, t_store_length=50,
                            saveloc="images/", savefreq=100,
                            plot_results=F, thresh_plot=0.0001, log_file_path="",
                            record_resids_path="")

  # calculate commuting surface
  commute<-calculate_commute2(a=a_dist, dist=dist, mask=mask, col_pos=col_pos, thresh_com=1*10^(-8),
                              save_loc="Images/",
                              save_freq=20, tmax=500, colonies_sf=colonies_sf,
                              plot_results=F, thresh_plot_a=0.0001)

  commuting<-commute[[1]]

  # NORMALISE surfaces
  norm<-normalise_foraging_commuting(a=a_dist, commuting=commuting)
  anorm<-norm[[1]]
  comnorm<-norm[[2]]

  # ESTIMATED HOME RANGE
  ud<-(1-w)*anorm+w*comnorm

  # add small probability to all space to prevent infinite LL calculation
  for (k in 1:length(ud)){
    ud[[k]]<-ud[[k]]+(min(ud[[k]][ud[[k]]>0], na.rm=T))/1000
  }

  # calculate log likelihoods
  LL1<-sum(dpois(x=a1_data,lambda=colonies_sf$median[[1]]*ud[[1]], log=TRUE))
  LL2<-sum(dpois(x=a2_data,lambda=colonies_sf$median[[2]]*ud[[2]], log=TRUE))

  # sum -ve LLs
  sum_LL<--LL1-LL2

  # write tried parameters to text file
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

# set orginal and desired crs
des_crs<-3035
og_crs<-4326

# set resolutions (units of des_crs)
res<-5000

# set year of interest
year<-112
forage<-200000

# define colonies to validate with simulated data
wales<-c("Grassholm","Great Saltee")

# select colonies
colonies_wales<- dplyr::filter(colonies, Colony %in% wales)

# run function to create simulated data
sim_data<-generate_simulated_data(colonies=colonies_wales, year=year, og_crs=og_crs, des_crs=des_crs, max_forage=forage, res=res)

# run function to generate objects needed as input to optimisation
return_list<-pre_optim(land_poly=land_poly, colonies=colonies_wales, year=year, og_crs=og_crs, des_crs=des_crs, max_forage=forage, res=res)

# extract simulated data
a1_data<-sim_data[[6]]
a2_data<-sim_data[[7]]

# extract objects from pre-optim needed for optimisation function
a<-return_list[[1]]
mask<-return_list[[2]]
col_pos<-return_list[[3]]
colonies.sf<-return_list[[4]]

# write simulated data to file
write.csv(a1_data,"optimise/a1_data_ap10.csv")
write.csv(a2_data,"optimise/a2_data_ap10.csv")

# if simulated data already produced, read in from file
#a1_data<-as.matrix(read.csv("Data/a1_data_ap10.csv"))[,2:93]
#a2_data<-as.matrix(read.csv("Data/a2_data_ap10.csv"))[,2:93]

################################################################################

# Run model fitting

################################################################################

# max iterations
tmax<-10000

# set spillover proportion
pspill<-0.6

# fit model to simulated data
# change number of cores(num_cores) as needed
opt_ARS<-optim_ARS(par=c(0.4, 10, 0.4, -0.85,-0.015), fn=run_simulation,
                   lower=c(0,0.01,0,-1,-0.1),upper=c(1,30,10,-0.82,0),
                   iter_adapt=40, num_cores=11,
                   trace_iter=1, max_run=140, parallel=T,
                   a_init=a, mask=mask,
                   colonies_sf=colonies.sf, col_pos=col_pos,
                   tmax=tmax, pspill=pspill,
                   a1_data=a1_data, a2_data=a2_data)



