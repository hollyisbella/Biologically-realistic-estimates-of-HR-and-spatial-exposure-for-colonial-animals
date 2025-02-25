################################################################################

# Validate model HRs with tracking data

################################################################################

# README #

# Northern gannet GPS tracking data used in this file openly available from the
# Seabird Tracking Database: https://www.seabirdtracking.org/
# Datasets: Ailsa Craig (ID: 716), Bass Rock (ID:718), Sule Skerry (ID: 719),
# Bull Rock (ID: 720), Lambay (ID: 724), Little Skellig (ID: 721),
# Les Etacs (ID:733), Ile Rouzic (ID:734)

# load libraries
library(stars)
library(dplyr)
library(gtools)

# source functions
source("Functions/process_posterior_run69.R")
source("Functions/select_colony_year.R")

################################################################################

# Read in data

################################################################################

# set og and des crs
og_crs<-4326
des_crs<-3035

# CHOOSE YEAR (1 = year 1900, i.e. 112= year 2011)
year<-112

# read in colonies data
colonies<-read.csv("Data/posterior_run69.csv")

# rename columns e.t.c for use in other functions
colonies<-process_posterior_run69(colonies)

# British Isles colonies
ukire<-c("Ailsa Craig", "Barra Head", "Bass Rock", "Bempton Cliff", "Copinsay",
         "Fair Isle", "Flannan Isles", "Foula", "Grassholm","Great Saltee",
         "Hermaness", "Holy Isle Arran","Isle of May","Les Etacs", "Lundy", "Noss", "Ortac", "Rockall", "Scar Rocks",
         "Shiant Islands", "St Kilda", "St. Margarets", "Sula Sgeir", "Sule Skerry",
         "Sule Stack", "Westray", "Troup Head", "Bull Rock", "Clare Island",
         "Great Saltee", "Ireland's Eye", "Lambay", "Little Skellig", "Rouzic")

# select colonies
colonies_uk<- dplyr::filter(colonies, Colony %in% ukire)

# select specific colony year
colonies.year<-select_colony_year(colonies=colonies_uk, year=year)

# order colonies by ID
colonies.year<-colonies.year%>%arrange(ID)

# extract IDs for all colonies
all_col_ID<-colonies.year$ID

# read in model HRs and sort
HR_2011_files<-list.files("Data/",
                          pattern="*home_range", full.names=TRUE)

# order HR files and read in
HR_2011<-read_stars(mixedsort(HR_2011_files))

# read in projected distributions
projected_files<-list.files("Data/", full.names=T)
projected_files_sort<-mixedsort(projected_files)
proj_dist<-read_stars(projected_files_sort)
st_crs(proj_dist)<-3035

# read in foraging radii
radii_files<-list.files("Data/", full.names=T)
radii_files_sort<-mixedsort(radii_files)
radius<-read_stars(radii_files_sort)


# read in tracking data
les_etacs<-read.csv("Data/")
les_etacs.sf<-st_as_sf(les_etacs, coords=c("longitude", "latitude"), crs=og_crs)
les_etacs.sf<-st_transform(les_etacs.sf, des_crs)

sule_skerry<-read.csv("Data/")
sule_skerry.sf<-st_as_sf(sule_skerry, coords=c("longitude", "latitude"), crs=og_crs)
sule_skerry.sf<-st_transform(sule_skerry.sf, des_crs)

bull_rock<-read.csv("Data/")
bull_rock.sf<-st_as_sf(bull_rock, coords=c("longitude", "latitude"), crs=og_crs)
bull_rock.sf<-st_transform(bull_rock.sf, des_crs)

bass_rock<-read.csv("Data/")
bass_rock.sf<-st_as_sf(bass_rock, coords=c("longitude", "latitude"), crs=og_crs)
bass_rock.sf<-st_transform(bass_rock.sf, des_crs)

ailsa_craig<-read.csv("Data/")
ailsa_craig.sf<-st_as_sf(ailsa_craig, coords=c("longitude", "latitude"), crs=og_crs)
ailsa_craig.sf<-st_transform(ailsa_craig.sf, des_crs)

lambay<-read.csv("Data/")
lambay.sf<-st_as_sf(lambay, coords=c("longitude", "latitude"), crs=og_crs)
lambay.sf<-st_transform(lambay.sf, des_crs)

ile_rouzic<-read.csv("Data/")
ile_rouzic.sf<-st_as_sf(ile_rouzic, coords=c("longitude", "latitude"), crs=og_crs)
ile_rouzic.sf<-st_transform(ile_rouzic.sf, des_crs)

little_skellig<-read.csv("Data/")
little_skellig.sf<-st_as_sf(little_skellig, coords=c("longitude", "latitude"), crs=og_crs)
little_skellig.sf<-st_transform(little_skellig.sf, des_crs)


################################################################################

# Create data structures and assign tracking locations

################################################################################

# create 26 by 26 matrix to run assignments to colonies
prob_assignment.mat<-matrix(data=0, nrow=dim(colonies.year)[1], ncol=dim(colonies.year)[1])
prob_assignment_radius.mat<-matrix(data=0, nrow=dim(colonies.year)[1], ncol=dim(colonies.year)[1])
prob_assignment_proj.mat<-matrix(data=0, nrow=dim(colonies.year)[1], ncol=dim(colonies.year)[1])

select_col_ID<-c(1,3,8,16,17,27,28,29,36,48)

# create data frame recording colony tracking data points with colony IDs
confusion_pop.df<-data.frame(ID=c(rep(select_col_ID[1], dim(ailsa_craig.sf)[1]), rep(select_col_ID[2], dim(bass_rock.sf)[1]),
                                  rep(select_col_ID[3], dim(bull_rock.sf)[1]), rep(select_col_ID[6], dim(lambay.sf)[1]),
                                  rep(select_col_ID[7], dim(les_etacs.sf)[1]),
                                  rep(select_col_ID[8], dim(little_skellig.sf)[1]),
                                  rep(select_col_ID[9], dim(ile_rouzic.sf)[1]), rep(select_col_ID[10], dim(sule_skerry.sf)[1])),
                             assigned_col=rep(NA, 298637))

confusion_proj.df<-data.frame(ID=c(rep(select_col_ID[1], dim(ailsa_craig.sf)[1]), rep(select_col_ID[2], dim(bass_rock.sf)[1]),
                                   rep(select_col_ID[3], dim(bull_rock.sf)[1]), rep(select_col_ID[6], dim(lambay.sf)[1]),
                                   rep(select_col_ID[7], dim(les_etacs.sf)[1]),
                                   rep(select_col_ID[8], dim(little_skellig.sf)[1]),
                                   rep(select_col_ID[9], dim(ile_rouzic.sf)[1]), rep(select_col_ID[10], dim(sule_skerry.sf)[1])),
                              assigned_col=rep(NA, 298637))


confusion_radius.df<-data.frame(ID=c(rep(select_col_ID[1], dim(ailsa_craig.sf)[1]), rep(select_col_ID[2], dim(bass_rock.sf)[1]),
                                     rep(select_col_ID[3], dim(bull_rock.sf)[1]), rep(select_col_ID[6], dim(lambay.sf)[1]),
                                     rep(select_col_ID[7], dim(les_etacs.sf)[1]),
                                     rep(select_col_ID[8], dim(little_skellig.sf)[1]),
                                     rep(select_col_ID[9], dim(ile_rouzic.sf)[1]), rep(select_col_ID[10], dim(sule_skerry.sf)[1])),
                                assigned_col=rep(NA, 298637))


# combine tracking locations into one object
combined_tracking<-rbind(ailsa_craig.sf[,18],bass_rock.sf[,18], bull_rock.sf[,18],
                         lambay.sf[,20], les_etacs.sf[,18], little_skellig.sf[,20],ile_rouzic.sf[,20], sule_skerry.sf[,18])


# record number of tracking locations for each colony
tracking_length<-c(dim(ailsa_craig.sf)[1],dim(bass_rock.sf)[1],dim(bull_rock.sf)[1],
                   dim(lambay.sf)[1], dim(les_etacs.sf)[1],
                   dim(little_skellig.sf)[1],dim(ile_rouzic.sf)[1],dim(sule_skerry.sf)[1])

# multiply 2011 HRs by colony size in 2011
HR_2011_pop<-HR_2011*0
for(i in 1:length(all_col_ID)){

  HR_2011_pop[[i]]<-HR_2011[[i]]*colonies.year$median[i]
}

# extract HR cell values at each tracking location
HR_point_values_pop<-st_extract(HR_2011_pop, combined_tracking)

HR_point_values_pop.df<-st_drop_geometry(HR_point_values_pop)

HR_point_values_proj<-st_extract(proj_dist, combined_tracking)

HR_point_values_proj.df<-st_drop_geometry(HR_point_values_proj)

HR_point_values_radius<-st_extract(radius, combined_tracking)

HR_point_values_radius.df<-st_drop_geometry(HR_point_values_radius)

# number of repetitions of assignment
nr<-4

# for each repetition, apportion tracking locations probabilistically
# from the proportion of usage of each colony in each cell
for (j in 1:nr){

  # loop through each tracking location
  for(i in 1:dim(combined_tracking)[1]){

    # check values in cell for each method
    HR_cell<-HR_point_values_pop.df[i,]
    HR_cell_radius<-HR_point_values_radius.df[i,]
    HR_cell_proj<-HR_point_values_proj.df[i,]

    # if there are any NA values, set these to 0
    HR_cell[is.na(HR_cell)]<-0
    HR_cell_radius[is.na(HR_cell_radius)]<-0
    HR_cell_proj[is.na(HR_cell_proj)]<-0

    # if all values are zero, assign every colony value of 1
    # (i.e. probability of assignment is the same for each colony)
    if(all(HR_cell==0)){HR_cell<-rep(1, dim(colonies.year)[1])}
    if(all(HR_cell_radius==0)){HR_cell_radius<-rep(1, dim(colonies.year)[1])}
    if(all(HR_cell_proj==0)){HR_cell_proj<-rep(1, dim(colonies.year)[1])}

    # name vectors of cell values
    names(HR_cell)<-seq_along(HR_cell)
    names(HR_cell_radius)<-seq_along(HR_cell_radius)
    names(HR_cell_proj)<-seq_along(HR_cell_proj)

    # calculate probabilities from row of values
    HR_prob<-HR_cell/sum(HR_cell, na.rm=T)
    HR_prob_radius<-HR_cell_radius/sum(HR_cell_radius, na.rm=T)
    HR_prob_proj<-HR_cell_proj/sum(HR_cell_proj, na.rm=T)

    # assignment to colonies based on probabilities
    sample_prob<-sample(HR_cell, size=1, prob=HR_prob)
    sample_prob_radius<-sample(HR_cell_radius, size=1, prob=HR_prob_radius)
    sample_prob_proj<-sample(HR_cell_proj, size=1, prob=HR_prob_proj)

    # index of value associated with selected colony
    sample_index<-which(names(HR_cell) %in% names(sample_prob))
    sample_index_radius<-which(names(HR_cell_radius) %in% names(sample_prob_radius))
    sample_index_proj<-which(names(HR_cell_proj) %in% names(sample_prob_proj))

    # find position of ID in all colony IDs
    ID_pos<-which(all_col_ID %in% confusion_pop.df[i,1])
    ID_pos_radius<-which(all_col_ID %in% confusion_radius.df[i,1])
    ID_pos_proj<-which(all_col_ID %in% confusion_proj.df[i,1])

    # assign point to colony in matrix i.e. add 1 to cell
    prob_assignment.mat[ID_pos,sample_index]<-(prob_assignment.mat[ID_pos,sample_index]+1)
    prob_assignment_radius.mat[ID_pos_radius,sample_index_radius]<-(prob_assignment_radius.mat[ID_pos_radius,sample_index_radius]+1)
    prob_assignment_proj.mat[ID_pos_proj,sample_index_proj]<-(prob_assignment_proj.mat[ID_pos_proj,sample_index_proj]+1)

  }

}

# calculate no. tracking locations correctly assigned
sum(diag(prob_assignment.mat))/(sum(prob_assignment.mat))
sum(diag(prob_assignment_radius.mat))/sum(prob_assignment_radius.mat)
sum(diag(prob_assignment_proj.mat))/sum(prob_assignment_proj.mat)


