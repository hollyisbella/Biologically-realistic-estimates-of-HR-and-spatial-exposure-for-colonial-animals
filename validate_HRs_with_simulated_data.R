################################################################################

# Validation with simulated data
# Test similarity of true and estimated HRs

################################################################################

# README #

# wind farm polygon data can be openly accessed from EMODnet Human Activities at (accessed on 2023-11-24):
# https://ows.emodnet-humanactivities.eu/geonetwork/srv/eng/catalog.search#/metadata/8201070b-4b0b-4d54-8910-abcea5dce57f

################################################################################

# load libraries
library(stars)
library(sf)
library(terra)
library(amt)
library(ggplot2)
library(ggspatial)

################################################################################

# Load surfaces and plotting info

################################################################################

# read in HRs
true_location<-""

est_location<-""

# true HRs
true_ud_list<-list.files("Data/",
                        pattern="*exposure_surface", full.names=TRUE)

# estimated HRs
est_ud_list<-list.files("Data/",
                        pattern="*exposure_surface", full.names=TRUE)

res<-5000

################################################################################

# HR similarity comparison

################################################################################

# read in true and estimated surfaces as stars objects
true_ud.stars<-read_stars(true_ud_list)
est_ud.stars<-read_stars(est_ud_list)

true_ud.stars[is.na(true_ud.stars)]<-0
est_ud.stars[is.na(est_ud.stars)]<-0

# create matrices to store values of bhattacharyya's coefficient (affinity) +
# probability of finding estimated HR within true HR (PHR)
bc<-matrix(nrow=1,ncol=length(est_ud.stars))
phr_true<-matrix(nrow=1, ncol=length(est_ud.stars))
phr_est<-matrix(nrow=1, ncol=length(est_ud.stars))

# calculate BA and PHR
for (i in 1:length(est_ud.stars)){
  bc[1,i]<-sum(sqrt(c(est_ud.stars[[i]]))*sqrt(c(true_ud.stars[[i]])))
  phr_true[1,i]<-sum(true_ud.stars[[i]][est_ud.stars[[i]]>0])
  phr_est[1,i]<-sum(est_ud.stars[[i]][true_ud.stars[[i]]>0])
}

# calculate quantiles and medians
quantile(phr_true,0.05)
quantile(phr_true,0.95)
median(phr_true)
quantile(bc,0.05)
quantile(bc,0.95)
median(bc)

# create df of similarity info
overlap.df<-data.frame(ba=c(bc), phr=c(phr_true))

# write overlap to file
write.csv(overlap.df,file="", row.names=F)

################################################################################

# Compare exposure calculated between true and estimated HRs

################################################################################

# select year for 'present day'
present<-2011

# read in windfarm polygons
eu<-read_sf("Data/")

# check status's
unique(eu$STATUS)

# make all polygons valid
eu<-st_make_valid(eu)

# transform to desired crs
eu<-st_transform(eu,crs=3035)

# filter countries for windfarms surrounding British Isles
eu<-eu %>%filter(COUNTRY %in% c("United Kingdom","Ireland","France","Belgium","Netherlands"))

# filter for planned, existing or in-construction windfarms
eu_planned<-filter(eu, eu$STATUS %in% c("Planned", "Approved") |
                     (eu$STATUS %in% c("Production","Test site") & eu$YEAR>present))

# rasterize planned windfarms
eu_all.stars<-st_rasterize(eu_planned, dx=res,dy=res,template=st_as_stars(st_bbox(est_ud.stars),
                                                                      dx=res, dy=res))

# set all cells with windfarm presence to 1
all_wind<-eu_all.stars[1,,]
all_wind[all_wind>=1]<-1

# create objects to store exposure to windfarms
wind_true_all<-est_ud.stars*0
wind_est_all<-est_ud.stars*0

sum_effect_true_all<-rep(0,(length(est_ud.stars)-1))
sum_effect_est_all<-rep(0,(length(est_ud.stars)-1))

for (k in 1:(length(est_ud.stars)-1)){

  # overlap of windfarms with HRs
  wind_true_all[[k]]<-true_ud.stars[[k]]*all_wind[[1]]
  wind_est_all[[k]]<-est_ud.stars[[k]]*all_wind[[1]]

  # sum of overlap with HRs
  sum_effect_true_all[k]<-sum(wind_true_all[[k]], na.rm=T)
  sum_effect_est_all[k]<-sum(wind_est_all[[k]], na.rm=T)

}

# plot against each other
plot(sum_effect_true_all, sum_effect_est_all)

# fit linear model
fit<-lm(sum_effect_true_all~sum_effect_est_all)

# summary output
fit_summary<-summary(fit)

# calculate CIs
ci<-confint(fit)
ci

# create .df
compare_exp<-data.frame(coefficients=fit$coefficients, ci=ci)

names(compare_exp)[2:3]<-c("2.5%","97.5%")
rownames(compare_exp)<-c("intercept","slope")

# save to file
write.csv(compare_exp, file="")

