# FUNCTION to run before model fitting
pre_optim<-function(land_poly, colonies, year, og_crs, des_crs, max_forage, res){

  # select specific colony year
  colonies.year<-select_colony_year(colonies=colonies, year=year)

  # create colonies .shp file
  colonies.sf<-create_colonies_sf(colonies.year, og_crs, des_crs)

  # create bbox
  bbox<-create_bbox(colonies_sf=colonies.sf, max_forage=forage)

  # create landseamask from polygons shapefile
  land_rast<-create_landseamask(bbox=bbox, des_crs=des_crs, res=res, land_poly=land_poly )

  # create colonies raster and get colony center cell locations
  list<-create_colonies_stars_array(colonies_sf=colonies.sf, landseamask=land_rast, res=res)

  a<-list[[1]]
  col_pos<-list[[2]]

  # move colonies off land
  updated_positions<-move_colonies2(landseamask=land_rast, a=a,col_pos=col_pos)

  # extract updated colony positions
  a_update<-updated_positions[[1]]
  col_pos_update<-updated_positions[[2]]

  return_list<-list(a_update, land_rast, col_pos_update, colonies.sf)

  return(return_list)
}
