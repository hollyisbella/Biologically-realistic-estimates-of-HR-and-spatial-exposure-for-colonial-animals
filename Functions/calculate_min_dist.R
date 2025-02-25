#' calculate_min_dist
#'
#' A function to calculate shortest distances from colonies to all permeable cells
#'
#' @param mask \code{stars} object landseamask - created from 'create_landseamask'
#' @param a \code{stars} object with colony populations (one colony per layer)- created from 'create_colonies_stars_array' function
#' @param col_pos \code{data frame} containing cell locations of colony centres
#' @param avoid_land \code{logical} should impermeable cells (land) be avoided in distance calculation
#'
#' @return Returns a \code{stars} raster object containing min distances from colony centres to all permeable cells
#' @export
#'
#' @examples
#'

calculate_min_dist<-function(a, mask, col_pos, avoid_land){

  for (k in 1:length(a)){

    # select colony centre location
    xpos<-col_pos$x[k]
    ypos<-col_pos$y[k]

    ak<-a[k,,] # extract kth colony layer

    ak[[1]][ak[[1]]>=0]<-1 # define cells to calculate distance to
    ak[[1]][xpos, ypos]<-2 # define colony centre (from which to calculate distance)

    # if calculating distance accounting for impermeable areas
    if (avoid_land==T){
      ak[[1]][mask$mask==-999] <-NA #set impermeable areas to NA
    }

    # convert kth layer to terra SpatRaster
    aterrak<-as(ak, "SpatRaster")

    # if first layer
    if (k==1){
      # calculate distance in km
      dist_rast<-terra::gridDist(aterrak, target=2, scale=1000)

      # convert to stars object
      dist<-stars::st_as_stars(dist_rast)

    } else{

      dist_temp_rast<-terra::gridDist(aterrak, target=2, scale=1000)
      dist_temp<-stars::st_as_stars(dist_temp_rast)

      # append distance for colony k to distance raster
      dist<-c(dist,dist_temp)

      # remove temporary objects
      rm(dist_temp_rast)
      rm(dist_temp)
    }

  }

  return(dist)
}
