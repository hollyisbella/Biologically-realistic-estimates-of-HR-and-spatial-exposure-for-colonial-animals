#' create_colonies_stars_array
#'
#' A function to create a colonies stars raster of specified km resolution, with
#' each colony size at colony centre in a separate layer.
#' Using landseamask (created from 'create_landseamask') as template
#'
#'  @param colonies_sf \code{sf} object created from 'create_colonies_sf' function containing colony locations, sizes and
#' IDs
#' @param landseamask \code{stars} object as template - created from 'create_landseamask'
#' @param res \code{numeric} resolution of raster in units of crs being used (des_crs)
#'
#' @return Returns a list containing \code{stars} raster object with layers for each colony by ID, colony sizes located at colony centre cell and \code{dataframe} containing cell coordinates of colonies
#' @export
#'
#' @examples
#'

create_colonies_stars_array<-function(colonies_sf, landseamask, res){

  # extract colony IDs
  pop_id<-colonies_sf$ID

  # create names of layers
  pop_names<-paste("a", pop_id, sep="")

  # select colony sizes
  pop<-colonies_sf %>% dplyr::select(median)

  # loop through all colonies
  for (i in 1:length(pop$median)){

    # for the first colony
    if (i==1){

      # select ith colony layer with info on colony size
      pop_temp<-pop[i,"median"]

      # rasterize colony sf object using landseamask template
      a<-stars::st_rasterize(pop_temp,template=stars::st_as_stars(st_bbox(landseamask),
                                                           dx=res, dy=res))

      # get cell coordinates of colony centre
      central_loc_x<-unname(which(a$median>0, arr.ind=T)[,1])
      central_loc_y<-unname(which(a$median>0, arr.ind=T)[,2])

      # create data frame to store cell coordinates
      col_pos<-data.frame(x=central_loc_x, y=central_loc_y)

      # for all other colonies
    } else{

      pop_temp<-pop[i,"median"]

      # create colony layer (temporary object)
      a_temp<-stars::st_rasterize(pop_temp,template=stars::st_as_stars(st_bbox(landseamask),
                                                                dx=res, dy=res))

      central_loc_x<-unname(which(a_temp$median>0, arr.ind=T)[,1])
      central_loc_y<-unname(which(a_temp$median>0, arr.ind=T)[,2])

      # append colony central cell coordinates to data frame
      col_pos<-rbind(col_pos, data.frame(x=central_loc_x, y=central_loc_y))

      # assign names to each colony
      assign(pop_names[i],a_temp)

      # append colony layer to array of colonies
      a<-c(a,a_temp)

      # remove temporary colony layer
      rm(a_temp)
    }
  }

  names(a)<-pop_names

  list_pop_pos<-list(a, col_pos)

  # return colony cell coordinates and array
  return(list_pop_pos)
}
