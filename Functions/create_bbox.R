#' create_bbox
#'
#' A function to create a bounding box encompassing all colonies' max foraging ranges.
#'
#' @param colonies_sf \code{sf} object created from 'create_colonies_sf' function containing colony locations, populations and
#' colony IDs
#' @param max_forage \code{numeric} value of maximum foraging range for colonial animal in unit of CRS of colonies_sf
#'
#' @return Returns a \code{bbox} object
#' @export
#'
#' @examples
#'
create_bbox<-function(colonies_sf, max_forage){

  # check colonies object of correct class
  if (class(colonies_sf)[1]!="sf" & class(colonies_sf)[2]!="data.frame") stop(paste("Error: colonies data inputted is
                                           not of sf data type. Data type of current colonies is",
                                           class(colonies.sf),", please see create_colonies_sf function for more info.", sep=""))


  # create buffers to decide on bounding box
  col_buffer<-sf::st_buffer(colonies_sf, max_forage)


  # extract bbox
  bbox<-sf::st_bbox(col_buffer)

  # return bbox
  return(bbox)
}
