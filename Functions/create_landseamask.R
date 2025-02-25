#' create_landseamask
#'
#' A function to create a land sea mask of specified resolution.
#'
#' @param bbox \code{bbox} bounding box (as created by create_bbox function) of study area
#' @param des_crs \code{numeric} as used in bounding box
#' @param res \code{numeric} resolution of raster in units of des_crs and bbox
#' @param land_poly \code{MULTIPOLYGON} sf object containing land polygons
#'
#' @return Returns a \code{stars} raster object specifying land (-999) and sea (0)
#' @export
#'
#' @examples
#'

create_landseamask<-function(bbox, des_crs, res, land_poly){

  # check bbox object of correct class
  if (class(bbox)!="bbox") stop(paste("Error: bbox inputted is not of bbox data type. Data type of current bbox is", class(bbox), sep=""))

  # check crs of land polygon against bbox, convert crs of land polygon as needed
  if (sf::st_crs(land_poly)!=sf::st_crs(bbox)) land_poly<-sf::st_transform(land_poly, des_crs)

  # select feature layer of land-sea mask
  land_poly<-land_poly["featurecla"]

  # create raster of land polygon based on bbox
  land_rast<-stars::st_rasterize(land_poly, st_as_stars(bbox, dx=res, dy=res, values=NA))

  # rename layer
  names(land_rast)<-"mask"

  # all values not equal to NA are land
  land_rast[!is.na(land_rast)]<--999

  # all NA values are sea
  land_rast[is.na(land_rast)]<-0

  return(land_rast)
}
