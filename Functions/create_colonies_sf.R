#' create_colonies_sf
#'
#' A function to create a Shapefile of colonies
#'
#' @param colonies_year \code{data.frame} containing colony locations, sizes and
#' colony IDs with columns names 'latitude', 'longitude', 'median', 'ID' for 1 year
#' @param og_crs \code{numeric} original crs of inputted colony location data (e.g 4326)
#' @param des_crs \code{numeric} desired crs for data
#'
#' @return Returns a \code{sf} object with colonies data
#' @export
#'
#' @examples
#'

create_colonies_sf<-function(colonies_year, og_crs, des_crs){

  # check for existence of columns for latitude, longitude, ID and pop
  if(sum(names(colonies_year) %in% c("latitude", "longitude","ID","median"))!=4) print("Check column
                                                                    names of inputted colonies dataframe.
                                                                    Colony locations should be specified
                                                                    under the column names 'latitude' and 'longitude',
                                                                    colony populations as 'median' and colony IDs as 'ID'.")

  # extract columns of IDs, colony sizes and locations (longitude, latitude)
  colonies.df<-data.frame(ID=colonies_year$ID, median=round(colonies_year$median), longitude=colonies_year$longitude, latitude=colonies_year$latitude)

  # convert df to sf object specifying original crs
  colonies.sf<-sf::st_as_sf(colonies.df, coords=c("latitude","longitude"), crs=og_crs)

  # transform crs to desired crs
  colonies.sft<-sf::st_transform(colonies.sf, des_crs)

  return(colonies.sft)
}
