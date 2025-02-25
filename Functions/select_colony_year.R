#' select_colony_year
#'
#' A function to filter a dataframe of colony sizes and locations for a specific colony year
#'
#' @param colonies \code{data.frame} containing colony locations, populations and
#' colony IDs with columns names 'latitude', 'longitude', 'median' and 'ID' and 'parameters'
#'
#' @param year \code{numeric} specify year to be selected (i.e year search term)
#'
#' @return Returns a \code{data.frame} object of colony year-specific information
#' @export
#'
#' @examples
#'

select_colony_year<-function(colonies, year){

  # specify colony IDs
  colony_ID<-unique(colonies$ID)

  # names of parameters to select
  name_params<-paste("P", colony_ID, ".",year, sep="" )

  # select single year of colony populations
  coloniesyear<-dplyr::filter(colonies, parameter %in% name_params)

  return(coloniesyear)
}
