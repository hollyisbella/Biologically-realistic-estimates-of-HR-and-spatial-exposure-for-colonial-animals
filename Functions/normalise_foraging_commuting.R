#' normalise_foraging_commuting
#'
#' A function to create a bounding box encompassing all colonies maximal foraging ranges.
#'
#' @param a \code{stars} foraging surface for each colony - output from 'calculate_foraging'
#' @param commuting \code{stars} commuting surface for each colony - output 'calculate_commute2'
#'
#' @return Returns normalised foraging and commuting distributions as \code{stars} raster objects
#' @export
#'
#' @examples
#'


normalise_foraging_commuting<-function(a, commuting){

  # normalise surfaces
  # dividing by sum of flux/usage
  for (k in 1:length(a)){
    suma<-sum(a[[k]], na.rm=T)
    sumcom<-sum(commuting[[k]], na.rm=T)

    a[[k]]<-a[[k]]/suma
    commuting[[k]]<-commuting[[k]]/sumcom
  }

  # return normalised surfaces
  comb<-list(a,commuting)
  return(comb)
}
