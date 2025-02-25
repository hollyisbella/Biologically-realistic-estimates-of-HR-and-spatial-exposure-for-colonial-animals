#' create_p_layer
#'
#' A function to create a permeability layer from a given mask specifying impermeable areas and permeable areas
#'
#' @param mask \code{stars} object landseamask - created from 'create_landseamask'
#' @param cap \code{num} giving carrying capacity value over space (homogeneous environment)
#' @return Returns a \code{stars} raster object containing permeability layer p
#' @export
#'
#' @examples
#'
create_p_layer<-function(mask, cap){

  # create p to be same size as mask
  p<-0*mask
  # name permeability layer
  names(p)<-"perm"

  # set land to be impermeable (infinitesimally small carrying capacity)
  p$perm[mask$mask==-999]<-cap*(1/100000) # set land to be impermeable

  p$perm[p$perm==0]<-cap # set sea to be homogeneously permeable (same carrying capacity)

  return(p)
}

