#'move_colonies2
#'
#' A function to check if colony centers are on land and if they are,
#' move to a permitted cell
#'
#' @param col_pos \code{data.frame} containing colony locations (by grid cell) created from "create_colonies_stars_array"
#' @param landseamask \code{stars} object with information on land and sea cells - created from 'create_landseamask'
#' @param a \code{stars} object with colony sizes - created from 'create_colonies_stars' function
#' @return Returns a \code{stars} object with colony centres in new locations and \code{data.frame} with new positions of colony centres
#' @export
#'
#' @examples
#'

move_colonies2<-function(landseamask, col_pos, a){

  # define neighbourhood
  nei<-cbind(c(1,0,-1,0,-1,1,-1,1), c(0,1,0,-1,1,1,-1,-1))

  # loop through colonies
  for (k in 1:length(a)){

    # define central coordinates of colony
    x_loc<-col_pos[k,1]
    y_loc<-col_pos[k,2]

    # define neighbours of cell centres
    nei_xy<-cbind(x_loc+nei[,1],y_loc+nei[,2])

    # check if colony on land
    check<-landseamask$mask[x_loc,y_loc]==-999

    # if colony located on a land cell
    if (check==1){

      # check if any surrounding cells are sea cells
      check_nei<-landseamask$mask[nei_xy]==0

      # if no surrounding cells are sea - expand nei definition
      if (sum(check_nei)==0){

        # new definition of neighbourhood
        nei<-cbind(rep(c(-2,-1,0,1,2),5),c(rep(2,5),rep(1,5),rep(0,5),rep(-1,5),rep(-2,5)))

        # expanded neighbourhood from specified colony centre
        nei_xy<-cbind(x_loc+nei[,1],y_loc+nei[,2])

        # check if any neighbours are sea cells
        check_nei<-landseamask$mask[nei_xy]==0

      }

      # if one available sea cell in neighbourhood
      # move colony to this cell
      if (sum(check_nei)==1){

        # relocate colony to permitted (non-land) location
        a[[k]][nei_xy]<-a[[k]][x_loc,y_loc]*check_nei

        # set old position to 0
        a[[k]][x_loc,y_loc]<-0

        # if more than one sea cell available
        # choose first cell in list
      } else if (sum(check_nei)>1){

        # loop through neighbourhood checked values
        for (l in 1:length(check_nei)){

          # check if neighbourhood cell is a sea cell
          check_check_nei<-check_nei[l]==1

          # if cell is a sea cell set all other neighbourhood cells to 0
          # (i.e select first sea cell in neighbourhood)
          if(check_check_nei==1){
            check_nei<-rep(0,length(check_nei))
            check_nei[l]<-1
          }

        }

        # reassigning colony centre to permitted cell
        a[[k]][nei_xy]<-a[[k]][x_loc,y_loc]*check_nei

        # remove old centre
        a[[k]][x_loc,y_loc]<-0


      }

        # update colony center position in grid
        col_pos[k,1]<-unname(which(a[[k]]>0, arr.ind=T)[,1])
        col_pos[k,2]<-unname(which(a[[k]]>0, arr.ind=T)[,2])

    }

  }
  return_list<-list(a, col_pos)

  return(return_list)
}
