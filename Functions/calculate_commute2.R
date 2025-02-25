#' calculate_commute2
#'
#' A function to calculate commuting surfaces from foraging surfaces (from 'calculate_foraging')
#'
#' @param dist \code{stars} object which holds distances from colony centres
#' @param a \code{stars} final foraging surface (output from 'calculate_foraging')
#' @param mask \code{stars} landseamask for plotting
#' @param saveloc \code{string} location to save simulation output images
#' @param save_freq \code{num} no. iterations between saving each image
#' @param tmax \code{num} giving maximum iteration value
#' @param colonies_sf \code{sf} object created from 'create_colonies_sf' function containing colony locations, populations and
#' colony IDs
#' @param thresh_com \code{num} threshold at which to stop running simulation (between 0 and 1)
#' @param plot_results \code{logic} value of T (TRUE) to plot, value of F (FALSE) not to plot
#' @param thresh_plot_a \code{num} value at which to threshold individuals being plotted
#'
#' @return Returns a \code{stars} raster object containing commuting surface
#' @export
#'
#' @examples
#'
calculate_commute2<-function(a, dist, mask, col_pos, save_loc, save_freq, tmax,
                             colonies_sf, thresh_com, plot_results, thresh_plot_a){

  # get size of grid
  x_size<-unname(dim(a))[1]
  y_size<-unname(dim(a))[2]

  # define neighbourhood
  nei<-cbind(c(1,0,-1,0,-1,1,-1,1), c(0,1,0,-1,1,1,-1,-1))

  # make blank raster to store commuting surface
  commute<-a*0

  # create vector to check stopping thresholds
  sum_diff_a<-rep(1,length(a))

  for (t in 1:tmax){
    # create temporary foraging surface
    ad<-a*0

    # loop through all cells in grid
    for (i in 2:(x_size-1)){
      for (j in 2:(y_size-1)){

        aij<-c(rep(0,length(a)))
        distij<-c(rep(0,length(a)))

        # for each colony
        for (k in 1:length(a)){
          # extract foraging usage in cell ij
          aij[k]<-a[[k]][i,j]
          # extract distance of cell ij to kth colony
          distij[k]<-dist[[k]][i,j]
        }

        # determine number of colonies in current cell
        ncol<-sum(aij>0, na.rm=T)

        # determine position of colonies present in cell in list
        col_no<-which(aij>0)

        # extract distances to all colonies from cell ij
        distij<-dist[,i,j]

        # if only one colony in cell
        if (ncol==1){

          # and usage not already in centre cell
          if (distij[[col_no]]!=0){

            # extract distances of neighbouring cells
            check_nei<-dist[[col_no]][cbind(i+nei[,1],j+nei[,2])]

            # check which of surrounding cells have a lower distance to colony
            # centre than current cell
            dec_dist<-rep(distij[[col_no]], length(check_nei))>check_nei

            # check no NAs in distance check
            if(any(is.na(dec_dist))) stop("NA in dec_dist")

            # split usage between number of cells to move to
            dn<-aij[[col_no]]/sum(dec_dist, na.rm=T)

            # move usage to neighbouring cells with lower distances to colony centers
            ad[[col_no]][cbind(i+nei[,1],j+nei[,2])]<-ad[[col_no]][cbind(i+nei[,1],j+nei[,2])]+dn*dec_dist

            # summing usage in current cell (flux)
            commute[[col_no]][i,j]<-commute[[col_no]][i,j]+aij[col_no]

            # move usage out of current cell
            ad[[col_no]][i,j]<-ad[[col_no]][i,j]

            # checks for NAs
            if(any(is.na(ad[[col_no]][cbind(i+nei[,1],j+nei[,2])]))) stop("NA in adnei")
            if(is.na(ad[[col_no]][i,j])) stop("NA in ad")

          } else{

            # if usage already in colony centre cell, leave in this cell
            ad[[col_no]][i,j]<-ad[[col_no]][i,j]+a[[col_no]][i,j]

          }

          # if more than one colony in cell
        } else if (ncol>1){

          # loop through all colonies
          for (ni in 1:length(col_no)){

            # colony position in list
            col<-col_no[ni]

            # if colony usage not at centre cell
            if (distij[[col]]!=0){
              check_nei<-dist[[col]][cbind(i+nei[,1],j+nei[,2])]  # check neighbouring cell distances
              dec_dist<-rep(distij[[col]], length(check_nei))>check_nei

              if(any(is.na(dec_dist))) stop("NA in dec_dist")

              dn<-aij[[col]]/sum(dec_dist, na.rm=T)

              ad[[col]][cbind(i+nei[,1],j+nei[,2])]<-ad[[col]][cbind(i+nei[,1],j+nei[,2])]+dn*dec_dist

              commute[[col]][i,j]<-commute[[col]][i,j]+aij[col]

              ad[[col]][i,j]<-ad[[col]][i,j]

              if(any(is.na(ad[[col]][cbind(i+nei[,1],j+nei[,2])]))) stop("NA in adnei")
              if(is.na(ad[[col]][i,j])) stop("NA in ad")

            } else{

              # if already in colony centre cell retain usage here
              ad[[col]][i,j]<-ad[[col]][i,j]+a[[col]][i,j]

            }
          }
        }

      }
    }

    # update reversing foraging distribution
    a<-ad

    if (t==1){
      astore<-a
    } else{

      for (k in 1:length(a)){

        # absolute residuals between current and previous foraging surface
        # which is being pulled back to the colony centre
        diff_a<-abs(a[[k]]-astore[[k]])/colonies_sf$median[k]

        # sum residuals
        sum_diff_a[k]<-sum(diff_a)
      }

    }

    # if within threshold for all colonies stop running calculation
    if(all(sum_diff_a<thresh_com)) break

    astore<-a # store t-1 reverse foraging surface for comparison in next iteration


    # if true plot
    if((t/save_freq==round(t/save_freq)| t==1) & plot_results==T) {
      aplot<-a
      complot<-commute
      suma<-rep(0,length(a))
      sumcom<-rep(0,length(a))
      for (k in 1:length(a)){

        # normalise reverse foraging surface and commuting surface
        suma[k]<-sum(a[[k]])
        sumcom[k]<-sum(commute[[k]])

        aplot[[k]][aplot[[k]]<=thresh_plot_a]<-NA
        aplot[[k]]<-aplot[[k]]/suma[k]

        complot[[k]][complot[[k]]==0]<-NA
        complot[[k]]<-complot[[k]]/sumcom[k]
      }

      # plot reverse foraging surface
      tmap::tmap_mode("plot")
      tmapsavea<-tmap::tm_shape(mask) + tmap::tm_raster(legend.show = FALSE, col="mask", palette=c("lightblue", "darkgreen")) +
        tmap::tm_shape(aplot) + tmap::tm_raster(style="cont", breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+tmap::tm_facets(as.layers=T)+
        tmap::tm_layout(title="Reverse Foraging Distribution", title.position = c("center","top"))

      # plot commuting surface
      tmapsavecom<-tmap::tm_shape(mask) + tmap::tm_raster(legend.show = FALSE, col="mask", palette=c("lightblue", "darkgreen")) +
        tmap::tm_shape(complot) + tmap::tm_raster(style="cont", breaks=c(0,0.005,0.01,0.015,0.02,0.025,0.03,0.035))+tmap::tm_facets(as.layers=T)+
        tmap::tm_layout(title="Commuting Distribution", title.position = c("center","top"))

      # save plots to specified file locations
      tmap::tmap_save(tmapsavea, paste(save_loc,"reverse_forage_dist_t_", t,".png", sep=""))
      tmap::tmap_save(tmapsavecom, paste(save_loc,"commute_dist_t_", t,".png", sep=""))

      # remove plot objects
      rm(complot, aplot)
    }



  }

  # set commute of colony centre to total colony usage
  for (k in 1:length(a)){
    xpos<-col_pos$x[k]
    ypos<-col_pos$y[k]
    commute[[k]][xpos,ypos]<-a[[k]][xpos,ypos]
  }

  # final plot
  if (plot_results==T){

    sumcom<-rep(0,length(a))
    for (k in 1:length(a)){
      sumcom[k]<-sum(commute[[k]])
      commute[[k]][commute[[k]]<=thresh_plot_a]<-NA
      commute[[k]]<-commute[[k]]/sumcom[k]
    }

    tmap::tmap_mode("plot")

    tmapsavecom<-tmap::tm_shape(mask) + tmap::tm_raster(legend.show = FALSE, col="mask", palette=c("lightblue", "darkgreen")) +
      tmap::tm_shape(commute) + tmap::tm_raster(style="cont", breaks=c(0,0.005,0.01,0.015,0.02,0.025,0.03,0.035))+tmap::tm_facets(as.layers=T)+
      tmap::tm_layout(title="Commuting Distribution", title.position = c("center","top"))

    tmap::tmap_save(tmapsavecom, paste(save_loc,"commute_dist_t_", t,".png", sep=""))

    for (k in 1:length(a)){
      commute[[k]][is.na(commute[[k]])]<-0
      commute[[k]]<-commute[[k]]*sumcom[k]
    }

  }

  # return commuting distribution and reverse foraging distribution
  comb<-list(commute,a)
  return(comb)
}
