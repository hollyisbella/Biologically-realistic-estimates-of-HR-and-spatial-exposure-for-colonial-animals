#' calculate_foraging
#'
#' A function to calculate the foraging surface
#'
#' @param mask \code{stars} object with landseamask for plotting purposes
#' @param a \code{stars} object with colony populations (one colony per layer)- created from 'create_colonies_stars_array' function
#' @param p \code{stars} object giving permeability of space - created with 'create_p_layer' function
#' @param dist \code{stars} object giving distance from each colony to accessible areas (sea), created from 'calculate_min_dist'
#' @param colonies_sf \code{sf} object created from 'create_colonies_sf' function containing colony locations, sizes and
#' IDs
#' @param col_pos \code{data frame} containing locations of colony centres, created from 'create_colonies_stars_array'
#' @param tmax \code{num} giving maximum iteration value
#' @param m \code{num} value controlling sliding scale between spatial sharing and segregation
#' @param pspill \code{num} proportion of animal usage to spill over
#' @param a2 \code{num} coefficient of density-dependent competition, a value from -0.82 to -1
#' @param a3 \code{num} value from 0 to 1 which constrains usage from moving too far from the colony centre
#' @param t_store_length \code{num} number of iterations to use to check gradient of residuals
#' @param saveloc \code{string} location to save model output images every savefreq iterations
#' @param savefreq \code{num} no. iterations between saving each image
#' @param plot_results \code{logic} value of T (TRUE) to plot, value of F (FALSE) not to plot
#' @param thresh_plot \code{num} value at which to threshold animal usage being plotted (units of probability density 0-1)
#' @param log_file_path \code{string} path to file to keep track of convergence test values
#' @param record_resids_path \code{string} path to file to keep track of residuals
#'
#' @return Returns a \code{stars} raster object containing calculated foraging surface
#' @export
#'
#' @examples
#'
calculate_foraging<-function(mask, a, p, dist, colonies_sf, col_pos, tmax,m,
                             pspill,a2, a3, t_store_length, saveloc,
                             savefreq, plot_results, thresh_plot, log_file_path,
                             record_resids_path){

  # record starting time
  start_time<-Sys.time()

  # define 9-cell neighbourhood including central cell
  nei<-cbind(c(1,0,-1,0,-1,1,-1,1,0), c(0,1,0,-1,1,1,-1,-1,0))

  # create matrix for neighbours of cell i,j
  neiij<-matrix(0,nrow=9, ncol=2)

  # specify dimensions of grid
  x_size<-unname(dim(mask))[1]
  y_size<-unname(dim(mask))[2]

  # create vector to store residuals
  residuals<-matrix(NA, ncol=t_store_length, nrow=length(a))

  # create objects to store confidence intervals of gradient of residuals for each colony
  confint_gradient_save<-matrix(data=NA,nrow=length(a), ncol=4)
  confint_gradient_save[,1]<-colonies_sf$ID

  # extract permeability layer
  pp<-p$perm

  # value of lowest permeability
  minp<-min(pp, na.rm=T)

  # select lowest permeability greater than land cells
  minpp<-min(pp[pp>1*10^-5])

  # select maximum permeability
  maxpp<-max(pp)

  # create file to record residuals
  write.csv(paste(record_resids_path,"resids.csv"))

  # create object to store results of convergence check
  check_thresh<-rep(0,length(a))

  # create objects to store values
  distnei<-matrix(0,nrow=length(a), ncol=9)
  aij<-rep(0,length(a))
  anei_store<-rep(0,length(a))

  # loop through time
  for(t in 1:tmax) {

    # convert stars object to array
    a_sum<-stars::as.tbl_cube.stars(a)
    # sum usage for each cell across colonies
    a_sum<-Reduce("+",a_sum$mets)

    # select coordinates where usage is over the carrying capacity
    coords<-which(a_sum>minpp, arr.ind=T)
    # only select cells within map (including neighbours)
    coords<-subset(coords, coords[,1]>2 &coords[,2]>2&coords[,1]<(x_size-1)&coords[,2]<(y_size-1))

    # create random order of rows of coords
    rand<-dqrng::dqsample.int(nrow(coords))

    # randomise coord rows
    coords<-coords[rand,]

    # convert to matrix
    coords<-matrix(coords, ncol=2)

    # select x and y coordinates
    coordsx<-coords[,1]
    coordsy<-coords[,2]

    # loop through all selected cells
    for(k in 1:length(coordsx)) {

      # assign i and j positions
      i<-coordsx[k]
      j<-coordsy[k]

      # define neighbourhood
      neiij[,1]<-nei[,1]+i
      neiij[,2]<-nei[,2]+j

      # permeability at current cell location
      pij<-pp[i,j]

      # create vectors to store values
      anei<-c(rep(0,(length(nei)/2)))

      for (k in 1:length(a)){
        distnei[k,]<-dist[[k]][neiij] # select distances from colony centre of cells in neighbourhood
        aij[k]<-a[[k]][i,j] # select colony usage in current cell
        anei<-anei+a[[k]][neiij] # colony usage in defined neighbourhood
        anei_store[k]<-sum(a[[k]][neiij]) # sum usage of full neighbourhood
      }

      # permeability at all neighbourhood cells
      pnei<-pp[neiij]

      # set small usage value for cells with 0 usage
      anei[anei==0]<-pij*(1/100)

      if(any(pnei<=0)) stop("problem with pnei/pp")

      # extract ratio between usage and cap
      # for impermeable cells, ratio should be very large (i.e anei should be infinitesimally small -see line 88)
      ratio<-anei/pnei

      # total usage in current cell
      Ntot<-sum(aij)

      # position of colonies in list that are present in cell ij
      col_no<-which(aij>0)

      # check if usage in cell exceeds carrying capacity
      spill<-(Ntot>(pij))*pspill

      # loop through all colonies
      for (ni in 1:length(col_no)){

        # assign colony number in list
        col<-col_no[ni]

        check_p<-exp(a2*ratio+a3*distnei[col,])

        if (all(check_p==0)){
          check_p<-pnei>minp
        }

        # dominance of current colony usage to other colonies' usage in neighbourhood
        # proportional advantage of colony controlled by segregation parameter, m
        prop<-(sum(anei_store[col]))^m/sum(anei_store^m)

        # calculate overflow of usage
        over<-aij[col]-pij*prop

        # check enough usage to supply proportion
        check_over<-(over>0)

        # current cell losing proportion of usage being passed onto other cells
        a[[col]][i,j]<-aij[col]*(1-spill)+check_over*spill*pij*prop+(1-check_over)*spill*aij[col]

        # outflux from current cell
        dn<-spill*over*check_over*(check_p/sum(check_p))

        # influx to neighbouring cells (from outflux from cell ij)
        a[[col]][neiij]<-a[[col]][neiij]+dn

      }
    }


    if (t==1){

      # if first iteration store a raster
      astore<-a

      # convert to 3 dimensions
      astore<-c(astore,along=3)

      # write blank dataframe in first time step
      resids.df<-data.frame(ID=colonies_sf$ID, resids=NA, timestep=t)

      # for subsequent iterations
    } else if (t %in% 2:(t_store_length+1)){

      # append surface to previous
      # increase dimensions of a
      acopy<-c(a, along=3)
      astore<-c(astore,acopy, along=3)

      # remove acopy object
      rm(acopy)


      for (k in 1:length(a)){

        sum_diff_a<-sum(((astore[,,,t][[k]]-astore[,,,(t-1)][[k]])/colonies_sf$median[k])^2)

        residuals[k,(t-1)]<-sum_diff_a
      }

      resids.df<-data.frame(ID=colonies_sf$ID, resids=residuals[,(t-1)], timestep=t)
      # write.table(resids.df,paste(record_resids_path,"resids.csv", sep=""),
      #             append=T,
      #             sep=",",
      #             col.names=F,
      #             row.names=F,
      #             quote=F)

    } else {


      acopy<-c(a, along=3)

      astore<-c(astore,acopy, along=3)

      # remove acopy object
      rm(acopy)

      # remove old surfaces
      astore<-astore[,,,2:(t_store_length+1)]

      # remove oldest residuals
      residuals[,1:(t_store_length-1)]<-residuals[,2:t_store_length]

      # create vector to store checked threshold values
      check_thresh<-rep(NA,length(a))

      # add new year of residuals
      for (k in 1:length(a)){

        # calculate residuals
        sum_diff_a<-sum(((astore[,,,(t_store_length-1)][[k]]-astore[,,,t_store_length][[k]])/colonies_sf$median[k])^2)

        residuals[k,t_store_length]<-sum_diff_a

        # creat data frame from which to run lm
        lm.df<-data.frame(timestep=1:t_store_length,residuals=residuals[k,1:t_store_length])

        # run linear model
        linear_model<-lm(timestep~residuals,lm.df)

        # check confidence intervals (is 0 included in gradient)
        confint<-confint(linear_model)

        confint_gradient_save[k,2]<-confint[2,1]
        confint_gradient_save[k,3]<-confint[2,2]
        confint_gradient_save[k,4]<-sum_diff_a

        # check if gradient significantly different from 0
        check<-(confint[2,1]<0&confint[2,2]>0)

        check_thresh[k]<-check

      }

      resids.df<-data.frame(ID=colonies_sf$ID, resids=residuals[,t_store_length], timestep=t)

    }

    if(all(check_thresh==1|is.na(check_thresh))) {

      sink(log_file_path, append=T)

      # print matrix
      print(confint_gradient_save)

      print(paste(maxpp, m, a2, a3))
      print(t)


      sink()

      a<-st_apply(astore, c("x", "y"), mean)

      break

    }

    #  check if plotting results
    if((t/savefreq==round(t/savefreq)| t==1) & plot_results==T) {

      sink(log_file_path, append=T)

      print(confint_gradient_save)

      sink()

      # write current residuals to file when plotting results
      write.table(resids.df,paste(record_resids_path,"resids.csv", sep=""),
                  append=T,
                  sep=",",
                  col.names=F,
                  row.names=F,
                  quote=F)

      aplot<-a
      for (k in 1:length(a)){
        # normalise foraging surface before plotting
        suma<-sum(a[[k]])
        aplot[[k]][aplot[[k]]<=thresh_plot]<-NA
        aplot[[k]]<-aplot[[k]]/suma

      }

      # plot
      tmap::tmap_mode("plot")
      tmapsavea<-tmap::tm_shape(mask) + tmap::tm_raster(legend.show = FALSE, col="mask", palette=c("lightblue", "darkgreen")) +
        tmap::tm_shape(aplot) +tmap::tm_facets(as.layers=T)+tmap::tm_raster(style="cont", breaks=c(0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01))+
        tmap::tm_layout(title="Foraging Distribution", title.position = c("center","top"))

      # save plot to specified location
      tmap::tmap_save(tmapsavea, paste(saveloc,"foraging_dist_t_", t,".png", sep=""))

      # remove plot object
      rm(aplot)

    }
  }

  # end time of foraging surface
  end_time<-Sys.time()

  # calculate duration
  time_diff<-as.numeric (end_time - start_time, units = "mins")

  sink(paste(record_resids_path,"time_taken.txt", sep=""))

  print(time_diff)

  sink()


  return(a)
}

