# FUNCTION to process colony data
# rename columns e.t.c

process_posterior_run69<-function(colonies){

  # change column names to those used in functions
  names(colonies)[2]<-"ID"

  names(colonies)[c(5,19,20)]<-c("median","longitude","latitude")

  return(colonies)

}
