rm(list=ls())
library(ncdf4)
library(rainfarmr)
library(rle)
library(data.table)
library(magrittr)
library(dplyr)


#########################
####for all grid point
setwd('../')
ImageDirectory="SPEI_SRFI_SWSI/DZD_event/"

duration_dry_nc   <- nc_open('reservoir/GRanD_Version_1_3/time_todry/reservoir_time_todry_duration.nc')
duration_dry      <- ncvar_get(duration_dry_nc,varid="trd") 

#############  ## CNRM selected for differents scenario ssp245 and ssp370
#for (nscen in c("ssp245", "ssp370","ssp585")){

#models <- c("CM6-1_",as.character(nscen),"_r1i1p1f2","CM6-1_",as.character(nscen),"_r2i1p1f2","CM6-1_",as.character(nscen),"_r3i1p1f2","CM6-1_",as.character(nscen),"_r4i1p1f2",       
#            "CM6-1_",as.character(nscen),"_r5i1p1f2","CM6-1_",as.character(nscen),"_r6i1p1f2","ESM2-1_",as.character(nscen),"_r14i1p1f2", "ESM2-1_",as.character(nscen),"_r15i1p1f2",
#            "ESM2-1_",as.character(nscen),"_r1i1p1f2","ESM2-1_",as.character(nscen),"_r2i1p1f2","ESM2-1_",as.character(nscen),"_r3i1p1f2","ESM2-1_",as.character(nscen),"_r4i1p1f2",
#            "ESM2-1_",as.character(nscen),"_r5i1p1f2","CM6-1-HR_",as.character(nscen),"_r1i1p1f2","CM6-1_",as.character(nscen),"_r1i1p1f2")

models <- 1:100    #### CESM2-LE

maxlength = 3012   ### timestep

for (nens in models) {
  srfi48_nc     <- nc_open(paste("SPEI_SRFI_SWSI/srfi48/srfi48_corrected_1850_2100.",as.character(nens),".nc",sep=""))
  spei48_nc     <- nc_open(paste("SPEI_SRFI_SWSI/spei48/spei4845_corrected_1850_2100.",as.character(nens),".nc",sep=""))
  swsi48_nc     <- nc_open(paste("SPEI_SRFI_SWSI/swsi48/swsi48_correctedInput_1850_2100.",as.character(nens),".nc",sep=""))
  
  srfi48                   <- ncvar_get(srfi48_nc,varid='srfi48')
  spei48                   <- ncvar_get(spei48_nc,varid='spei48')
  swsi48                   <- ncvar_get(swsi48_nc,varid='swsi48')
  
  lats <- ncvar_get(spei48_nc,"lat")
  lons <- ncvar_get(spei48_nc, "lon")
  
  nlats <- dim(lats)
  nlons <- dim(lons)

  ### create array
  
  boolean_spei4815_srfi4815_swsi4806         <- array(NA,dim(spei48))
  DZD_event_spei4815_srfi4815_swsi4806       <- array(NA,dim(spei48))
  
  ################################################################################
  ###################### SPEI48 <= -1.5 && SRFI48 <= -1.5 && SWSI48 <= 0.6
  #########################################################################################
  
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      mean_spei48 <-  mean(spei48[jlon,jlat,], na.rm=TRUE)
      mean_srfi48 <-  mean(srfi48[jlon,jlat,], na.rm=TRUE)
      mean_swsi48 <-  mean(swsi48[jlon,jlat,], na.rm=TRUE)
      
      if (!is.na(mean_spei48) && !is.na(mean_srfi48) && !is.na(mean_swsi48)){
        
        tser_spei48 <- spei48[jlon,jlat,]
        tser_srfi48 <- srfi48[jlon,jlat,]
        tser_swsi48 <- swsi48[jlon,jlat,]
        
        tser_spei48[is.infinite(tser_spei48)] <- NA
        tser_spei48[is.nan(tser_spei48)] <- NA
        
        tser_srfi48[is.infinite(tser_srfi48)] <- NA
        tser_srfi48[is.nan(tser_srfi48)] <- NA
        
        tser_swsi48[is.infinite(tser_swsi48)] <- NA
        tser_swsi48[is.nan(tser_swsi48)] <- NA
        
        ###########
        
        list_joint <- list()
        for (k in 1:3012){
          list_joint[k] <- ifelse(tser_spei48[k] <= -1.5 && tser_srfi48[k] <= -1.5 && tser_swsi48[k] <= 0.6, 1, 0)
        }
        unlist_joint  <- unlist(list_joint)
        unlist_joint2 <- as.numeric(unlist_joint)
        boolean_spei4815_srfi4815_swsi4806[jlon,jlat,] <- unlist_joint2
      }
    }
    
    #print(jlon)
  } 
  
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      mean_bool <-  mean(boolean_spei4815_srfi4815_swsi4806[jlon,jlat,], na.rm=TRUE)
      if (!is.na(mean_bool)){
        ts_boolean_spei4815_srfi4815_swsi4806 <- boolean_spei4815_srfi4815_swsi4806[jlon,jlat,]
        ts_boolean_spei4815_srfi4815_swsi4806[is.infinite(ts_boolean_spei4815_srfi4815_swsi4806)] <- NA
        ts_boolean_spei4815_srfi4815_swsi4806[is.nan(ts_boolean_spei4815_srfi4815_swsi4806)] <- NA
        ts_boolean_spei4815_srfi4815_swsi4806[is.na(ts_boolean_spei4815_srfi4815_swsi4806)] <- 0
        
        
        rle_boolean_spei4815_srfi4815_swsi4806    <- rle(ts_boolean_spei4815_srfi4815_swsi4806)
        df_rle_boolean_spei4815_srfi4815_swsi4806 <- data.frame(rle_boolean_spei4815_srfi4815_swsi4806$values,rle_boolean_spei4815_srfi4815_swsi4806$lengths)
        
        #########rleid() generates the ids or repeated group of equal length to the original vector
        
        id_bool <- rleid(ts_boolean_spei4815_srfi4815_swsi4806)
        
        rep_id_bool <- list() 
        for (i in 1:3012){
          rep_id_bool[i] <- df_rle_boolean_spei4815_srfi4815_swsi4806$rle_boolean_spei4815_srfi4815_swsi4806.lengths[id_bool[i]]
        }
        
        unlist_rep_id_bool        <- unlist(rep_id_bool)
        df_rle_SPEI48_RFI48lessminus1   <- data.frame(ts_boolean_spei4815_srfi4815_swsi4806,unlist_rep_id_bool) 
        
        eventDZD <- list()
        for (k in 1:3012){
          if(df_rle_SPEI48_RFI48lessminus1$ts_boolean_spei4815_srfi4815_swsi4806[k]==1 && !is.na(duration_dry[jlon,jlat]) && df_rle_SPEI48_RFI48lessminus1$unlist_rep_id_bool[k] >= duration_dry[jlon,jlat]){
            eventDZD[k] <- 1
          }
          else if(df_rle_SPEI48_RFI48lessminus1$ts_boolean_spei4815_srfi4815_swsi4806[k]==1 && is.na(duration_dry[jlon,jlat])){
            
            eventDZD[k] <- 1
          }
          else{
            eventDZD[k] <- 0
            
          }
        }
        
        unlist_eventDZD          <- unlist(eventDZD)
        DZD_event_spei4815_srfi4815_swsi4806[jlon,jlat,] <- unlist_eventDZD
      }
    }
    #print(jlon)
  } 
  
  #######################################
  ##################Save to netcdf 
  ###################################
  
  
  time        <- ncvar_get(spei48_nc, "time")
  tunits      <- ncatt_get(spei48_nc,"time","units")
  calendar    <- ncatt_get(spei48_nc,"time","calendar")$value
  
  xvals <- lons
  yvals <- lats
  
  nx <- length(xvals)
  ny <- length(yvals)
  nt <- length(time)    # Imagine this is actually some very large number of timesteps
  
  xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
  ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
  tdim <- ncdim_def("time", tunits$value, as.numeric(time),unlim=T,calendar = calendar)
  
  mv <- NA     # missing value
  
  var_prec1 <- ncvar_def("dzd", 'dzd', list(xdim,ydim,tdim), mv )
  
  
  ######Save to netcdf 
  ncout1 <- nc_create(paste(ImageDirectory,"DZD_event_spei4815_srfi4815_swsi4806_TRD.",as.character(nens),".nc", sep = ""), list(var_prec1),force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout1, var_prec1, DZD_event_spei4815_srfi4815_swsi4806, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  ncatt_put(ncout1,"time","axis","T")
  
  # get global attributes
  title        <- ncatt_get(spei48_nc,0,"DZD event")
  institution  <- ncatt_get(spei48_nc,0,"institution")
  datasource   <- ncatt_get(spei48_nc,0,"source")
  references   <- ncatt_get(spei48_nc,0,"references")
  history      <- ncatt_get(spei48_nc,0,"history")
  Conventions <- ncatt_get(spei48_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout1,0,"title",paste("DZD event p(SPEI48 <= -1.5,SRFI48 <= -1.5, SWSI48 <= 0.6,with and without TRD of reservoir) from model",as.character(nens), sep = ""))
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","ERA5")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  
  nc_close(srfi48_nc)
  nc_close(spei48_nc)
  nc_close(swsi48_nc)
  
  print(nens)
}
