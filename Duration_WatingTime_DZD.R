rm(list=ls())
library(ncdf4)
library(dplyr)

setwd('../')

ImageDirectory_dur="SPEI_RFI/DZD/DZD_withSWSI/durationDZD_ToFE_Mix/duration_DZDspei4815_srfi4815_swsi4806/"
ImageDirectory_wait="SPEI_RFI/DZD/DZD_withSWSI/waitingtimeDZD_ToFE_Mix/waitingtime_DZDspei4815_srfi4815_swsi4806/"

landfrac_nc <- nc_open('LANDFRAC4.nc')
landfrac <- ncvar_get(landfrac_nc,varid="LANDFRAC",start=c(1,1,1),count=c(288,192,1))

lats <- ncvar_get(landfrac_nc,"lat")
lons <- ncvar_get(landfrac_nc, "lon")

nlats <- dim(lats)
nlons <- dim(lons)

ToFE_spei4815_srfi4815_swsi48065_nc   <- nc_open("SPEI_RFI/DZD/DZD_withSWSI/ToFE/ToFE_FAR_spei4815_srfi4815_swsi4806_MIXno_with_TRDreservoir_afterYearReservoir.nc")

ToFE_spei4815_srfi4815_swsi48065     <- ncvar_get(ToFE_spei4815_srfi4815_swsi48065_nc,varid="ToFE")

maxlength <- 3012

for (nens in 1:100) {
  boolean_spei4815_srfi4815_swsi48065_nc    <- nc_open(paste("SPEI_RFI/DZD/DZD_withSWSI/DZDevent_mixLongerTimeTRD_and_Noreservoir/DZD_spei4815_srfi4815_swsi4806/eventDZD_spei4815_srfi4815_swsi48065_MIXno_with_TRDreservoir.",as.character(nens),".nc",sep=""))
  
  boolean_spei4815_srfi4815_swsi48065       <- ncvar_get(boolean_spei4815_srfi4815_swsi48065_nc,varid='event')
  
  waitingTime_spei4815_srfi4815_swsi48065    <- array(NA,c(nlons,nlats,3012))
  duration_spei4815_srfi4815_swsi48065       <- array(NA,c(nlons,nlats,3012))
  
  mean_waitingTime_spei4815_srfi4815_swsi48065    <- array(NA,c(nlons,nlats))
  mean_duration_spei4815_srfi4815_swsi48065       <- array(NA,c(nlons,nlats))
  
  ################################################
  ######  Waiting time spei4815_srfi4815_swsi48065
  #########################################
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      if (landfrac[jlon,jlat]>0.5 && !is.na(ToFE_spei4815_srfi4815_swsi48065[jlon,jlat])){
        
        year_ToFE <- ToFE_spei4815_srfi4815_swsi48065[jlon,jlat]
        ts_bool  <- ts(boolean_spei4815_srfi4815_swsi48065[jlon,jlat,],start=c(1850,01),frequency=12)
        data     <- window(ts_bool,start = c(year_ToFE,01),frequency=12)
        
        #        time <- (2100 - year_ToFE + 1) * 12
        #       mois <- seq(1,time,1)
        mois <- seq(1,length(data),1)
        
        df_joint <- data.frame(mois,data)
        j=0
        list_joint <- list()
        
        for (k in 1:length(data)){
          if (!is.na(df_joint[k,]$data)  && df_joint[k,]$data == 1){
            j=j+length(k)
            list_joint[j] <- df_joint[k,]$mois
            #          print(j)
          }
        }
        unlist_joint  <- unlist(list_joint)
        difs <- diff(unlist_joint)
        ###If the waiting time is 1, that means we have consecutive events
        waitingTime <- replace(difs, difs == 1, 0)
        waitingTime[is.infinite(waitingTime)] <- NA
        waitingTime[is.nan(waitingTime)] <- NA
        
        df_waitingTime <- data.frame(waitingTime)
        new_waittime <- df_waitingTime[apply(df_waitingTime!=0, 1, all),]      ### remove all 0 and consider only the value different of 0.
        new_waittime[is.infinite(new_waittime)] <- NA
        new_waittime[is.nan(new_waittime)] <- NA
        
        
        waitingTime_spei4815_srfi4815_swsi48065[jlon,jlat,]   <- c(new_waittime,rep(NA, maxlength - length(new_waittime)))
        mean_waitingTime_spei4815_srfi4815_swsi48065[jlon,jlat] <- mean(new_waittime)
        #print(jlon)
      }
    }
  }
  
  
  ##########################################
  ######### Duration spei4815_srfi4815_swsi48065
  #########################################
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      if (landfrac[jlon,jlat]>0.5 && !is.na(ToFE_spei4815_srfi4815_swsi48065[jlon,jlat])){
        
        year_ToFE <- ToFE_spei4815_srfi4815_swsi48065[jlon,jlat]
        ts_bool  <- ts(boolean_spei4815_srfi4815_swsi48065[jlon,jlat,],start=c(1850,01),frequency=12)
        data     <- window(ts_bool,start = c(year_ToFE,01),frequency=12)
        
        duration <- rle(as.numeric(data))
        durLeng <- duration$lengths
        durVal  <- duration$values
        df_dur <- data.frame(durVal,durLeng)
        
        list_duration <- list()
        for (k in 1:length(durLeng)){
          if (!is.na(df_dur[k,]$durVal)  && df_dur[k,]$durVal==1){
            j=j+length(k)
            list_duration[j] <- df_dur[k,]$durLeng
            #          print(j)
          }
        }
        duration <- unlist(list_duration)
        duration_spei4815_srfi4815_swsi48065[jlon,jlat,]    <- c(duration,rep(NA, maxlength - length(duration)))
        mean_duration_spei4815_srfi4815_swsi48065[jlon,jlat] <- mean(duration)
        #print(jlon)
      }
    }
  }
  
  
  
  ##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##################
  ##### save mean waiting time and Duration 
  ##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##################
  
  xvals <- lons
  yvals <- lats
  
  nx <- length(xvals)
  ny <- length(yvals)
  
  xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
  ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
  
  mv <- NA     # missing value
  
  var_prec1 <- ncvar_def("duration", 'duration', list(xdim,ydim), mv )
  var_prec2 <- ncvar_def("waitingtime", 'waitingtime', list(xdim,ydim), mv )
  
  ######Save to netcdf 
  ncout1 <- nc_create(paste(ImageDirectory_dur,"mean_duration_spei4815_srfi4815_swsi48065_MIXno_with_TRDreservoir.",as.character(nens),".nc", sep = ""), list(var_prec1),force_v4=TRUE)
  ncout2 <- nc_create(paste(ImageDirectory_wait,"mean_waitingTime_spei4815_srfi4815_swsi48065_MIXno_with_TRDreservoir.",as.character(nens),".nc", sep = ""), list(var_prec2),force_v4=TRUE)
  # put variables
  
  ncvar_put(ncout1, var_prec1, mean_duration_spei4815_srfi4815_swsi48065, start=c(1,1), count=c(nx,ny), verbose=TRUE)
  ncvar_put(ncout2, var_prec2, mean_waitingTime_spei4815_srfi4815_swsi48065, start=c(1,1), count=c(nx,ny), verbose=TRUE)
  
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  
  ncatt_put(ncout2,"lon","axis","X") 
  ncatt_put(ncout2,"lat","axis","Y")
  
  
  # get global attributes
  title        <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"DZD Duration and Waiting Time")
  institution  <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"institution")
  datasource   <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"source")
  references   <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"references")
  history      <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"history")
  Conventions  <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout1,0,"title","DZD Duration")
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","SPEI48_RFI48")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  ncatt_put(ncout2,0,"title","DZD Waiting Time")
  ncatt_put(ncout2,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout2,0,"source",datasource$value)
  ncatt_put(ncout2,0,"references","SPEI48_RFI48")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout2,0,"history",history)
  ncatt_put(ncout2,0,"Conventions",Conventions$value)
  nc_close(ncout2)
  
  
  ##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##################
  ##### save waiting time and Duration time series 
  ##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##################
  
  ########
  time <- ncvar_get(boolean_spei4815_srfi4815_swsi48065_nc, "time")
  tunits <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,"time","units")
  calendar <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,"time","calendar")$value
  
  tdim <- ncdim_def("time", tunits$value, as.numeric(time),unlim=T,calendar = calendar)
  
  mv <- NA     # missing value
  var_prec11 <- ncvar_def("duration", 'duration', list(xdim,ydim,tdim), mv )
  var_prec12 <- ncvar_def("waitingtime", 'waitingtime', list(xdim,ydim,tdim), mv )
  
  nt <- length(time) # Imagine this is actually some very large number of timesteps
  
  ######Save to netcdf 
  ncout11 <- nc_create(paste(ImageDirectory_dur,"duration_spei4815_srfi4815_swsi48065_MIXno_with_TRDreservoir.",as.character(nens),".nc", sep = ""), list(var_prec11),force_v4=TRUE)
  ncout12 <- nc_create(paste(ImageDirectory_wait,"waitingTime_spei4815_srfi4815_swsi48065_MIXno_with_TRDreservoir.",as.character(nens),".nc", sep = ""), list(var_prec12),force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout11, var_prec11, duration_spei4815_srfi4815_swsi48065, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  ncvar_put(ncout12, var_prec12, waitingTime_spei4815_srfi4815_swsi48065, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout11,"lon","axis","X") 
  ncatt_put(ncout11,"lat","axis","Y")
  ncatt_put(ncout11,"time","axis","T")
  
  ncatt_put(ncout12,"lon","axis","X") 
  ncatt_put(ncout12,"lat","axis","Y")
  ncatt_put(ncout12,"time","axis","T")
  
  # get global attributes
  title        <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"DZD")
  institution  <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"institution")
  datasource   <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"source")
  references   <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"references")
  history      <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"history")
  Conventions  <- ncatt_get(boolean_spei4815_srfi4815_swsi48065_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout11,0,"title","DZD Duration")
  ncatt_put(ncout11,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout11,0,"source",datasource$value)
  ncatt_put(ncout11,0,"references","SPEI48_RFI48")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout11,0,"history",history)
  ncatt_put(ncout11,0,"Conventions",Conventions$value)
  nc_close(ncout11)
  
  ncatt_put(ncout12,0,"title","DZD Waiting Time")
  ncatt_put(ncout12,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout12,0,"source",datasource$value)
  ncatt_put(ncout12,0,"references","SPEI48_RFI48")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout12,0,"history",history)
  ncatt_put(ncout12,0,"Conventions",Conventions$value)
  nc_close(ncout12)
  
  nc_close(boolean_spei4815_srfi4815_swsi48065_nc)
  print(nens)
}

