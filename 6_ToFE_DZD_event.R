rm(list=ls())
library(ncdf4)
library(zoo)
library(rainfarmr)
library(MBC)
library(SPEI)
#library(fitur)
#library(actuar)
library(dplyr)
library(tidyr)

setwd('./')

ImageDirectory="SPEI_SRFI_SWSI/ToFE_DZD_event/"

####### we use the ensemble sum of the DZD event (period 1850-2100)
enssum_DZD_event_nc      <- nc_open("SPEI_SRFI_SWSI/DZD_event/enssum_DZD_event_spei4815_srfi4815_swsi4806_TRD.nc")
enssum_DZD_event         <- ncvar_get(enssum_DZD_event_nc,varid="dzd")

lats <- ncvar_get(enssum_DZD_event_nc,"lat")
lons <- ncvar_get(enssum_DZD_event_nc, "lon")

nlats <- dim(lats)
nlons <- dim(lons)


#####################
reservervoir_year_nc <- nc_open("reservoir/GRanD_Version_1_3/time_todry/reservoir_year_cap.nc")
reservervoir_year    <- ncvar_get(reservervoir_year_nc,varid="year_reser")

##############################################################################################################################################################################################
############################################### FAR FAR FAR spei4815_srfi4815_swsi4806 ################################################################################################################################
##############################################################################################################################################################################################

FAR_DZD_event         <- array(NA,c(nlons,nlats,20))
ToFE_DZD_event        <- array(NA,c(nlons,nlats))

for (jlon in 1:nlons){
  for (jlat in 1:nlats){
    
    mean_event <- mean(enssum_DZD_event[jlon,jlat,],na.rm=TRUE)
    if (!is.na(mean_event)){
      
      ts_enssum_DZD_event     <- ts(enssum_DZD_event[jlon,jlat,],start=c(1850,1),frequency=12)
      ww_enssum_DZD_event     <- window(ts_enssum_DZD_event,start=c(1850,1), end=c(1899,12), frequency=12)
      
      prob_couterfactual_enssum_DZD_event         <- sum(ww_enssum_DZD_event, na.rm=TRUE)/60000
      
      annee <- seq(1900,2099,10)
      list_enssum_DZD_event_decad <- list()
      jj=0
      for (decad1 in seq(1900,2099,10)){
        decad2 <- decad1 + 9
        jj = jj +1
        list_enssum_DZD_event_decad[jj]     <- sum(window(ts_enssum_DZD_event,start=c(decad1,1), end=c(decad2,12), frequency=12), na.rm=TRUE)
        
      }
      enssum_DZD_event_decad        <- unlist(list_enssum_DZD_event_decad)
      
      #### calculation of probability
      prob_factual_enssum_DZD_event_decad         <- enssum_DZD_event_decad/12000
      
      ### prob infinite to NA
      prob_couterfactual_enssum_DZD_event[is.infinite(prob_couterfactual_enssum_DZD_event)] <- NA
      prob_factual_enssum_DZD_event_decad[is.infinite( prob_factual_enssum_DZD_event_decad)] <- NA
      
      
      ##### Change nan to NA
      prob_couterfactual_enssum_DZD_event[is.nan(prob_couterfactual_enssum_DZD_event)] <- NA
      prob_factual_enssum_DZD_event_decad[is.nan(prob_factual_enssum_DZD_event_decad)] <- NA
      
      
      ###### FAR
      cal_FAR_enssum_DZD_event       <- 1- (prob_couterfactual_enssum_DZD_event/prob_factual_enssum_DZD_event_decad)
      
      ####### to NA
      cal_FAR_enssum_DZD_event[is.infinite(cal_FAR_enssum_DZD_event)] <- NA
      cal_FAR_enssum_DZD_event[is.nan(cal_FAR_enssum_DZD_event)] <- NA
      
      
      ########## to array
      FAR_DZD_event[jlon,jlat,]         <- cal_FAR_enssum_DZD_event
      
    }
  }
  #print(jlon)
  
}  
#################
#######################################
##################Save FAR to netcdf
###################################

time        <- ncvar_get(enssum_DZD_event_nc, "time")
tunits      <- ncatt_get(enssum_DZD_event_nc,"time","units")
calendar    <- ncatt_get(enssum_DZD_event_nc,"time","calendar")$value

xvals <- lons
yvals <- lats
tvals <- seq(1900,2090,10)

nx <- length(xvals)
ny <- length(yvals)
nt <- length(tvals)

xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
tdim <- ncdim_def( 'time', tunits$value, tvals)
#tdim <- ncdim_def("time", tunits$value, as.numeric(1),unlim=T,calendar = calendar)

mv <- NA     # missing value
var_prec <- ncvar_def("FAR", 'FAR', list(xdim,ydim,tdim), mv )


######
ncout1 <- nc_create(paste(ImageDirectory,"FAR_enssum_DZD_event.nc", sep = ""), list(var_prec),force_v4=TRUE)

# put variables
ncvar_put(ncout1, var_prec, FAR_DZD_event, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
# put additional attributes into dimension and data variables
ncatt_put(ncout1,"lon","axis","X") 
ncatt_put(ncout1,"lat","axis","Y")
ncatt_put(ncout1,"time","axis","T")


# get global attributes
title        <- ncatt_get(enssum_DZD_event_nc,0,"ToFE")
institution  <- ncatt_get(enssum_DZD_event_nc,0,"institution")
datasource   <- ncatt_get(enssum_DZD_event_nc,0,"source")
references   <- ncatt_get(enssum_DZD_event_nc,0,"references")
history      <- ncatt_get(enssum_DZD_event_nc,0,"history")
Conventions <- ncatt_get(enssum_DZD_event_nc,0,"Conventions")

# add global attributes
ncatt_put(ncout1,0,"title","FAR p(SPEI48 <= -1.5,srfi48 <= -1.5, swsi48 <= 0.6,no and with TRD reservoir)")
ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
ncatt_put(ncout1,0,"source",datasource$value)
ncatt_put(ncout1,0,"references","ERA5")
history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
ncatt_put(ncout1,0,"history",history)
ncatt_put(ncout1,0,"Conventions",Conventions$value)
nc_close(ncout1)



##############################################################################################################################################################################################
###############################################ToFE ToFE ToFE spei4815_srfi4815_swsi4806 ################################################################################################################################
##############################################################################################################################################################################################

for (jlon in 1:nlons){
  for (jlat in 1:nlats){
    mean_FAR_enssum_DZD_event <- mean(FAR_DZD_event[jlon,jlat,],na.rm=TRUE)
    if ( !is.na(mean_FAR_enssum_DZD_event)){
      
      vec_FAR_DZD_event     <- FAR_DZD_event[jlon,jlat,]
      annee                 <- seq(1900,2099,10)
      df_FAR_DZD_event      <- data.frame(annee,vec_FAR_DZD_event)
      
      
      j=0
      list_year_emergence <- list()
      for (i in 1:20){
        if (!is.na(df_FAR_DZD_event[i,]$vec_FAR_DZD_event) && df_FAR_DZD_event[i,]$vec_FAR_DZD_event >= 0.99){
          j=j+length(i)
          list_year_emergence[j] <- df_FAR_DZD_event[i,]$annee
          
        }
      }
      unlist_year_emergence <- unlist(list_year_emergence)
      
      
      
      if (!is.na(reservervoir_year[jlon,jlat])){
        year_resr <- reservervoir_year[jlon,jlat]
        year_DZDevent  <- unlist_year_emergence[unlist_year_emergence > year_resr]
        
        if (!is.null(year_DZDevent[1])){
               
          ToFE_DZD_event[jlon,jlat] <- year_DZDevent[1]
          
        }
        else if(is.null(year_DZDevent[1]) && !is.null(unlist_year_emergence[1])){
          ToFE_DZD_event[jlon,jlat] <- unlist_year_emergence[1]
        }
        
      }
      else if(is.na(reservervoir_year[jlon,jlat]) && !is.null(unlist_year_emergence[1])) {
        ToFE_DZD_event[jlon,jlat] <- unlist_year_emergence[1]
      }
      else{ToFE_DZD_event[jlon,jlat] <- NA}
      #print(jlon)
      
    }
  }
}

#######################################
##################Save to netcdf
###################################

time        <- ncvar_get(enssum_DZD_event_nc, "time")
tunits      <- ncatt_get(enssum_DZD_event_nc,"time","units")
calendar    <- ncatt_get(enssum_DZD_event_nc,"time","calendar")$value

xvals <- lons
yvals <- lats

nx <- length(xvals)
ny <- length(yvals)
xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
#tdim <- ncdim_def("time", tunits$value, as.numeric(1),unlim=T,calendar = calendar)

mv <- NA     # missing value
var_prec1 <- ncvar_def("ToFE", 'ToFE', list(xdim,ydim), mv )

######Save to netcdf 
ncout1 <- nc_create(paste(ImageDirectory,"ToFE_DZD_event.nc", sep = ""), list(var_prec1),force_v4=TRUE)

# put variables
ncvar_put(ncout1, var_prec1, ToFE_DZD_event, start=c(1,1), count=c(nx,ny), verbose=TRUE)

# put additional attributes into dimension and data variables
ncatt_put(ncout1,"lon","axis","X") 
ncatt_put(ncout1,"lat","axis","Y")


# get global attributes
title        <- ncatt_get(enssum_DZD_event_nc,0,"ToFE")
institution  <- ncatt_get(enssum_DZD_event_nc,0,"institution")
datasource   <- ncatt_get(enssum_DZD_event_nc,0,"source")
references   <- ncatt_get(enssum_DZD_event_nc,0,"references")
history      <- ncatt_get(enssum_DZD_event_nc,0,"history")
Conventions <- ncatt_get(enssum_DZD_event_nc,0,"Conventions")

# add global attributes
ncatt_put(ncout1,0,"title","ToFE p(SPEI48 <= -1.5,srfi48 <= -1.5, swsi48 <= 0.6,no and with TRD reservoir)")
ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
ncatt_put(ncout1,0,"source",datasource$value)
ncatt_put(ncout1,0,"references","ERA5")
history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
ncatt_put(ncout1,0,"history",history)
ncatt_put(ncout1,0,"Conventions",Conventions$value)
nc_close(ncout1)

#nc_close(freq_couterfactual_spei48_srfi4816_swsi4806_nc)
#nc_close(enssum_DZD_event_nc)
