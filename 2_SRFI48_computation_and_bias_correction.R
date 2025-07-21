rm(list=ls())
library(ncdf4)
library(MBC)
library(SPEI)
library(dplyr)
library(zoo)

setwd('./')

ImageDirectory="SPEI_SRFI_SWSI/srfi48/"


#############  ## CNRM selected for differents scenario ssp245 and ssp370
#for (nscen in c("ssp245", "ssp370","ssp585")){

#models <- c("CM6-1_",as.character(nscen),"_r1i1p1f2","CM6-1_",as.character(nscen),"_r2i1p1f2","CM6-1_",as.character(nscen),"_r3i1p1f2","CM6-1_",as.character(nscen),"_r4i1p1f2",       
#            "CM6-1_",as.character(nscen),"_r5i1p1f2","CM6-1_",as.character(nscen),"_r6i1p1f2","ESM2-1_",as.character(nscen),"_r14i1p1f2", "ESM2-1_",as.character(nscen),"_r15i1p1f2",
#            "ESM2-1_",as.character(nscen),"_r1i1p1f2","ESM2-1_",as.character(nscen),"_r2i1p1f2","ESM2-1_",as.character(nscen),"_r3i1p1f2","ESM2-1_",as.character(nscen),"_r4i1p1f2",
#            "ESM2-1_",as.character(nscen),"_r5i1p1f2","CM6-1-HR_",as.character(nscen),"_r1i1p1f2","CM6-1_",as.character(nscen),"_r1i1p1f2")

models <- 1:100    #### CESM2-LE

maxlength = 3012   ### timestep 
##########################

for (nens in models) {
  river_nc   <- nc_open(paste("River/river_monmean",as.character(nens),".nc",sep=""))
  river      <- ncvar_get(river_nc,varid='river')
  
  #####ERA data
  river_era_nc    <- nc_open("ERA/monRiverEra_GloFAS_mosartgrid_1979_2022.nc")
  river_era       <- ncvar_get(river_era_nc,varid='dis')
  
  #############

  lats <- ncvar_get(river_nc,"lat")
  lons <- ncvar_get(river_nc, "lon")
  
  nlats <- dim(lats)
  nlons <- dim(lons)
  
  maxlength = 3012   ### timestep
  
  ###################### create array
  
  ssi48_era       <- array(NA,dim(river_era))
  ssi48_cesm      <- array(NA,dim(river))
  ssi48_calib     <- array(NA,dim(river))
  na_ssi48_calib  <- array(NA,dim(river))
  
  
  ########################################
  
  ############# calibration
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      land_riv_era <-  mean(river_era[jlon,jlat,], na.rm=TRUE)
      land_riv_cesm <-  mean(river[jlon,jlat,], na.rm=TRUE)
      
      if (!is.na(land_riv_era) && !is.na(land_riv_cesm)){
        
        ###########################
        ts_river_era      <- ts(river_era[jlon,jlat,],start=c(1979,1),frequency=12)
        ts_river          <- ts(river[jlon,jlat,],start=c(1850,1),frequency=12)
        
        ssi_era48         <- spi(ts_river_era,48,ref.start=c(1979,1),ref.end=c(2014,12),frequency=12,na.rm=TRUE,verbose = FALSE)
        ssi_cesm48        <- spi(ts_river,48,ref.start=c(1979,1),ref.end=c(2014,12),frequency=12,na.rm=TRUE,verbose = FALSE)
        
        fssi_era48        <- as.numeric(ssi_era48$fitted)
        fssi_cesm48       <- as.numeric(ssi_cesm48$fitted)
        
        fssi_era48[is.infinite(fssi_era48)]<-NA
        fssi_cesm48[is.infinite(fssi_cesm48)]<-NA
        
        ssi48_era[jlon,jlat,]   <- fssi_era48
        ssi48_cesm[jlon,jlat,]  <- fssi_cesm48
        
        
        #######
        ###QDM bias correction
        
        in_ssi48_cesm           <- ts(ssi48_cesm[jlon,jlat,],start=c(1850,1),frequency=12)  
        in_calib_ssi48_cesm     <- window(in_ssi48_cesm,start=c(1979,1),end=c(2022,12),frequency=12)  
        in_ssi48_era            <- ts(ssi48_era[jlon,jlat,],start=c(1979,1),frequency=12)  
        
        if (sum(!is.na(in_ssi48_cesm))>2 && sum(!is.na(in_calib_ssi48_cesm))>2 && sum(!is.na(in_ssi48_era))>2){
          proj48  <- na.approx(in_ssi48_cesm)
          obs48   <- na.approx(in_ssi48_era)
          calib48  <- na.approx(in_calib_ssi48_cesm)
          qdm_calib48 <- QDM(obs48, calib48, proj48, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,jitter.factor=0, n.tau=NULL, ratio.max=2,ratio.max.trace=10*trace, ECBC=FALSE, ties='first',subsample=NULL, pp.type=7)
        }
        
        
        ssi48_calib    <- qdm_calib48$mhat.p
        
        ##### Complete to have the same length
        na_ssi48_calib[jlon,jlat,]    <- c(rep(NA, maxlength - length(ssi48_calib)),ssi48_calib)
        
        
        #############################
        #print(jlon)
      }
    }
  }
  
  
  #######################################
  ##################Save to netcdf
  ###################################
  
  time <- ncvar_get(river_nc, "time")
  tunits <- ncatt_get(river_nc,"time","units")
  calendar <- ncatt_get(river_nc,"time","calendar")$value
  
  xvals <- lons
  yvals <- lats
  
  nx <- length(xvals)
  ny <- length(yvals)
  xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
  ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
  tdim <- ncdim_def("time", tunits$value, as.numeric(time),unlim=T,calendar = calendar)
  
  mv <- NA     # missing value
  var_prec1 <- ncvar_def("srfi48", 'srfi48', list(xdim,ydim,tdim), mv )
  
  nt <- length(time ) # Imagine this is actually some very large number of timesteps
  
  ######Save to netcdf 
  ncout1 <- nc_create(paste(ImageDirectory,"srfi48_corrected_1850_2100.",as.character(nens),".nc", sep = ""), list(var_prec1),force_v4=TRUE)
  # put variables
  
  ncvar_put(ncout1, var_prec1, na_ssi48_calib, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  ncatt_put(ncout1,"time","axis","T")
  
  
  # get global attributes
  title        <- ncatt_get(river_nc,0,"SRFI48 Bias corrected")
  institution  <- ncatt_get(river_nc,0,"institution")
  datasource   <- ncatt_get(river_nc,0,"source")
  references   <- ncatt_get(river_nc,0,"references")
  history      <- ncatt_get(river_nc,0,"history")
  Conventions  <- ncatt_get(river_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout1,0,"title",paste("Standardized River Flow Index from model ",as.character(nens), sep = ""))
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","ERA5")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  nc_close(river_nc)
  nc_close(river_era_nc)
  print(nens)
  
}


