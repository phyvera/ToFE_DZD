rm(list=ls())
library(ncdf4)
library(MBC)
library(SPEI)
library(zoo)
library(dplyr)
library(tidyr)
library(stats)

setwd('./')

ImageDirectory_spei48 ="SPEI_SRFI_SWSI/spei48/"

#############################################################
#####ERA data

prec_era_nc <- nc_open("ERA/monthly_total_precipitation_era5_sum_remapcon.nc")
tmin_era_nc <- nc_open("ERA/monTmin.era5_remapcon.nc")
tmax_era_nc <- nc_open("ERA/monTmax.era5_remapcon.nc")

prec_era <- ncvar_get(prec_era_nc,varid="tp")
tmin_era <- ncvar_get(tmin_era_nc,varid="mn2t")
tmax_era <- ncvar_get(tmax_era_nc,varid="mx2t")


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
  
  prec_nc <- nc_open(paste("precipitation/pr_monsum",as.character(nens),".nc",sep=""))
  tmax_nc <- nc_open(paste("Tmax/tmax_mon",as.character(nens),".nc",sep=""))
  tmin_nc <- nc_open(paste("Tmin/tmin_mon",as.character(nens),".nc",sep=""))
  
  
  prec <- ncvar_get(prec_nc,varid="pr")
  tmin <- ncvar_get(tmin_nc,varid="tmin")
  tmax <- ncvar_get(tmax_nc,varid="tmax")
  
  lats <- ncvar_get(prec_nc,"lat")
  lons <- ncvar_get(prec_nc, "lon")
  
  nlats <- dim(lats)
  nlons <- dim(lons)
  
  ###### create array
  
  pet_hgera         <- array(NA,dim(tmin_era))
  wb_hgera          <- array(NA,dim(tmin_era))
  pet_hgcmip        <- array(NA,dim(tmin))
  wb_hgcmip         <- array(NA,dim(tmin))
  spei48_hgcmip     <- array(NA,dim(tmin))
  spei48_hgera      <- array(NA,dim(tmin_era))
  spei48_hgcalib    <- array(NA,dim(tmin))
  na_spei48_hgcalib <- array(NA,dim(tmin))
  
  ########################################
  
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      
      prec_mean_era       <-  mean(prec_era[jlon,jlat,], na.rm=TRUE)
      prec_mean           <-  mean(prec[jlon,jlat,], na.rm=TRUE)
      
      if (!is.na(prec_mean_era) && !is.na(prec_mean)){
        
        pet_hgera[jlon,jlat,]  <- hargreaves(Tmin = tmin_era[jlon,jlat,], Tmax = tmax_era[jlon,jlat,],lat = lats[jlat],na.rm =TRUE)
        pet_hgcmip[jlon,jlat,] <- hargreaves(Tmin = (tmin[jlon,jlat,] - 273.15), Tmax = (tmax[jlon,jlat,]-273.15),lat = lats[jlat],na.rm =TRUE) ### temeperature to degree C
        
        wb_hgera[jlon,jlat,]  <- prec_era[jlon,jlat,]  - pet_hgera[jlon,jlat,]
        wb_hgcmip[jlon,jlat,] <- (prec[jlon,jlat,] * 86400) - pet_hgcmip[jlon,jlat,] # precipitation to mm/day
        
        
        ####################SPEI#######
        fwb_hgera     <- ts(wb_hgera[jlon,jlat,],start=c(1979,1),frequency=12)
        fwb_hgcmip    <- ts(wb_hgcmip[jlon,jlat,],start=c(1850,1),frequency=12)
        
        ##cmip
        spei_hgcmip         <- spei(fwb_hgcmip,48,ref.start=c(1979,1),ref.end=c(2014,12),frequency=12,na.rm=TRUE)
        fspei48_hgcmip      <- as.numeric(spei_hgcmip$fitted)
        fspei48_hgcmip[is.infinite(fspei48_hgcmip)]<-NA
        spei48_hgcmip[jlon,jlat,]  <- fspei48_hgcmip
        
        
        ##ERA     
        spei_hgera48      <- spei(fwb_hgera,48,ref.start=c(1979,1),ref.end=c(2014,12),frequency=12,verbose=FALSE, na.rm=TRUE)
        fspei48_hgera    <- as.numeric(spei_hgera48$fitted)
        fspei48_hgera[is.infinite(fspei48_hgera)]<-NA
        spei48_hgera[jlon,jlat,]  <- fspei48_hgera
        
        
        ####### #QDM correction SPEI
        
        
        in_spei48_hgcmip        <- ts(spei48_hgcmip[jlon,jlat,],start=c(1850,1),frequency=12)  
        in_calib_spei48_hgcmip  <- window(in_spei48_hgcmip,start=c(1979,1),end=c(2022,12),frequency=12)  
        in_spei48_hgera          <- ts(spei48_hgera[jlon,jlat,],start=c(1979,1),frequency=12)  
        
        if (sum(!is.na(in_spei48_hgcmip))>4 && sum(!is.na(in_calib_spei48_hgcmip))>4 && sum(!is.na(in_spei48_hgera))>4){
          proj48_hg       <- na.approx(in_spei48_hgcmip)
          obs48_hg        <- na.approx(in_spei48_hgera)
          calib48_hg      <- na.approx(in_calib_spei48_hgcmip)
          qdm_calib48_hg <- QDM(obs48_hg, calib48_hg, proj48_hg, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,jitter.factor=0, n.tau=NULL, ratio.max=2,ratio.max.trace=10*trace, ECBC=FALSE, ties='first',subsample=NULL, pp.type=7)
        }
        
        spei48_hgcalib <- qdm_calib48_hg$mhat.p
        
        
        ##### Complete to have the same length
        na_spei48_hgcalib[jlon,jlat,]    <- c(rep(NA, maxlength - length(spei48_hgcalib)),spei48_hgcalib)
        
        
      }
    }
    print(jlon)
  }
  
  
  #######################################
  ##################Save to netcdf
  ###################################
  
  time     <- ncvar_get(tmin_nc, "time")
  tunits   <- ncatt_get(tmin_nc,"time","units")
  calendar <- ncatt_get(tmin_nc,"time","calendar")$value
  
  xvals <- lons
  yvals <- lats
  
  nx <- length(xvals)
  ny <- length(yvals)
  xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
  ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
  tdim <- ncdim_def("time", tunits$value, as.numeric(time),unlim=T,calendar = calendar)
  
  mv <- NA    # missing value
  var_prec1 <- ncvar_def("spei48", 'spei48', list(xdim,ydim,tdim), mv )
  
  nt <- length(time ) # Imagine this is actually some very large number of timesteps
  
  ######Save to netcdf 
  ncout1 <- nc_create(paste(ImageDirectory_spei48,"spei4845_corrected_1850_2100.",as.character(nens),".nc", sep = ""), list(var_prec1),force_v4=TRUE)

  # put variables
  ncvar_put(ncout1, var_prec1, na_spei48_hgcalib, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  ncatt_put(ncout1,"time","axis","T")
  
  # get global attributes
  title        <- ncatt_get(tmin_nc,0,"SPEI48 Bias corrected")
  institution  <- ncatt_get(tmin_nc,0,"institution")
  datasource   <- ncatt_get(tmin_nc,0,"source")
  references   <- ncatt_get(tmin_nc,0,"references")
  history      <- ncatt_get(tmin_nc,0,"history")
  Conventions <- ncatt_get(tmin_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout1,0,"title",paste(" SPEI48: QDM from model ",as.character(nens), sep = ""))
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","ERA5")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  nc_close(prec_era_nc)
  nc_close(tmin_era_nc)
  nc_close(tmax_era_nc)
  nc_close(prec_nc)
  nc_close(tmax_nc)
  nc_close(tmin_nc)
  
  print(nens)
  
}

#}