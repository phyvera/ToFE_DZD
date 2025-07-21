rm(list=ls())

library(ncdf4)
library(zoo)
library(MBC)
library(SPEI)
library(dplyr)
library(tidyr)
library(stats)

setwd('./')

ImageDirectory ='SPEI_SRFI_SWSI/swsi48/'

model_watco <- rep_len(c('gfdl','hadgem','ipsl','miroc','noresm'), 100) ### to distribute along CESM memebers 


imod <- 0

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
  imod <- imod + 1
  
  ###### atmospheric variables
  prec_nc <- nc_open(paste("precipitation/pr_monsum",as.character(nens),".nc",sep=""))
  tmax_nc <- nc_open(paste("Tmax/tmax_mon",as.character(nens),".nc",sep=""))
  tmin_nc <- nc_open(paste("Tmin/tmin_mon",as.character(nens),".nc",sep=""))
  
  
  prec <- ncvar_get(prec_nc,varid="pr")
  tmin <- ncvar_get(tmin_nc,varid="tmin")
  tmax <- ncvar_get(tmax_nc,varid="tmax")
  
  lats <- ncvar_get(prec_cmip_nc,'lat')
  lons <- ncvar_get(prec_cmip_nc,'lon')
  
  nlats <- length(lats)
  nlons <- length(lons)
  
  ##### river variable
  river_nc   <- nc_open(paste("River/river_monmean",as.character(nens),".nc",sep=""))
  river      <- ncvar_get(river_nc,varid='river')
  
  lats_river <- ncvar_get(river_nc,'lat')
  lons_river <- ncvar_get(river_nc,'lon')

  
  ###### 
  Twat_model_nc      <- nc_open(paste('waterDemand/water_consumption/',as.character(model[imod]),'/total_consumption_',as.character(model[imod]),'_sectors_monthly_1850_2100.nc',sep=''))
  Twat_model_grid    <- ncvar_get(Twat_model_nc,  varid='water')
  
  lats_Twat_model     <- ncvar_get(Twat_model_nc,'lat')
  lons_Twat_model     <- ncvar_get(Twat_model_nc,'lon')
  

  ###### ERA realaysis
  prec_era_nc     <- nc_open("ERA/monthly_total_precipitation_era5_sum_remapcon.nc")
  tmin_era_nc     <- nc_open("ERA/monTmin.era5_remapcon.nc")
  tmax_era_nc     <- nc_open("ERA/monTmax.era5_remapcon.nc")
  river_era_nc    <- nc_open("ERA/monRiverEra_GloFAS_mosartgrid_1979_2022.nc")
  
  prec_era        <- ncvar_get(prec_era_nc,varid="tp")
  tmin_era        <- ncvar_get(tmin_era_nc,varid="mn2t")
  tmax_era        <- ncvar_get(tmax_era_nc,varid="mx2t")
  river_era       <- ncvar_get(river_era_nc,varid='dis')
  

  ########## create array
  swsi48          <- array(NA,c(nlons,nlats,3012))
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      land_riv_era <-  mean(remap_river_era[jlon,jlat,], na.rm=TRUE)
      
      if (!is.na(land_riv_era) && land_riv_era > 0 ){
        
        ### projections    
        ts_river              <- river[jlon,jlat,]         ### in m3/s
        ts_Twat_model         <- Twat_model_grid[jlon,jlat,]    ### in km3
        ts_pet_hgcmip         <- hargreaves(Tmin = (tmin_cmip[jlon,jlat,] - 273.15), Tmax = (tmax_cmip[jlon,jlat,]-273.15),lat = lats[jlat],na.rm =TRUE)
        ts_PE_hgcmip          <- (prec_cmip[jlon,jlat,] * 86400) - ts_pet_hgcmip
        
        vol_waterBalan        <-  ts_PE_hgcmip * 0.001               #### m/m2 = m3 /month           
        river_m3              <-  (ts_river * 30*24*60*60)        ### in m3/month
        #river_m3             <-  (river *24*60*60)
        Twat_m3               <-  ts_Twat_model * 1e9              #### m3
        
        #### ERA5 
        ts_river_era           <- river_era[jlon,jlat,]         ### in m3/s
        pet_hgera              <- hargreaves(Tmin = tmin_era[jlon,jlat,], Tmax = tmax_era[jlon,jlat,],lat = lats[jlat],na.rm =TRUE)
        PE_hgera               <- prec_era[jlon,jlat,]  - pet_hgera
        
        vol_waterBalan_era     <-  PE_hgera * 0.001               #### m/m2 = m3 /month           
        river_era_m3           <-  (ts_river_era * 30*24*60*60)        ### in m3/month
        
        ##############  SWSI48 
        rol_watsup48           <- rollapply((vol_waterBalan + river_m3),48,sum,align = "right",fill = NA)
        rol_watsup48_era       <- rollapply((vol_waterBalan_era+ river_era_m3),48,sum,align = "right",fill = NA)
        
        #QDM water supply
        in_watsuply        <- ts(rol_watsup48,start=c(1850,1),frequency=12)  
        in_calib_watsuply  <- window(in_watsuply,start=c(1979,1),end=c(2022,12),frequency=12)
        
        if (sum(!is.na(rol_watsup48))>4 && sum(!is.na(rol_watsup48_era))>4 && sum(!is.na(in_calib_watsuply))>4){
          proj48_watsup       <- na.approx(rol_watsup48)
          obs48_watsup        <- na.approx(rol_watsup48_era)
          calib48_watsup      <- na.approx(in_calib_watsuply)
          qdm_calib48_watsup  <- QDM(obs48_watsup, calib48_watsup, proj48_watsup, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,jitter.factor=0, n.tau=NULL, ratio.max=2,ratio.max.trace=10*trace, ECBC=FALSE, ties='first',subsample=NULL, pp.type=7)
          watsuply48_calib    <- qdm_calib48_watsup$mhat.p
        }

        
        ##### Complete to have the same length
        watsuply48_corrected    <- c(rep(NA, maxlength - length(watsuply48_calib)),watsuply48_calib)
        
        rol_watdemand48       <- rollapply(Twat_m3,48,sum,align = "right",fill = NA)
        water_inx48           <- watsuply48_corrected / rol_watdemand48
        
        water_inx48[is.infinite(water_inx48)] <- NA
        water_inx48[is.nan(water_inx48)] <- NA

        swsi48[jlon,jlat,]     <- water_inx48
        
      }
    }
    #print(jlon)
  }
 
  #######################################
  ##################Save to netcdf 
  ###################################
   
  time                 <- ncvar_get(prec_cmip_nc, "time")
  tunits               <- ncatt_get(prec_cmip_nc,"time","units")
  calendar             <- ncatt_get(prec_cmip_nc,"time","calendar")$value
  ntime                <- dim(time)
  
  xvals <- lons
  yvals <- lats
  
  nx   <- length(xvals)
  ny   <- length(yvals)
  xdim <- ncdim_def( 'lon', 'degrees_east', xvals)
  ydim <- ncdim_def( 'lat', 'degrees_north', yvals)
  tdim <- ncdim_def("time", tunits$value, as.numeric(time),unlim=T,calendar = calendar)
  nt <- length(time) # Imagine this is actually some very large number of timesteps
  
  
  mv <- NA     # missing value
  var_prec1 <- ncvar_def("swsi48", 'swsi48', list(xdim,ydim,tdim), mv)
  
  ncout1 <- nc_create(paste0(ImageDirectory,'swsi48_correctedInput_1850_2100.',as.character(nens),'.nc',sep=''), list(var_prec1),force_v4=TRUE)
  ncvar_put(ncout1, var_prec1, swsi48, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  ## put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  ncatt_put(ncout1,"time","axis","T")
  
  ## get global attributes
  title        <- ncatt_get(prec_cmip_nc,0,"SWSI48 bias corrected")
  institution  <- ncatt_get(prec_cmip_nc,0,"institution")
  datasource   <- ncatt_get(prec_cmip_nc,0,"source")
  references   <- ncatt_get(prec_cmip_nc,0,"references")
  history      <- ncatt_get(prec_cmip_nc,0,"history")
  Conventions  <- ncatt_get(prec_cmip_nc,0,"Conventions")

  ncatt_put(ncout1,0,"title",paste("Bias corrected SWSI48: QDM from nodel",as.character(nens), sep = ""))
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VIQEAB")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  nc_close(prec_era_nc)
  nc_close(tmin_era_nc)
  nc_close(tmax_era_nc)
  nc_close(river_era_nc)
  
  nc_close(prec_cmip_nc)
  nc_close(tmax_cmip_nc)
  nc_close(tmin_cmip_nc)
  nc_close(river_nc)
  
  print(nens) 
}


