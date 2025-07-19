rm(list=ls())

library(ncdf4)
library(stats)
# You can install aomisc from GitHub
#install.packages("devtools")
#devtools::install_github("onofriAndreaPG/aomisc")
### https://www.statforbiology.com/nonlinearregression/usefulequations
#or
#install.packages("remotes")
#remotes::install_github("OnofriAndreaPG/aomisc")
library(drc)
library(nlme)
library(aomisc)
#library(data.table)



model <- c('gfdl','hadgem','ipsl','miroc','noresm') #https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VIQEAB


for (ens in 1:length(model)){

  ####################################################################################
  ####################################################################################
  
  Twat_model_nc          = nc_open(paste(as.character(model[ens]),'/total_consumption_ssp2_rcp45_',as.character(model[ens]),'_sectors_monthly.nc',sep=''))
  Twat_model_grid        = ncvar_get(Twat_model_nc,  varid='water')

  
  lats <- ncvar_get(Twat_model_nc,'lat')
  lons <- ncvar_get(Twat_model_nc,'lon')
  
  nlats <- length(lats)
  nlons <- length(lons)

  ssp2_rcp45_waterdemand_past              <- array(NA,c(nlons,nlats,1920))     ### 1850 -2009
  
  for (jlon in 1:nlons){
    for (jlat in 1:nlats){
      mean_Twat    <- mean(Twat_model_grid[jlon,jlat,],na.rm=T)
      
      ##############################
      Twat_model  <- Twat_model_grid[jlon,jlat,]
      #val_Twat_all    <- Twat_model[13:1092]*1e12
      val_Twat_all    <- Twat_model*1e12  ##### km3 to Liter
      val_Twat_all[is.infinite(val_Twat_all)] <- 0
      val_Twat_all[is.nan(val_Twat_all)] <- 0
      
      #length(val_Twat_all)
      if ( !is.na(mean_Twat) && mean_Twat > 0){
        
        val_Twat <- val_Twat_all
        ################################################ reverse the data to be 2100 to 1850 to project the value of the past using exponantial decay function ########################################

        month <- rep_len(1:12, length(val_Twat))
        
        df_data <- data.frame(x=rev(month), y=rev(val_Twat))
        head(df_data)
        length(month)
        length(val_Twat)
        
        wt_jan0 <- df_data$y[df_data$x==1]
        
        wt_feb0 <- df_data$y[df_data$x==2]
        
        wt_mar0 <- df_data$y[df_data$x==3]
        
        wt_apr0 <- df_data$y[df_data$x==4]
        
        wt_mai0 <- df_data$y[df_data$x==5]
        
        wt_jun0 <- df_data$y[df_data$x==6]
        
        wt_jul0 <- df_data$y[df_data$x==7]
        
        wt_aug0 <- df_data$y[df_data$x==8]
        
        wt_sep0 <- df_data$y[df_data$x==9]
        
        wt_oct0 <- df_data$y[df_data$x==10]
        
        wt_nov0 <- df_data$y[df_data$x==11]
        
        wt_dec0 <- df_data$y[df_data$x==12]
        
        ######################################### x for each month
        month0 <- seq(1,91,1)
        month1 <- seq(92,252,1)
        length(month)
        ##################### Projection
        dati_jan <- data.frame(x=month0 , y=wt_jan0)
        model_jan <- drm(y ~ x, fct = DRC.logCurve(),data = dati_jan)
        pred_jan <- predict(model_jan , data.frame(month1), se = FALSE)
        real_past_predict_jan <- rev(pred_jan)
        
        dati_feb <- data.frame(x=month0 , y=wt_feb0)
        model_feb <- drm(y ~ x, fct = DRC.logCurve(),data = dati_feb)
        pred_feb <- predict(model_feb , data.frame(month1), se = FALSE)
        real_past_predict_feb <- rev(pred_feb)
        
        dati_mar <- data.frame(x=month0 , y=wt_mar0)
        model_mar <- drm(y ~ x, fct = DRC.logCurve(),data = dati_mar)
        pred_mar <- predict(model_mar , data.frame(month1), se = FALSE)
        real_past_predict_mar <- rev(pred_mar)
        
        dati_apr <- data.frame(x=month0 , y=wt_apr0)
        model_apr <- drm(y ~ x, fct = DRC.logCurve(),data = dati_apr)
        pred_apr <- predict(model_apr , data.frame(month1), se = FALSE)
        real_past_predict_apr <- rev(pred_apr)
        
        dati_mai <- data.frame(x=month0 , y=wt_mai0)
        model_mai <- drm(y ~ x, fct = DRC.logCurve(),data = dati_mai)
        pred_mai <- predict(model_mai , data.frame(month1), se = FALSE)
        real_past_predict_mai <- rev(pred_mai)
        
        dati_jun <- data.frame(x=month0 , y=wt_jun0)
        model_jun <- drm(y ~ x, fct = DRC.logCurve(),data = dati_jun)
        pred_jun <- predict(model_jun , data.frame(month1), se = FALSE)
        real_past_predict_jun <- rev(pred_jun)
        
        dati_jul <- data.frame(x=month0 , y=wt_jul0)
        model_jul <- drm(y ~ x, fct = DRC.logCurve(),data = dati_jul)
        pred_jul <- predict(model_jul , data.frame(month1), se = FALSE)
        real_past_predict_jul <- rev(pred_jul)
        
        dati_aug <- data.frame(x=month0 , y=wt_aug0)
        model_aug <- drm(y ~ x, fct = DRC.logCurve(),data = dati_aug)
        pred_aug <- predict(model_aug , data.frame(month1), se = FALSE)
        real_past_predict_aug <- rev(pred_aug)
        
        dati_sep <- data.frame(x=month0 , y=wt_sep0)
        model_sep <- drm(y ~ x, fct = DRC.logCurve(),data = dati_sep)
        pred_sep <- predict(model_sep , data.frame(month1), se = FALSE)
        real_past_predict_sep <- rev(pred_sep)
        
        dati_oct <- data.frame(x=month0 , y=wt_oct0)
        model_oct <- drm(y ~ x, fct = DRC.logCurve(),data = dati_oct)
        pred_oct <- predict(model_oct , data.frame(month1), se = FALSE)
        real_past_predict_oct <- rev(pred_oct)
        
        dati_nov <- data.frame(x=month0 , y=wt_nov0)
        model_nov <- drm(y ~ x, fct = DRC.logCurve(),data = dati_nov)
        pred_nov <- predict(model_nov , data.frame(month1), se = FALSE)
        real_past_predict_nov <- rev(pred_nov)
        
        dati_dec <- data.frame(x=month0 , y=wt_dec0)
        model_dec <- drm(y ~ x, fct = DRC.logCurve(),data = dati_dec)
        pred_dec <- predict(model_dec , data.frame(month1), se = FALSE)
        real_past_predict_dec <- rev(pred_dec)
        
        past_pred_total_wat_per_mon <- list()
        year=0
        #for(year in 1:91){
        for(num in seq(0,1932,12)){
          year= year + length(num)
          past_pred_total_wat_per_mon[num+1] <- real_past_predict_jan[year] 
          past_pred_total_wat_per_mon[num+2] <- real_past_predict_feb[year]
          past_pred_total_wat_per_mon[num+3] <- real_past_predict_mar[year] 
          past_pred_total_wat_per_mon[num+4] <- real_past_predict_apr[year] 
          past_pred_total_wat_per_mon[num+5] <- real_past_predict_mai[year] 
          past_pred_total_wat_per_mon[num+6] <- real_past_predict_jun[year] 
          past_pred_total_wat_per_mon[num+7] <- real_past_predict_jul[year] 
          past_pred_total_wat_per_mon[num+8] <- real_past_predict_aug[year] 
          past_pred_total_wat_per_mon[num+9] <- real_past_predict_sep[year] 
          past_pred_total_wat_per_mon[num+10] <- real_past_predict_oct[year] 
          past_pred_total_wat_per_mon[num+11] <- real_past_predict_nov[year] 
          past_pred_total_wat_per_mon[num+12] <- real_past_predict_dec[year] 
          #print(year)
        }
        #}
        
        past_pred_total_wat_per_mon2 <- unlist(past_pred_total_wat_per_mon)
        length(past_pred_total_wat_per_mon2)

        past_pred_total_wat_per_mon3 <- past_pred_total_wat_per_mon2[1:1920] ##### 1850 -2009
        
        
        past_pred_total_wat_per_mon3[is.infinite(past_pred_total_wat_per_mon3)] <- NA
        past_pred_total_wat_per_mon3[is.nan(past_pred_total_wat_per_mon3)] <- NA
        past_pred_total_wat_per_mon4 <- ifelse(past_pred_total_wat_per_mon3 < 0, 0, past_pred_total_wat_per_mon3)
        
        ssp2_rcp45_waterdemand_past[jlon,jlat,] <-  past_pred_total_wat_per_mon4 *1e-12 ##### liter to km3
        
        
      }
      
    }
    #print(jlon)
  }

  
  pop_past_nc          <- nc_open(paste("../ssp360/total_waterDeman_models/water_consumption/past_total_water_consumption_ssp3_rcp60_miroc_sectors_monthly_1850_2009_Log_time_rev.nc",sep=""))
  time                 <- ncvar_get(pop_past_nc, "time")
  tunits               <- ncatt_get(pop_past_nc,"time","units")
  calendar             <- ncatt_get(pop_past_nc,"time","calendar")$value
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
  var_prec1 <- ncvar_def("water", 'water', list(xdim,ydim,tdim), mv)
  
  ncout1 <- nc_create(paste0(as.character(model[ens]),'/past_total_water_demand_ssp2_rcp45_',as.character(model[ens]),'_sectors_monthly_1850_2010_Log_liter_time_rev.nc',sep=''), list(var_prec1),force_v4=TRUE)
  #    ncvar_put(ncout1, var_prec1, ssp2_rcp45_waterdemand_past, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  ncvar_put(ncout1, var_prec1, ssp2_rcp45_waterdemand_past, start=c(1,1,1), count=c(nx,ny,nt), verbose=TRUE)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout1,"lon","axis","X") 
  ncatt_put(ncout1,"lat","axis","Y")
  ncatt_put(ncout1,"time","axis","T")
  
  # get global attributes
  title        <- ncatt_get(Twat_model_nc,0,"reservoirs")
  institution  <- ncatt_get(Twat_model_nc,0,"institution")
  datasource   <- ncatt_get(Twat_model_nc,0,"source")
  references   <- ncatt_get(Twat_model_nc,0,"references")
  history      <- ncatt_get(Twat_model_nc,0,"history")
  Conventions  <- ncatt_get(Twat_model_nc,0,"Conventions")
  
  # add global attributes
  ncatt_put(ncout1,0,"title","Water demand in km3")
  ncatt_put(ncout1,0,"institution","ICCP & PNU, Busan")
  ncatt_put(ncout1,0,"source",datasource$value)
  ncatt_put(ncout1,0,"references","https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VIQEAB")
  history <- paste("Ravinandrasana.P.V.", date(), sep=", ")
  ncatt_put(ncout1,0,"history",history)
  ncatt_put(ncout1,0,"Conventions",Conventions$value)
  nc_close(ncout1)
  
  nc_close(Twat_model_nc)
  nc_close(pop_past_nc)
  
}


