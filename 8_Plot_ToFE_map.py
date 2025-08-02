#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 20:26:10 2025

@author: vecchia
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cmocean
import matplotlib.patches as mpatches
import os
import scipy
import h5netcdf
import netCDF4
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cartopy.feature as cfeature
from shapely.geometry import box, Polygon, MultiPolygon
from matplotlib import ticker, cm 
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import rcParams

os.chdir('/data/')


ds_ToFE                          = xr.open_dataset('ToFE/ToFE_FAR_spei4815_srfi4815_swsi4806_nc')
ds_reservoirs                    = xr.open_dataset('../04_04_2025/reservoir_year_cap.nc')

data_ToFE                       = ds_ToFE['ToFE']
data_reservoirs                 = ds_reservoirs['year_reser']

#################################################################################################################
############################ ToFE od DZD with CESM2-LE SSP370 ###################################################
#################################################################################################################

#-- create figure and axes object
fig = plt.figure(figsize=(20,20))

#-- choose map projection
ax = plt.axes(projection=ccrs.Robinson())
ax.set_extent([-180, 180, -60, 90])

#-- add coastlines, country border lines, and grid lines
ax.coastlines(zorder=1)


#-- create states outlines
states_provinces = cf.NaturalEarthFeature(category='cultural',
                                        name='admin_1_states_provinces_lines',
                                        scale='50m',
                                        facecolor='red')
                                        #facecolor='gray')

ax.add_feature(cf.BORDERS, linewidth=0.8, edgecolor='black', zorder=2)
ax.add_feature(cf.OCEAN.with_scale('50m'), facecolor='white')
ax.add_feature(cf.LAND.with_scale('50m'), facecolor='gainsboro', zorder=0)

levels= 1900,1910,1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100

#cmap= (mpl.colors.ListedColormap(['#FFB266','red','#660000']).with_extremes(under='yellow', over='#191970'))

#newcmp = mpl.colors.ListedColormap(['#000ecd', '#008dff', '#00e0ff', '#00fdd1', '#00fb73', '#13fb00', '#60ce00', '#7af200', '#a9ff2c', '#dbff20',
#      '#ffee00', '#ffd303', '#ffaa0c', '#ff5002', '#ff1800', '#ff009f', '#d615ff', '#b64af9', '#ed89ef', '#f5c0f7'])

newcmp = mpl.colors.ListedColormap(['#000ecd', '#008dff', '#00e0ff', '#00fb73', 'yellowgreen', '#60ce00', '#7af200', '#a9ff2c', '#dbff20',
      'yellow','#ffee00', '#ffd303', '#ffaa0c', 'crimson', '#ff1800', '#ff009f', '#d615ff', '#b64af9', '#ed89ef', '#f5c0f7'])
norm = mpl.colors.BoundaryNorm(levels, newcmp.N)



cnplot = ax.pcolormesh(ds_ToFE.lon, ds_ToFE.lat, data_ToFE,shading='auto', #0 means first and only time step
                                cmap=newcmp, #'plasma', #'gist_ncar', #'YlOrRd', #'coolwarm', #'seismic', #'viridis', #'YlOrRd', #'cmo.balance', #'Spectral_r',
                                norm =norm, zorder=0, transform=ccrs.PlateCarree())

cnp = ax.contourf(ds_reservoirs.lon, ds_reservoirs.lat,data_reservoirs.where(data_ToFE > 0),
                                hatches=['****'], 
                                levels=levels, 
                                colors = 'none',
                                transform=ccrs.PlateCarree())


# Create new axes according to image position
cax = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.02, ax.get_position().width, 0.02])

#-- add colorbar
cbar = plt.colorbar(cnplot,orientation="horizontal",ticks=levels, cax=cax)
cbar.ax.tick_params(labelsize=25)
cbar.ax.set_title('Time [decade]', fontweight='bold', size=20, pad=17) 
cbar.ax.set_xticklabels(labels = levels ,fontweight='bold', size=22)  # vertically oriented colorbar

#-- save graphic output to PNG file
plt.savefig('ToFE_DZD.png',bbox_inches='tight',dpi=300)

