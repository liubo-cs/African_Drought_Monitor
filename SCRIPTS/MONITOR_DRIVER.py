#################################################################
# Drought Monitor V2
# Written by: Nathaniel W. Chaney
# Date: 12 July 2013
# Location: Princeton University
# Contact: nchaney@princeton.edu
# Purpose: Main driver of the drought monitor
##################################################################

import MASTER_LIBRARY as ml
import datetime
import numpy as np
import grads
import os
import dateutil.relativedelta as relativedelta
import time

dataset_info = {
        'PGF':{'ctl':"../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl",'type':'xdfopen'},
        '3B42RT_BC':{'ctl':"../DATA/3B42RT_BC/DAILY/3B42RT_daily_0.25deg.ctl",'type':"xdfopen"},
        'GFS_ANL_BC':{'ctl':"../DATA/GFS_ANL_BC/DAILY/gfsanl_daily_0.25deg.ctl",'type':"xdfopen"},
        'VIC_PGF':{'ctl':"../DATA/VIC_PGF/DAILY/vic_daily_0.25deg.ctl",'type':"open"},
        'VIC_3B42RT':{'ctl':"../DATA/VIC_3B42RT/DAILY/vic_daily_0.25deg.ctl",'type':"open"},
        'ROUTING_VIC_PGF':{'ctl':"../DATA/ROUTING_VIC_PGF/DAILY/Streamflow.ctl",'type':"open"},
        'ROUTING_VIC_3B42RT':{'ctl':"../DATA/ROUTING_VIC_3B42RT/DAILY/Streamflow.ctl",'type':"open"},
        'MOD09_NDVI_MA':{'ctl':"../DATA/MOD09_NDVI_MA/DAILY/MOD09CMG_daily_0.25deg.ctl",'type':"xdfopen"},
        'MOD09_NDVI_MA_DERIVED':{'ctl':"../DATA/MOD09_NDVI_MA_DERIVED/DAILY/MOD09CMG_daily_0.25deg.ctl",'type':"xdfopen"},
        'VIC_DERIVED':{'ctl':"../DATA/VIC_DERIVED/DAILY/vic_derived_daily_0.25deg.ctl",'type':"xdfopen"},
        'ROUTING_VIC_DERIVED':{'ctl':"../DATA/ROUTING_VIC_DERIVED/DAILY/routing_vic_derived_daily_0.25deg.ctl",'type':"xdfopen"},
        'SPI':{'ctl':"../DATA/SPI/DAILY/spi_daily_0.25deg.ctl",'type':"xdfopen"},
        }

def Download_and_Process(date):

 #################################################
 #DOWNLOAD AND PREPROCESS ALL THE REQUIRED DATA
 #################################################

 ml.print_info_to_command_line("Downloading for %d/%d/%d" % (date.day,date.month,date.year))

 #Reprocess available pgf data (historical)
 #ml.Reprocess_PGF(date,dims)

 #Download and process the gfs final analysis data
 #ml.Download_and_Process_NCEP_FNL_Analysis(date,dims,idate,fdate,False)

 #Download and process the 3b42rt precipitation data
 ml.Download_and_Process_3b42RT(date,dims,False)

 #Download gfs forecast
 ml.Download_and_Process_GFS_forecast(date,dims,False)

 #Download and process modis NDVI
 ml.Download_and_Process_NDVI(date,dims,False)

 #Download and process the seasonal forecast
 ml.Download_and_Process_Seasonal_Forecast(date,False)

 #Download and process the gfs analysis data
 ml.Download_and_Process_GFS_Analysis(date,dims,False)

 return

def BiasCorrect(date):

 #################################################
 #BIAS CORRECT THE DOWNLOADED DATA
 #################################################

 ml.print_info_to_command_line("Bias Correcting for %d/%d/%d" % (date.day,date.month,date.year))

 #Regrid and downscale 3b42rt (ensure it complies with the pgf grid)
 ml.Regrid_and_Output_3B42rt(date,dims,False)

 #Bias correct the 3b42rt precipitation product
 ml.BiasCorrect_and_Output_Forcing_3B42RT_Daily(date,dims,False)

 #Bias correct the gfs final analysis product
 #ml.BiasCorrect_and_Output_Forcing_FNL_Daily(date,dims)

 #Bias correct the gfs forecast
 ml.BiasCorrect_and_Output_Forcing_GFS_Daily(date,dims,False)
 
 #Compute different moving averages of the ndvi product
 ml.Compute_NDVI_moving_average(date,dims,False)

 #Bias correct the gfs analysis product
 ml.BiasCorrect_and_Output_GFSANL_Daily(date,dims,False)

 return

def Compute_Indices(date):

 #################################################
 #COMPUTE INDICES
 #################################################

 Reprocess_Flag = False

 ml.print_info_to_command_line("Computing indices for %d/%d/%d" % (date.day,date.month,date.year))

 ml.Calculate_and_Output_SPI(date,dims,'monitor',date,Reprocess_Flag)

 ml.Calculate_and_Output_NDVI_Percentiles(date,dims,Reprocess_Flag)

 ml.Calculate_and_Output_SM_Percentiles(date,dims,'monitor',date,Reprocess_Flag)

 ml.Calculate_and_Output_Streamflow_Percentiles(date,dims,Reprocess_Flag,'monitor',date)

 return

def Compute_Forecast_Indices(date,idate):

 Reprocess_Flag = False

 ml.print_info_to_command_line("Computing forecast indices for %d/%d/%d" % (date.day,date.month,date.year))

 ml.Calculate_and_Output_Streamflow_Percentiles(date,dims,Reprocess_Flag,'forecast',idate)

 ml.Calculate_and_Output_SM_Percentiles(date,dims,'forecast',idate,Reprocess_Flag)

 ml.Calculate_and_Output_SPI(date,dims,'forecast',idate,Reprocess_Flag)
 
 ml.BiasCorrect_and_Compute_Seasonal_Forecast_Products(date,dims,False)

 return

def Compute_Averages(date,idate):

 #################################################
 #COLLECT DATASET BOUNDS
 #################################################

 if date == idate:
  for dataset in dataset_info:
   type = dataset_info[dataset]['type']
   ctl = dataset_info[dataset]['ctl']
   dataset_info[dataset] = ml.Collect_Dataset_Bounds(dataset,ctl,dataset_info[dataset],type)

 #################################################
 #COMPUTE MONTHLY AND ANNUAL PRODUCTS
 #################################################

 Averages_Reprocess_Flag = False#True
 for dataset in dataset_info:
  type = dataset_info[dataset]['type']
  ctl = dataset_info[dataset]['ctl']
  itime = dataset_info[dataset]['itime']
  ftime = dataset_info[dataset]['ftime']
  ml.Compute_Monthly_Yearly_Averages(date,dims,dt,dataset,ctl,type,Averages_Reprocess_Flag,itime,ftime)

 return

dims = {}
dims['minlat'] = -34.875000 #-89.8750
dims['minlon'] = -18.875000 #0.1250
dims['nlat'] = 292 #720
dims['nlon'] = 296 #1440
dims['res'] = 0.250
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)
dt = datetime.timedelta(days=1)
fdate = datetime.datetime.today() - 3*dt
#idate = datetime.datetime(fdate.year,fdate.month,1)
fdate = datetime.datetime(fdate.year,fdate.month,fdate.day)
idate = fdate - datetime.timedelta(days=20)
#idate = fdate
print idate
print fdate
date = datetime.datetime(2013,9,15)
#print date
ml.Download_and_Process_Seasonal_Forecast(date,True)
exit()
#idate = datetime.datetime(2003,1,31)
#fdate = datetime.datetime(2003,12,31)
#Download and bias correct data
date = idate
while date <= fdate:
 Download_and_Process(date)
 BiasCorrect(date)
 date = date + dt

#Run VIC
idate_model = datetime.datetime(fdate.year,fdate.month,1)
ml.Run_VIC(idate_model,fdate,dims,'3b42rt',False)
#ml.Run_VIC(idate,fdate,dims,'pgf')
ml.Run_VIC(idate_model,fdate,dims,'gfs_forecast',False)

#Run the Routing model
#ml.Run_VDSC(idate,fdate+datetime.timedelta(days=7),dims,'pgf')
ml.Run_VDSC(idate_model,fdate,dims,'3b42rt')
ml.Run_VDSC(idate_model,fdate,dims,'gfs_forecast')

#Compute Monitor Indices
date = idate
while date <= fdate:
 #print date
 Compute_Indices(date) 
 date = date + dt

#Compute Forecast Indices
date = fdate + datetime.timedelta(days=1)
while date <= fdate + 7*dt:
 Compute_Forecast_Indices(date,fdate + dt)
 date = date + dt

#Compute Averages
date = idate
while date <= fdate:
 print date
 Compute_Averages(date,idate)
 date = date + dt

#Finalize GFS forecast
ml.Finalize_GFS_forecast(fdate+datetime.timedelta(days=1),dims)

#ml.Finalize_NCEP_FNL_Connection()
