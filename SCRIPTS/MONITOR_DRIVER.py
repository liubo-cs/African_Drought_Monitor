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

def Download_and_BiasCorrect(date):

 print date
 #################################################
 #DOWNLOAD AND PREPROCESS ALL THE REQUIRED DATA
 #################################################

 #Reprocess available pgf data (historical)
 #ml.Reprocess_PGF(date,dims)

 #Download and process the gfs final analysis data
 #ml.Download_and_Process_NCEP_FNL_Analysis(date,dims,idate,fdate)

 #Download and process the 3b42rt precipitation data
 #ml.Download_and_Process_3b42RT(date,dims,'standard')

 #Download gfs forecast
 #ml.Download_and_Process_GFS_forecast(date,dims)

 #Download and process modis NDVI
 #ml.Download_and_Process_NDVI(date,dims)

 #Download and process the seasonal forecast
 #ml.Download_and_Process_Seasonal_Forecast(date)

 #Download and process the gfs analysis data
 #ml.Download_and_Process_GFS_Analysis(date,dims)

 #################################################
 #BIAS CORRECT THE DOWNLOADED DATA
 #################################################

 #Regrid and downscale 3b42rt (ensure it complies with the pgf grid)
 #ml.Regrid_and_Output_3B42rt(date,dims)

 #Bias correct the 3b42rt precipitation product
 #ml.BiasCorrect_and_Output_Forcing_3B42RT_Daily(date,dims)

 #Bias correct the gfs final analysis product
 #ml.BiasCorrect_and_Output_Forcing_FNL_Daily(date,dims)

 #Bias correct the gfs forecast
 #ml.BiasCorrect_and_Output_Forcing_GFS_Daily(date,dims)

 #Bias correct the seasonal forecast
 
 #Compute different moving averages of the ndvi product
 #ml.Compute_NDVI_moving_average(date,dims)

 #Bias correct the gfs analysis product
 #ml.BiasCorrect_and_Output_GFSANL_Daily(date,dims)

 #################################################
 #COMPUTE INDICES
 #################################################

 #ml.Calculate_and_Output_SPI(date,dims)

 #ml.Calculate_and_Output_NDVI_Percentiles(date,dims)

 #ml.Calculate_and_Output_SM_Percentiles(date,dims)

 #################################################
 #COMPUTE MONTHLY AND ANNUAL PRODUCTS
 #################################################

 #ml.Compute_Averages_3b42RT_BC(date,dims,dt)

 #ml.Compute_Averages_PGF(date,dims,dt,'standard')

 #ml.Compute_Averages_SM_Percentiles(date,dims,dt)

 #ml.Compute_Averages_SPI(date,dims,dt)


 #1. Determine the period that needs to be updated

 #2. Download all the relevant data
 
 #3. Process all the relevant data

 #4. Run all the relevant models

 #Create the runoff files for routing
 #ml.Extract_VIC_Baseflow_and_Runoff(date,dims)

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
date = datetime.datetime.today()
idate = datetime.datetime(date.year,date.month,date.day) - 6*dt
idate = datetime.datetime(2003,1,1)
fdate = datetime.datetime(2006,12,31)
date = idate
dates = []
#while date <= fdate:
# dates.append(date)
# date = date + dt

#SERVER SECTION
tic = time.clock()
while date <= fdate:

 Download_and_BiasCorrect(date)
 date = date + dt
toc = time.clock()
print toc - tic

#ml.Finalize_NCEP_FNL_Connection()

#################################################
#RUN THE VIC MODEL
#################################################

#Run the model

#Run VIC
#ml.Run_VIC(idate,fdate,dims,'3b42rt')
#ml.Run_VIC(idate,fdate,dims,'pgf')
#ml.Run_VIC(idate,fdate,dims,'gfsanl')

 
#################################################
#RUN ROUTING MODEL
#################################################

#Run the VDSC model (Josh Roundy)

ml.Run_VDSC(idate,fdate,dims,'3b42rt')
