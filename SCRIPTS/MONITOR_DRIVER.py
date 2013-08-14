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

 #Download gfs analysis data
 #ml.Download_and_Process_GFS_Historical(date,dims)

 #Download and process the 3b42rt precipitation data
 #ml.Download_and_Process_3b42RT(date,dims,'standard')

 #Download gfs forecast
 #ml.Download_and_Process_GFS_forecast(date,dims)

 #Download and process modis NDVI
 #ml.Download_and_Process_NDVI(date,dims)

 #Download and process the seasonal forecast
 #ml.Download_and_Process_Seasonal_Forecast(date)

 #################################################
 #BIAS CORRECT THE DOWNLOADED DATA
 #################################################

 #Regrid and downscale 3b42rt (ensure it complies with the pgf grid)
 #ml.Regrid_and_Output_3B42rt(date,dims)

 #Bias correct the 3b42rt precipitation product
 #ml.BiasCorrect_and_Output_Forcing_3B42RT_Daily(date,dims)

 #Bias correct the gfs final analysis product
 ml.BiasCorrect_and_Output_Forcing_FNL_Daily(date,dims)

 #Bias correct the gfs forecast
 #ml.BiasCorrect_and_Output_Forcing_GFS_Forecast

 #Bias correct the seasonal forecast
 
 #Compute different moving averages of the ndvi product
 #ml.Compute_NDVI_moving_average(date,dims)

 #################################################
 #COMPUTE MONTHLY AND ANNUAL PRODUCTS
 #################################################

 #ml.Compute_Averages_3b42RT_BC(date,dims,dt)

 #ml.Compute_Averages_PGF(date,dims,dt,'standard')

 #################################################
 #COMPUTE INDICES
 #################################################

 #ml.Calculate_and_Output_SPI(date,dims)

 #ml.Calculate_and_Output_NDVI_Percentiles(date,dims)

 #ml.Calculate_and_Output_SM_Percentiles(date,dims)

 return

#1. Determine the period that needs to be updated

#2. Download all the relevant data

#3. Process all the relevant data

#4. Run all the relevant models

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
idate = datetime.datetime(2010,1,1)
fdate = datetime.datetime(2012,12,31)
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

#ml.Run_VIC(forcing_file,idate,fdate)

#Prepare the VIC global parameter file
#ml.Prepare_VIC_Global_Parameter_File(idate,fdate,dims)

#Prepare the VIC forcings
#forcing_file = ml.Prepare_VIC_Forcings_Historical(idate,fdate,dims)

#Run VIC
#ml.Run_VIC(idate,fdate,dims)


 
#################################################
#RUN ROUTING MODEL
#################################################
'''
#################################################
#FIT DISTR
#################################################

#PGF dataset
#ctl_in = "../DATA/PGF/MONTHLY/pgf_monthly_0.25deg.ctl"
#dt = relativedelta.relativedelta(years=1)
#idate = datetime.datetime(1948,1,1)
#fdate = datetime.datetime(2008,12,31)
#dt_down = datetime.timedelta(days=0)
#dt_up = relativedelta.relativedelta(days=0)
#var = "prec"
#type = "all"
##Extract the desired data
#data = ml.Extract_Data_Period_Average(idate,fdate,dt_down,dt_up,dt,ctl_in,var,type)
#Calculate percentiles
#pct = ml.Calculate_Percentiles(data)
#dt = datetime.timedelta(days=1)
#date = datetime.datetime(1950,1,1)
#fdate = datetime.datetime(2008,12,31)
#while date <= fdate:
# ml.Calculate_and_Output_SPI(date,dims,'rp')
# date = date + dt
#Calculate the distribution parameters
#file_out = "../WORKSPACE/parameters.nc"
#ml.Fit_Distribution('gamma',data,file_out,dims)

#################################################
#BIAS CORRECT AND PREPARE ALL THE REQUIRED DATA
#################################################

#CLIENT SECTION

#################################################
#RUN ALL THE MODELS
#################################################

#################################################
#COMPUTE INDICES
#################################################

#Fit a gamma distribution to the 3B42RT product
ctl_in = "../DATA/3B42RT/MONTHLY/3B42RT_monthly_0.25deg.ctl"
file_out = "../WORKSPACE/parameters.nc"
dt = datetime.timedelta(days = 31)
ml.Fit_Distribution('gamma',ctl_in,idate,fdate,dt,'prec',file_out,dims)

date = datetime.datetime(2013,1,1)
#Download and process the seasonal forecast
ml.Download_and_Process_Seasonal_Forecast(date)

#Download and process modis NDVI
date = idate
dt = datetime.timedelta(days=1)
while date <= fdate:
 print date
 ml.Download_and_Process_NDVI(date,dims)
 date = date + dt

idate = datetime.datetime(2011,1,1)
fdate = datetime.datetime(2011,12,31)
date= idate
dt_ma = [datetime.timedelta(days=60),datetime.timedelta(days=30),datetime.timedelta(days=20),datetime.timedelta(days=10),datetime.timedelta(days=5),datetime.timedelta(days=1)]
dt = datetime.timedelta(days=1)
while date <= fdate:
 print date
 ml.Compute_NDVI_moving_average(date,dt_ma,dims)
 date = date + dt


while date <= fdate:
 #End of month routines
 ndate = date + dt
 print date,fdate
 if date.month != ndate.month:
  print "Tomorrow is a new month. We have some extra work to do today."
  idate_in = datetime.datetime(date.year,date.month,1)
  fdate_in = date
  #3B42RT
  ctl_in = "../DATA/3B42RT/3B42RT_daily_0.25deg.ctl"
  file_out = "../DATA/3B42RT/MONTHLY/3B42RT_%04d%02d_daily_0.250deg.nc" % (idate_in.year,idate_in.month)
  Compute_and_Output_Averages(ctl_in,file_out,idate_in,fdate_in,dims)
 date = date + dt

 #FNL
 #ctl_in = "../DATA/FNL_ANALYSIS/fnlanl_daily_0.25deg.ctl"
 #file_out = "../DATA/FNL_ANALYSIS/MONTHLY/fnlanl_%04d%02d_daily_0.250deg.nc" % (idate.year,idate.month)
 #Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims)

#End of year routines
ndate = date + dt
if date.year != ndate.year:
 print "Tomorrow is a new year. We have a lot of extra work to do today."
 idate = datetime.datetime(date.year,1,1)
 fdate = date
 #3B42RT
 ctl_in = "../DATA/3B42RT/3B42RT_daily_0.25deg.ctl"
 file_out = "../DATA/3B42RT/YEARLY/3B42RT_%04d_daily_0.250deg.nc" % (idate.year)
 Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims)
 #FNL
 ctl_in = "../DATA/FNL_ANALYSIS/fnlanl_daily_0.25deg.ctl"
 file_out = "../DATA/FNL_ANALYSIS/YEARLY/fnlanl_%04d_daily_0.250deg.nc" % (idate.year)
 Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims)
'''
