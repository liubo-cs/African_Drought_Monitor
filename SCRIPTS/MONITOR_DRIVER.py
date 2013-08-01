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

'''
def Download_and_Process_Seasonal_Forecast(date):

 models = ('CMC1-CanCM3','CMC2-CanCM4','COLA-RSMAS-CCSM3','GFDL-CM2p1-aer04','NASA-GMAO-062012')
 nensembles = (10,10,6,10,12)
 type = ('.FORECAST/','.FORECAST/','','','')
 vars = ('prec','tref')
 month = date.strftime('%b') 
 year = date.year
 dir0 = '../DATA/SEASONAL_FORECAST/%04d%02d' % (date.year,date.month)
 if not os.path.exists(dir0):
  os.makedirs(dir0)
 for var in vars:
  dir = dir0 + '/' + var
  if not os.path.exists(dir):
   os.makedirs(dir)
  i = 0 
  for model in models:
   for ens in xrange(1,nensembles[i]+1):
    root = ''
    http_file = 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.{0}/{3}.MONTHLY/.{5}/S/%280000%201%20{1}%20{2}%29%280000%201%20{1}%20{2}%29RANGEEDGES/M/%28{4}.%29VALUES/M/removeGRID/S/removeGRID/-9.9900000E08/setmissing_value/data.cdf'.format(model,month,year,type[i],ens,var)
    file_out = '{0}/{1}_{2}.nc'.format(dir,model,ens)
    os.system('wget -O {0} {1}'.format(file_out,http_file))
   i = i + 1
'''
 
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
'''
date = datetime.datetime.today()
date = datetime.datetime(date.year,date.month,date.day) - 6*dt

#Initialize connection for FNL analysis
connection_info = ml.Initialize_NCEP_FNL_Connection()

#Download and process the gfs analysis data
ml.Download_and_Process_NCEP_FNL_Analysis(date,dims,connection_info)

#Download and process the 3b42rt precipitation data
ml.Download_and_Process_3b42RT(date,dims)

#Finalize connection for FNL analysis
ml.Finalize_NCEP_FNL_Connection()

#Download gfs forecast
ml.Download_and_Process_GFS_forecast(date,dims)

#Restart the gds server
os.system('../LIBRARIES/gds-2.0/rebootserver')
'''
#grads_exe = '../LIBRARIES/grads-2.0.1.oga.1/Contents/grads'
#ga = grads.GaNum(Bin=grads_exe,Window=False,Echo=False)
idate = datetime.datetime(2001,1,1)
fdate = datetime.datetime(2008,12,31)
date = idate

'''
#Fit a gamma distribution to the 3B42RT product
ctl_in = "../DATA/3B42RT/MONTHLY/3B42RT_monthly_0.25deg.ctl"
file_out = "../WORKSPACE/parameters.nc"
dt = datetime.timedelta(days = 31)
ml.Fit_Distribution('gamma',ctl_in,idate,fdate,dt,'prec',file_out,dims)
'''
date = datetime.datetime(2013,1,1)
#Download and process the seasonal forecast
ml.Download_and_Process_Seasonal_Forecast(date)


'''
print date
#Download and process modis NDVI
ml.Download_and_Process_NDVI(date,dims)
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
