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
def datetime2gradstime(date):

 #Convert datetime to grads time
 str = date.strftime('%HZ%d%b%Y')
 
 return str

def Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims):

 #Open file in grads
 ga("xdfopen %s" % ctl_in)

 #What variables exist in the dataset?
 qh = ga.query("file")

 #Create and open access to netcdf file
 nt = 1
 tstep = 'hours'
 fp = ml.Create_NETCDF_File(dims,file_out,qh.vars,qh.var_info,idate,tstep,nt) #Create and open access to netcdf file

 #Iterate through the variables
 for var in qh.vars:
 
   #Determine what type of operation
   type = "ave"
   if var == "prec":
    type = "sum"

   #Convert times
   t1 = datetime2gradstime(idate)
   t2 = datetime2gradstime(fdate)

   #Perform operation
   ga("data = %s(%s,time=%s,time=%s)" % (type,var,t1,t2))

   #Write to file
   data = np.ma.getdata(ga.exp("data"))
   fp.variables[var][0] = data

 #Close access to the grads file
 ga("close 1")



 
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
grads_exe = '../LIBRARIES/grads-2.0.1.oga.1/Contents/grads'
ga = grads.GaNum(Bin=grads_exe,Window=False,Echo=False)
idate = datetime.datetime(2013,6,1)
fdate = datetime.datetime(2013,6,30)
date = idate
 
while date <= fdate:

 print date
 #Download and process modis NDVI
 ml.Download_and_Process_NDVI(date,dims)
 date = date + dt
'''
#End of month routines
ndate = date + dt
if date.month != ndate.month:
 print "Tomorrow is a new month. We have some extra work to do today."
 idate = datetime.datetime(date.year,date.month,1)
 fdate = date
 #3B42RT
 ctl_in = "../DATA/3B42RT/3B42RT_daily_0.25deg.ctl"
 file_out = "../DATA/3B42RT/MONTHLY/3B42RT_%04d%02d_daily_0.250deg.nc" % (idate.year,idate.month)
 Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims)
 #FNL
 ctl_in = "../DATA/FNL_ANALYSIS/fnlanl_daily_0.25deg.ctl"
 file_out = "../DATA/FNL_ANALYSIS/MONTHLY/fnlanl_%04d%02d_daily_0.250deg.nc" % (idate.year,idate.month)
 Compute_and_Output_Averages(ctl_in,file_out,idate,fdate,dims)
 
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
