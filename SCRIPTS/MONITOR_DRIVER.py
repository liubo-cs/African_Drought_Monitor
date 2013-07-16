#################################################################
# Drought Monitor V2
# Written by: Nathaniel W. Chaney
# Date: 12 July 2013
# Location: Princeton University
# Contact: nchaney@princeton.edu
# Purpose: Main driver of the drought monitor
##################################################################

import datetime
import os
import grads
import numpy as np

#import matplotlib.pyplot as plt

def Grads_Regrid(var_in,var_out,dims):

 ga("%s = re(%s,%d,linear,%f,%f,%d,linear,%f,%f)" % (var_out,var_in,dims['nlon'],dims['minlon'],dims['res'],dims['nlat'],dims['minlat'],dims['res']))

 return

def Download_NCEP_FNL_Analysis(date):

 pswd = 'ZlWBqFNK'
 user = 'chaneyna@gmail.com'

 #Login to server
 syscmd = "wget --no-check-certificate -O /dev/null --save-cookies auth.rda_ucar_edu --post-data=\"email=%s&passwd=%s&action=login\" =\"email=%s&passwd=%s&action=login\" https://rda.ucar.edu/cgi-bin/login" % (user,pswd,user,pswd)
 os.system(syscmd)

 #Download data
 dwncmd = 'wget -N --no-check-certificate --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds083.2/'
 file = "grib1/2013/2013.07/fnl_20130715_00_00_c"
 os.system("%s%s" % (dwncmd,file))

 #Disconnect from server
 os.system('rm -f auth.rda_ucar_edu')

def Download_and_Process_GFS_Analysis(date,dims):

 #Analysis (Complete)
 workspace = '../WORKSPACE'
 gfs_analysis_root = '../DATA/GFS_ANALYSIS'
 hours = [0,6,12,18]

 #Download all the analysis files for the day (+6 hours)
 dt = datetime.timedelta(hours=6)
 idate = date + dt
 fdate = idate + 3*dt
 tmp = idate
 while tmp < fdate:
  ftp_root = 'ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/%04d%02d/%04d%02d%02d' % (tmp.year,tmp.month,tmp.year,tmp.month,tmp.day)
  grb2_file = 'gfsanl_4_%04d%02d%02d_%02d00_000.grb2'  % (tmp.year,tmp.month,tmp.day,tmp.hour)
  if os.path.isfile('%s/%s' % (workspace,grb2_file)) == False:
   os.system('wget -P %s %s/%s' % (workspace,ftp_root,grb2_file))
  tmp = tmp + dt

 #Create index and control file for the entire period
 ctl_file = 'gfsanl_4_%04d%02d%02d_000.ctl' % (date.year,date.month,date.day)
 grb2_file = 'gfsanl_4_%y4%m2%d2_%h200_000.grb2'
 print 'perl ../LIBRARIES/g2ctl -0 %s/%s > %s/%s' % (workspace,grb2_file,workspace,ctl_file)
 os.system('perl ../LIBRARIES/g2ctl -0 %s/%s > %s/%s' % (workspace,grb2_file,workspace,ctl_file))
 os.system('gribmap -0 -i %s/%s' % (workspace,ctl_file))

 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))
 '''
 #3hr Forecast (3 hours after each analysis period)
 workspace = '../WORKSPACE'
 gfs_analysis_root = '../DATA/GFS_ANALYSIS'
 hours = [0,6,12,18]
 ftp_root = 'ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/%04d%02d/%04d%02d%02d' % (date.year,date.month,date.year,date.month,date.day)

 #Download all the analysis files for the day (+6 hours)
 for hour in hours:
  grb2_file = 'gfsanl_4_%04d%02d%02d_%02d00_003.grb2'  % (date.year,date.month,date.day,hour)
  if os.path.isfile('%s/%s' % (workspace,grb2_file)) == False:
   os.system('wget -P %s %s/%s' % (workspace,ftp_root,grb2_file))

 #Create index and control file for the entire period
 ctl_file = 'gfsanl_4_%04d%02d%02d_003.ctl' % (date.year,date.month,date.day)
 grb2_file = 'gfsanl_4_%04d%02d%02d_%sh200_003.grb2' % (date.year,date.month,date.day,'%')
 os.system('perl ../LIBRARIES/g2ctl -0 %s/%s > %s/%s' % (workspace,grb2_file,workspace,ctl_file))
 os.system('gribmap -0 -i %s/%s' % (workspace,ctl_file))

 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))

 #6hr Forecast (6 hours after each analysis period)
 workspace = '../WORKSPACE'
 gfs_analysis_root = '../DATA/GFS_ANALYSIS'
 hours = [0,6,12,18]
 ftp_root = 'ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/%04d%02d/%04d%02d%02d' % (date.year,date.month,date.year,date.month,date.day)
 
 #Download all the analysis files for the day (+6 hours)
 for hour in hours:
  grb2_file = 'gfsanl_4_%04d%02d%02d_%02d00_006.grb2'  % (date.year,date.month,date.day,hour)
  if os.path.isfile('%s/%s' % (workspace,grb2_file)) == False:
   os.system('wget -P %s %s/%s' % (workspace,ftp_root,grb2_file))
 
 #Create index and control file for the entire period
 ctl_file = 'gfsanl_4_%04d%02d%02d_006.ctl' % (date.year,date.month,date.day)
 grb2_file = 'gfsanl_4_%04d%02d%02d_%sh200_006.grb2' % (date.year,date.month,date.day,'%')
 os.system('perl ../LIBRARIES/g2ctl -0 %s/%s > %s/%s' % (workspace,grb2_file,workspace,ctl_file))
 os.system('gribmap -0 -i %s/%s' % (workspace,ctl_file))

 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))
 '''
 exit()

 #Define the output variables
 ga("set t 1")
 tmp2m = []
 for t in xrange(1,5):
  tmp2m.append(ga.exp("tmp2m"))
 #tmax = np.max(tmp2m,axis=0)
 #tmin = np.min(tmp2m,axis=0)
 #print tmax
 ga.imp("tmax",np.max(tmp2m,axis=0))
 #ga.imp("tmin",tmin)
 
 
 ga("tmax = 0")
 ga("tmin = 0")
 ga("prec = sum(apcpsfc.2,t=1,t=4) + sum(apcpsfc.2,t=1,t=4)")
 ga("wind = 0")

 #Calculate the output variables for each time step
 #for ifile in xrange(1,len(hours)+1):
 # ga("tmax = tmax 
 #Regrid to 1/4 degree
 var_out = "prec"
 var_in = "prec"
 Grads_Regrid(var_in,var_out,dims)

 #Save to netcdf file
 netcdf_file = 'gfsanl_4_%04d%02d%02d_006_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 ga("set sdfwrite -flt -zip %s/%s" % (gfs_analysis_root,netcdf_file))
 ga("sdfwrite prec")

 #Make daily file with desired variables
 
 #for ifile in xrange(1,len(hours)+1):
  
  

 #Close access to all files in grads
 ga("close 3")
 ga("close 2")
 ga("close 1")
  

#def Download_and_Process_3b42RT():

 #Download all precipitation for the day

 #Create control file

 #Make daily file   
 

 #1. Download global GFS analysis fields


#1. Determine the period that needs to be updated

#2. Download all the relevant data

#3. Process all the relevant data

#4. Run all the relevant models

date = datetime.datetime(2013,6,12)
dims = {}
dims['minlat'] = -89.8750
dims['minlon'] = 0.1250
dims['nlat'] = 720
dims['nlon'] = 1440
dims['res'] = 0.250

#Open connection to grads through pygrads
#ga = grads.GaNum(Bin='grads',Window=False,Echo=False)
ga = grads.GrADS(Bin='grads',Window=False,Echo=False)

#Download and process the gfs analysis data

Download_and_Process_GFS_Analysis(date,dims)

