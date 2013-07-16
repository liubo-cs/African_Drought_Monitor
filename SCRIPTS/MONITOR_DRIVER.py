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
import fileinput

#import matplotlib.pyplot as plt

def Grads_Regrid(var_in,var_out,dims):

 ga("%s = re(%s,%d,linear,%f,%f,%d,linear,%f,%f)" % (var_out,var_in,dims['nlon'],dims['minlon'],dims['res'],dims['nlat'],dims['minlat'],dims['res']))

 return

def Initialize_NCEP_FNL_Connection(user,pswd):

 #Login to server
 syscmd = "wget --no-check-certificate -O /dev/null --save-cookies auth.rda_ucar_edu --post-data=\"email=%s&passwd=%s&action=login\" =\"email=%s&passwd=%s&action=login\" https://rda.ucar.edu/cgi-bin/login" % (user,pswd,user,pswd)
 os.system(syscmd)

 return

def Finalize_NCEP_FNL_Connection():

 #Disconnect from server
 os.system('rm -f auth.rda_ucar_edu')

 return

def Download_NCEP_FNL_Analysis(date,root):

 #Download data
 grb_file = "fnl_%04d%02d%02d_%02d_00_c" % (date.year,date.month,date.day,date.hour)
 dwncmd = 'wget -N --no-check-certificate -P %s --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds083.2/grib1/%04d/%04d.%02d/' % (root,date.year,date.year,date.month)
 if os.path.isfile('%s/%s' % (root,grb_file)) == False:
  os.system("%s%s" % (dwncmd,grb_file))

def Download_and_Process_NCEP_FNL_Analysis(date,dims):

 workspace = '../WORKSPACE'
 fnl_analysis_root = '../DATA/FNL_ANALYSIS'
 dt = datetime.timedelta(hours=6)
 idate = date

 #Download date of NCEP Final Analysis (http://rda.ucar.edu/datasets/ds083.2/)
 for i in xrange(0,4):
  date = date + dt
  Download_NCEP_FNL_Analysis(date,workspace)
  
 #Create index and control file for the entire period
 ctl_file = 'fnl_%04d%02d%02d_00_c.ctl' % (idate.year,idate.month,idate.day)
 grb_file = 'fnl_%04d%02d%s_%s_00_c'% (idate.year,idate.month,'%d2','%h2')
 os.system('perl ../LIBRARIES/grib2ctl.pl %s/%s > %s/%s' % (workspace,grb_file,workspace,ctl_file))

 #Correct errors in ctl file
 os.system("sed -i 's/18hr/6hr/g' %s/%s" % (workspace,ctl_file))
 os.system("sed -i 's/tdef 2/tdef 4/g' %s/%s" % (workspace,ctl_file))

 #Create index file
 os.system('gribmap -0 -i %s/%s' % (workspace,ctl_file))
 
 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))

 #Define the output variables
 ga("tmax = max(TMAX2m,t=1,t=4)")
 ga("tmin = min(TMIN2m,t=1,t=4)")
 ga("prec = 3600*24*sum(pratesfc,t=1,t=4)/4")
 ga("wind = pow(pow(UGRD10m,2) + pow(VGRD10m,2),0.5)")

 #Regrid to 1/4 degree
 Grads_Regrid("tmax","tmax",dims)
 Grads_Regrid("tmin","tmin",dims)
 Grads_Regrid("prec","prec",dims)
 Grads_Regrid("wind","wind",dims)

 #Save to netcdf file
 netcdf_file = 'fnlanl_%04d%02d%02d_006_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 ga("set sdfwrite -flt -zip %s/%s" % (fnl_analysis_root,netcdf_file))
 ga("sdfwrite prec")
 ga("sdfwrite tmax")
 ga("sdfwrite tmin")
 ga("sdfwrite wind")

 #Close access to all files in grads
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

#Initialize connection for FNL analysis
pswd = 'ZlWBqFNK'
user = 'chaneyna@gmail.com'
Initialize_NCEP_FNL_Connection(user,pswd)

#Download and process the gfs analysis data
Download_and_Process_NCEP_FNL_Analysis(date,dims)

#Finalize connection for FNL analysis
Finalize_NCEP_FNL_Connection()
