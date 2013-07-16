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
import netCDF4 as netcdf

def Create_NETCDF_File(dims,file,vars,vars_info,tinitial,tstep,nt):

 nlat = dims['nlat']
 nlon = dims['nlon']
 res = dims['res']
 minlon = dims['minlon']
 minlat = dims['minlat']
 t = np.arange(0,nt)

 #Prepare the netcdf file
 #Create file
 f = netcdf.Dataset(file, 'w')

 #Define dimensions
 f.createDimension('lon',nlon)
 f.createDimension('lat',nlat)
 f.createDimension('t',len(t))

 #Longitude
 f.createVariable('lon','d',('lon',))
 f.variables['lon'][:] = np.linspace(minlon,minlon+res*(nlon-1),nlon)
 f.variables['lon'].units = 'degrees_east'
 f.variables['lon'].long_name = 'Longitude'
 f.variables['lon'].res = res

 #Latitude
 f.createVariable('lat','d',('lat',))
 f.variables['lat'][:] = np.linspace(minlat,minlat+res*(nlat-1),nlat)
 f.variables['lat'].units = 'degrees_north'
 f.variables['lat'].long_name = 'Latitude'
 f.variables['lat'].res = res

 #Time
 times = f.createVariable('t','d',('t',))
 f.variables['t'][:] = t
 f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (tstep,tinitial.year,tinitial.month,tinitial.day,tinitial.hour)
 f.variables['t'].long_name = 'Time'

 #Data
 i = 0
 for var in vars:
  f.createVariable(var,'f',('t','lat','lon'),fill_value=-9.99e+08)
  f.variables[var].long_name = vars_info[i]
  i = i + 1

 return f

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

def Download_and_Process_NCEP_FNL_Analysis(date,dims,connection_info):

 workspace = '../WORKSPACE'
 fnl_analysis_root = '../DATA/FNL_ANALYSIS'
 dt = datetime.timedelta(hours=6)
 idate = date

 #Initialize connection for FNL analysis
 pswd = connection_info['NCEP_FNL']['password']
 user = connection_info['NCEP_FNL']['username']
 Initialize_NCEP_FNL_Connection(user,pswd)

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

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("tmax","tmax",dims)
 Grads_Regrid("tmin","tmin",dims)
 Grads_Regrid("prec","prec",dims)
 Grads_Regrid("wind","wind",dims)

 #Create and open access to netcdf file
 netcdf_file = 'fnlanl_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 file = '%s/%s' % (fnl_analysis_root,netcdf_file)
 vars = ['tmax','tmin','prec','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily total precip (mm)','daily mean wind speed (m/s)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var))
  fp.variables[var][0] = data

 #Close access to all files in grads
 ga("close 1")

 #Finalize connection for FNL analysis
 Finalize_NCEP_FNL_Connection()

 #Finalize and close NETCDF file
 fp.close()

 #Remove files from the workspace
 os.system('rm -f %s/fnl_*' % workspace)

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
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)

#Load connection info (usernames and passwords)
data = np.loadtxt('CONNECTION_INFO.txt',skiprows=1,dtype={'names': ('domain','username','password'), 'formats': ('S100','S100','S100')})
connection_info = {}
connection_info['NCEP_FNL'] = {}
connection_info['NCEP_FNL']['password'] = data['password']
connection_info['NCEP_FNL']['username'] = data['username']

#Open connection to grads through pygrads
ga = grads.GrADS(Bin='grads',Window=False,Echo=False)

#Download and process the gfs analysis data
Download_and_Process_NCEP_FNL_Analysis(date,dims,connection_info)
