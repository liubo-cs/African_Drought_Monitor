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
  f.createVariable(var,'f',('t','lat','lon'),fill_value=-99999)
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
 date_change = datetime.datetime(2008,9,30,12)
 new_file = "fnl_%04d%02d%02d_%02d_00" % (date.year,date.month,date.day,date.hour)
 grb_file = new_file
 if date >= date_change:
  grb_file = "fnl_%04d%02d%02d_%02d_00_c" % (date.year,date.month,date.day,date.hour)
 dwncmd = 'wget -N --no-check-certificate -O %s/%s --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds083.2/grib1/%04d/%04d.%02d/' % (root,new_file,date.year,date.year,date.month)
 os.system("%s%s" % (dwncmd,grb_file))

 return 

def Download_and_Process_NCEP_FNL_Analysis(date,dims,connection_info):

 workspace = '../WORKSPACE'
 fnl_analysis_root = '../DATA/FNL_ANALYSIS'
 dt = datetime.timedelta(hours=6)
 idate = date

 #Download date of NCEP Final Analysis (http://rda.ucar.edu/datasets/ds083.2/)
 for i in xrange(0,4):
  date = date + dt
  Download_NCEP_FNL_Analysis(date,workspace)
  
 #Create index and control file for the entire period
 ctl_file = 'fnl_%04d%02d%02d_00.ctl' % (idate.year,idate.month,idate.day)
 grb_file = 'fnl_%04d%02d%s_%s_00'% (idate.year,idate.month,'%d2','%h2')
 os.system('perl ../LIBRARIES/grib2ctl.pl %s/%s > %s/%s' % (workspace,grb_file,workspace,ctl_file))

 #Correct errors in ctl file
 os.system("sed -i 's/18hr/6hr/g' %s/%s" % (workspace,ctl_file))
 os.system("sed -i 's/tdef 2/tdef 4/g' %s/%s" % (workspace,ctl_file))

 #Create index file
 gribmap = '../LIBRARIES/grads-2.0.1.oga.1/Contents/gribmap'
 os.system('%s -0 -i %s/%s' % (gribmap,workspace,ctl_file))
 
 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))

 #Set region
 ga("set lat -89.5 89.5")
 ga("set lon -179.5 179.5")

 #Define the output variables
 ga("tmax = max(TMAX2m,t=1,t=4)")
 ga("tmin = min(TMIN2m,t=1,t=4)")
 ga("prec = 3600*24*sum(pratesfc,t=1,t=4)/4")
 ga("wind = pow(pow(ave(UGRD10m,t=1,t=4),2) + pow(ave(VGRD10m,t=1,t=4),2),0.5)")

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("tmax","tmax",dims)
 Grads_Regrid("tmin","tmin",dims)
 Grads_Regrid("prec","prec",dims)
 Grads_Regrid("wind","wind",dims)

 #Create and open access to netcdf file
 netcdf_file = 'fnlanl_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
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
 #Finalize_NCEP_FNL_Connection()

 #Finalize and close NETCDF file
 fp.close()

 #Remove files from the workspace
 os.system('rm -f %s/fnl_*' % workspace)

def Download_and_Process_3b42RT(date,dims):

 workspace = '../WORKSPACE'
 tmpa_3b42rt_root = '../DATA/3B42RT'
 dt = datetime.timedelta(hours=3)
 idate = date

 #Download date of TMPA 3b42rt product
 for i in xrange(0,9):
  ftp_root = 'ftp://disc2.nascom.nasa.gov/data/TRMM/Gridded/3B42RT/%04d%02d' % (date.year,date.month)
  #ftp_root = 'ftp://trmmopen.gsfc.nasa.gov/pub/merged/mergeIRMicro/%04d' % date.year
  file = "3B42RT.%04d.%02d.%02d.%02dz.bin" % (date.year,date.month,date.day,date.hour)
  ftp_file = '%s/%s' % (ftp_root,file)
  dwncmd = 'wget -P %s %s' % (workspace,ftp_file)
  #if os.path.isfile('%s/%s' % (workspace,file)) == False:
  os.system(dwncmd)
  date = date + dt

 #Construct the control file
 ctl_file = '3B42RT.%04d.%02d.%02d.ctl' % (idate.year,idate.month,idate.day)
 f = open('%s/%s' % (workspace,ctl_file), 'w')
 f.write('dset ^3B42RT.%y4.%m2.%d2.%h2z.bin\n')
 f.write('options template byteswapped\n')
 f.write('title Real-Time Three Hourly TRMM and Other Satellite Rainfall (3B42RT)\n')
 f.write('undef -99999.0\n')
 f.write('xdef 1440 linear 0.1250 0.25\n')
 f.write('ydef 480  linear -59.8750 0.25\n')
 f.write('zdef 1 levels 10000\n')
 f.write('tdef 9 linear 00Z%02d%s%04d 3hr\n' % (idate.day,idate.strftime('%b'),idate.year))
 f.write('vars 1\n')
 f.write('p 0 99 Precipitation (mm/hr)\n')
 f.write('endvars\n')
 f.close()

 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))

 #Set region
 ga("set lat -89.875 89.875")
 ga("set lon -179.875 179.875")

 #Define the output variables
 ga("tmp = 1.5*p(t=1) + 1.5*p(t=9) + 3*sum(p,t=2,t=8)") #mm/day
 ga("prec = maskout(tmp,tmp)")

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("prec","prec",dims)

 #Create and open access to netcdf file
 netcdf_file = '3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
 file = '%s/%s' % (tmpa_3b42rt_root,netcdf_file)
 vars = ['prec']
 vars_info = ['daily total precip (mm)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var))
  fp.variables[var][0] = data

 #Close access to all files in grads
 ga("close 1")

 #close output file
 fp.close()

 #Remove files from the workspace
 os.system('rm -f %s/3B42RT*' % workspace)

 return

def Download_and_Process_GFS_forecast(date,dims):

 print "Downloading and processing the gfs 7-day forecast"
 gfs_hd_root = '../DATA/GFS'

 #Establish connection to 00 forecast
 gds_file = 'http://nomads.ncdc.noaa.gov:80/dods/NCEP_GFS/%04d%02d/%04d%02d%02d/gfs_3_%04d%02d%02d_0000_fff' % (date.year,date.month,date.year,date.month,date.day,date.year,date.month,date.day)
 #gds_file = 'http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/gfs_hd%04d%02d%02d/gfs_hd_00z' % (date.year,date.month,date.day)
 ga("sdfopen %s" % gds_file)

 #Create and open access to netcdf file
 nt = 7
 netcdf_file = 'gfs_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
 file = '%s/%s' % (gfs_hd_root,netcdf_file)
 vars = ['tmax','tmin','prec','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily total precip (mm)','daily mean wind speed (m/s)']
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,tstep,nt)

 #Process each day
 t1 = 2
 t2 = 9
 for t in xrange(0,nt):
 
  #Set region
  ga("set lat -89.5 89.5")
  ga("set lon -179.5 179.5")

  #Define the output variables
  ga("tmax = max(TMAX2m,t=%d,t=%d)" % (t1,t2))
  ga("tmin = min(TMIN2m,t=%d,t=%d)" % (t1,t2))
  ga("prec = 3600*24*sum(oprate,t=%d,t=%d)/8" % (t1,t2))
  ga("wind = pow(pow(ave(UGRD10m,t=%d,t=%d),2) + pow(ave(VGRD10m,t=%d,t=%d),2),0.5)" % (t1,t2,t1,t2))

  #Set to new region
  ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
  ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

  #Regrid to 1/4 degree
  Grads_Regrid("tmax","tmax",dims)
  Grads_Regrid("tmin","tmin",dims)
  Grads_Regrid("prec","prec",dims)
  Grads_Regrid("wind","wind",dims)

  #Write to file
  for var in vars:
   data = np.ma.getdata(ga.exp(var))
   fp.variables[var][t] = data

  #Update time
  t1 = t1 + 8
  t2 = t2 + 8

 #Close access to all files in grads
 ga("close 1")

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
idate = datetime.datetime(2013,7,21)
fdate = datetime.datetime(2009,12,31)

#Load connection info (usernames and passwords)
data = np.loadtxt('CONNECTION_INFO.txt',skiprows=1,dtype={'names': ('domain','username','password'), 'formats': ('S100','S100','S100')})
connection_info = {}
connection_info['NCEP_FNL'] = {}
connection_info['NCEP_FNL']['password'] = data['password']
connection_info['NCEP_FNL']['username'] = data['username']

#Initialize connection for FNL analysis
pswd = connection_info['NCEP_FNL']['password']
user = connection_info['NCEP_FNL']['username']
Initialize_NCEP_FNL_Connection(user,pswd)

#Open connection to grads through pygrads
grads_exe = '../LIBRARIES/grads-2.0.1.oga.1/Contents/opengrads'
ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)

date = idate
#dt = datetime.timedelta(days=1)
#while date <= fdate:
# print date

 #Download and process the gfs analysis data
Download_and_Process_NCEP_FNL_Analysis(date,dims,connection_info)

 #Download and process the 3b42rt precipitation data
Download_and_Process_3b42RT(date,dims)
# date = date + dt

#Finalize connection for FNL analysis
Finalize_NCEP_FNL_Connection()

#Download gfs forecast
Download_and_Process_GFS_forecast(date,dims)

#Reboot the gds server

