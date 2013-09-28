#################################################################
# Drought Monitor V2
# Written by: Nathaniel W. Chaney
# Date: 12 July 2013
# Location: Princeton University
# Contact: nchaney@princeton.edu
# Purpose: Library of routines for the "master" of the monitor
##################################################################

import datetime
import os
import grads
import numpy as np
import fileinput
import netCDF4 as netcdf
import pyhdf.SD as sd
import library_f90
import scipy.stats as ss
import subprocess
import dateutil.relativedelta as relativedelta
import time
import random
import cPickle as pickle
grads_exe = '../LIBRARIES/grads-2.0.1.oga.1/Contents/grads'
ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)

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

def Initialize_NCEP_FNL_Connection():

 print "Initializing FNL connection"
 data = np.loadtxt('CONNECTION_INFO.txt',skiprows=1,dtype={'names': ('domain','username','password'), 'formats': ('S100','S100','S100')})
 connection_info = {}
 connection_info['NCEP_FNL'] = {}
 connection_info['NCEP_FNL']['password'] = data['password']
 connection_info['NCEP_FNL']['username'] = data['username']

 #Initialize connection for FNL analysis
 pswd = connection_info['NCEP_FNL']['password']
 user = connection_info['NCEP_FNL']['username']

 #Login to server
 syscmd = "wget -nv --no-check-certificate -O /dev/null --save-cookies auth.rda_ucar_edu --post-data=\"email=%s&passwd=%s&action=login\" https://rda.ucar.edu/cgi-bin/login" % (user,pswd)
 subprocess.call(syscmd,shell=True)

 return

def Finalize_NCEP_FNL_Connection():

 print "Finalizing FNL connection"
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
 dwncmd = 'wget -nv --no-check-certificate -O %s/%s --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds083.2/grib1/%04d/%04d.%02d/' % (root,new_file,date.year,date.year,date.month)
 os.system("%s%s" % (dwncmd,grb_file))

 return 

def Download_and_Process_NCEP_FNL_Analysis(date,dims,idate,fdate,Reprocess_Flag):
 
 workspace = '../WORKSPACE'
 fnl_analysis_root = '../DATA/FNL_ANALYSIS/DAILY'
 dt = datetime.timedelta(hours=6)

 if date == idate:
  #Initialize connection for FNL analysis
  connection_info = Initialize_NCEP_FNL_Connection()

 #If the date is before the product's start date:
 if date < datetime.datetime(2010,1,1):
  return

 idate = date
 #If the file already exists exit:
 file = 'fnlanl_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
 netcdf_file = '%s/%s' % (fnl_analysis_root,file)
 if os.path.exists(netcdf_file) == True and Reprocess_Flag == False:
  return
 
 print_info_to_command_line("Downloading and processing the NCEP GFS final analysis")

 #Download date of NCEP Final Analysis (http://rda.ucar.edu/datasets/ds083.2/)
 for i in xrange(0,4):
  date = date + dt
  Download_NCEP_FNL_Analysis(date,workspace)
  
 #Create index and control file for the entire period
 ctl_file = 'fnl_%04d%02d%02d_00.ctl' % (idate.year,idate.month,idate.day)
 grb_file = 'fnl_%04d%02d%s_%s_00'% (idate.year,idate.month,'%d2','%h2')
 os.system('perl ../LIBRARIES/grib2ctl.pl %s/%s 1> %s/%s 2> /dev/null' % (workspace,grb_file,workspace,ctl_file))

 #Correct errors in ctl file
 os.system("sed -i 's/18hr/6hr/g' %s/%s" % (workspace,ctl_file))
 os.system("sed -i 's/tdef 2 linear 18Z/tdef 4 linear 06Z/g' %s/%s" % (workspace,ctl_file))
 os.system("sed -i 's/_18_/_6_/g' %s/%s" % (workspace,ctl_file))

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
 ga("prec = 3600*24*ave(pratesfc,t=1,t=4)")
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
 #netcdf_file = 'fnlanl_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
 #file = '%s/%s' % (fnl_analysis_root,netcdf_file)
 vars = ['tmax','tmin','prec','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily total precip (mm)','daily mean wind speed (m/s)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,netcdf_file,vars,vars_info,idate,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var))
  fp.variables[var][0] = data

 #Close access to all files in grads
 ga("close 1")

 #Finalize and close NETCDF file
 fp.close()

 #Remove files from the workspace
 os.system('rm -f %s/fnl_*' % workspace)

def Download_and_Process_GFS_Analysis(date,dims,Reprocess_Flag):

 gfs_anl_root = '../DATA/GFS_ANL/DAILY'
 workspace = '../WORKSPACE'
 #If the data is before the product's start date return:
 itime = datetime.datetime(2008,1,1)
 if date < itime:
  return

 #If the file already exists exit:
 file = 'GFS_ANL_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 netcdf_file = '%s/%s' % (gfs_anl_root,file)
 if os.path.exists(netcdf_file) == True and Reprocess_Flag == False:
  return
 
 print_info_to_command_line("Downloading and processing the NCEP GFS analysis")

 #Download GFS Analysis for the entire day
 idate = date
 dt = datetime.timedelta(hours=6)
 for i in xrange(0,4):
  date2 = date
  while os.path.exists('../WORKSPACE/gfsanl_3_%04d%02d%02d_%02d00_000.grb' % (date.year,date.month,date.day,date.hour)) == False:
   ftp_file = 'ftp://nomads.ncdc.noaa.gov/GFS/analysis_only/%04d%02d/%04d%02d%02d/gfsanl_3_%04d%02d%02d_%02d00_000.grb' % (date2.year,date2.month,date2.year,date2.month,date2.day,date2.year,date2.month,date2.day,date2.hour)
   cmd = 'wget -nv -P ../WORKSPACE %s' % ftp_file 
   os.system(cmd)
   if (date - date2).days > 0:
    os.system("mv ../WORKSPACE/gfsanl_3_%04d%02d%02d_%02d00_000.grb ../WORKSPACE/gfsanl_3_%04d%02d%02d_%02d00_000.grb" % (date2.year,date2.month,date2.day,date2.hour,date.year,date.month,date.day,date.hour))
   date2 = date2 - datetime.timedelta(days=1)
   if (date - date2).days > 10: 
    print "Missing GFS analysis data for more than 10 days"
    return
  date = date + dt

 #Create index and control file for the entire period
 ctl_file = 'gfsanl_3_%04d%02d%02d_0000_000.ctl' % (idate.year,idate.month,idate.day)
 grb_file = 'gfsanl_3_%04d%02d%02d_%s00_000.grb'% (idate.year,idate.month,idate.day,'%h2')
 os.system('perl ../LIBRARIES/grib2ctl.pl %s/%s 1> %s/%s 2> /dev/null' % (workspace,grb_file,workspace,ctl_file))

 #Create index file
 gribmap = '../LIBRARIES/grads-2.0.1.oga.1/Contents/gribmap'
 os.system('%s -0 -i %s/%s' % (gribmap,workspace,ctl_file))

 #Open access to the file
 ga("open %s/%s" % (workspace,ctl_file))

 #Set region
 ga("set lat -89.5 89.5")
 ga("set lon -179.5 179.5")

 #Define the output variables
 ga("tmax0 = max(tmp2m,t=1,t=4)")
 ga("tmin0 = min(tmp2m,t=1,t=4)")
 ga("wind0 = pow(pow(ave(UGRD10m,t=1,t=4),2) + pow(ave(VGRD10m,t=1,t=4),2),0.5)")

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("tmax0","tmax0",dims)
 Grads_Regrid("tmin0","tmin0",dims)
 Grads_Regrid("wind0","wind0",dims)

 #Create and open access to netcdf file
 vars = ['tmax','tmin','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily mean wind speed (m/s)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,netcdf_file,vars,vars_info,idate,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var+'0'))
  fp.variables[var][0] = data

 #Close access to all files in grads
 ga("close 1")

 #Remove files from the workspace
 os.system('rm -f %s/gfsanl_3*' % workspace)

 #Update the control file
 ctl_out = '../DATA/GFS_ANL/DAILY/gfsanl_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('GFS_ANL','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 return

def Download_and_Process_3b42RT(date,dims,Reprocess_Flag):

 workspace = '../WORKSPACE'
 tmpa_3b42rt_root = '../DATA/3B42RT/DAILY'
 dt = datetime.timedelta(hours=3)
 idate = date

 #If the date is before the product's start date:
 itime = datetime.datetime(2000,3,1)
 if date < itime:
  return

 #If the file already exists exit:
 file = '3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,dims['res'])
 netcdf_file = '%s/%s' % (tmpa_3b42rt_root,file)
 if os.path.exists(netcdf_file) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Downloading and processing the 3b42RT product")

 #Download date of TMPA 3b42rt product
 for i in xrange(0,9):
  ftp_root = 'ftp://disc2.nascom.nasa.gov/data/TRMM/Gridded/3B42RT/%04d%02d' % (date.year,date.month)
  #ftp_root = 'ftp://trmmopen.gsfc.nasa.gov/pub/merged/mergeIRMicro/%04d' % date.year
  file = "3B42RT.%04d.%02d.%02d.%02dz.bin" % (date.year,date.month,date.day,date.hour)
  ftp_file = '%s/%s' % (ftp_root,file)
  dwncmd = 'wget -nv -P %s %s' % (workspace,ftp_file)
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
 ga("prec0 = const(maskout(tmp,tmp),-9.99e+08,-u)")

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("prec0","prec0",dims)

 #Create and open access to netcdf file
 vars = ['prec']
 vars_info = ['daily total precip (mm)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,netcdf_file,vars,vars_info,idate,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var+'0'))
  fp.variables[var][0] = data

 #Close access to all files in grads
 ga("close 1")

 #close output file
 fp.close()

 #Remove files from the workspace
 os.system('rm -f %s/3B42RT*' % workspace)

 #Update the control file
 ctl_out = '../DATA/3B42RT/DAILY/3B42RT_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('3B42RT','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 return

def Download_and_Process_GFS_forecast(date,dims,Reprocess_Flag):

 gfs_hd_root = '../DATA/GFS'
 date = date + datetime.timedelta(days=1) #7-day forecast starting tomorrow
 idate = date
 dir = '../DATA/GFS/%04d%02d%02d' % (date.year,date.month,date.day)
 dir_derived = '../DATA/GFS_DERIVED/%04d%02d%02d' % (idate.year,idate.month,idate.day)
 if os.path.exists(dir_derived) == False:
  os.mkdir(dir_derived)

 #If the date is before the product's start date:
 date_tmp = datetime.datetime.today() - relativedelta.relativedelta(years=1)
 if date < datetime.datetime(date_tmp.year,date_tmp.month,date_tmp.day):
  return

 #If the directory already exists exit:
 #print os.stat.getsize(dir)
 #dir_size = os.popen('du %s' % dir).read().split()[0]
 if os.path.exists(dir) == False and Reprocess_Flag == False:
  os.system("mkdir %s" % dir)

 print_info_to_command_line("Downloading and processing the gfs 7-day forecast")

 #Establish connection to 00 forecast
 gds_file = 'http://nomads.ncdc.noaa.gov:80/dods/NCEP_GFS/%04d%02d/%04d%02d%02d/gfs_3_%04d%02d%02d_0000_fff' % (date.year,date.month,date.year,date.month,date.day,date.year,date.month,date.day)
 gds_file = 'http://nomads.ncep.noaa.gov:9090/dods/gfs/gfs%04d%02d%02d/gfs_00z' % (date.year,date.month,date.day)
 
 #gds_file = 'http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/gfs_hd%04d%02d%02d/gfs_hd_00z' % (date.year,date.month,date.day)
 try:
  ga("sdfopen %s" % gds_file)
 except:
  return

 #Define info for the netcdf files
 vars = ['tmax','tmin','prec','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily total precip (mm)','daily mean wind speed (m/s)']

 #Process each day
 t1 = 2
 t2 = 9
 nt = 7
 for t in xrange(0,nt):
  
  file = dir + '/gfs_%04d%02d%02d_daily_%.3fdeg_day%d.nc' % (idate.year,idate.month,idate.day,dims['res'],t+1)
  #Determine if we skip the time step
  if os.path.exists(file) == True and Reprocess_Flag == False:
   continue
 
  #Set region
  ga("set lat -89.5 89.5")
  ga("set lon -179.5 179.5")

  #Define the output variables
  ga("tmax0 = max(TMAX2m,t=%d,t=%d)" % (t1,t2))
  ga("tmin0 = min(TMIN2m,t=%d,t=%d)" % (t1,t2))
  #ga("prec = 3600*24*ave(oprate,t=%d,t=%d)" % (t1,t2))
  ga("prec0 = 3600*24*ave(pratesfc,t=%d,t=%d)" % (t1,t2))
  ga("wind0 = pow(pow(ave(UGRD10m,t=%d,t=%d),2) + pow(ave(VGRD10m,t=%d,t=%d),2),0.5)" % (t1,t2,t1,t2))

  #Set to new region
  ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
  ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

  #Regrid to 1/4 degree
  Grads_Regrid("tmax0","tmax0",dims)
  Grads_Regrid("tmin0","tmin0",dims)
  Grads_Regrid("prec0","prec0",dims)
  Grads_Regrid("wind0","wind0",dims)

  #Write to file
  fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,'days',1)
  for var in vars:
   data = np.ma.getdata(ga.exp(var+"0")) 
   fp.variables[var][0] = data
  fp.close()

  #Update control file
  date_ctl = datetime.datetime(2012,8,14) 
  ndays = (date - date_ctl).days + 1
  ctl_new = '../DATA/GFS/gfs_daily_0.250deg_day%d.ctl' % (t+1)
  fp = open(ctl_new,'w')
  fp.write('dset ^%s%s%s/gfs_%s%s%s_daily_0.250deg_day%d.nc\n' % ('%y4','%m2','%d2','%y4','%m2','%d2',t+1))
  fp.write('options template\n')
  fp.write('dtype netcdf\n')
  fp.write('tdef t %d linear 14aug2012 1dy\n' % ndays)
  fp.close()

  #Update time
  t1 = t1 + 8
  t2 = t2 + 8

 #Close access to all files in grads
 ga("close 1")

 return

def Download_and_Process_NDVI(date,dims,Reprocess_Flag):

 workspace = '../WORKSPACE'
 modis_root = '../DATA/MOD09_NDVI/DAILY' 

 #If the date is before the product's start date:
 itime = datetime.datetime(2003,1,1)
 if date < itime:
  return

 #If the file already exists exit:
 file = 'MOD09CMG_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 netcdf_file = '%s/%s' % (modis_root,file)
 if os.path.exists(netcdf_file) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Downloading and processing the MODIS MOD09CMG product (NDVI)")

 http_dir = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD09CMG.005/%04d.%02d.%02d/' % (date.year,date.month,date.day)
 os.system('wget -nv -P %s -r -l2 --no-parent -A.hdf %s' % (workspace,http_dir))
 os.system("mv ../WORKSPACE/e4ftl01.cr.usgs.gov/MOLT/MOD09CMG.005/%04d.%02d.%02d/* ../WORKSPACE/MOD09CMG.005_%04d.%02d.%02d.hdf" % (date.year,date.month,date.day,date.year,date.month,date.day))
 os.system("rm -rf ../WORKSPACE/e4ftl01.cr.usgs.gov")

 #HDF file and SDS names
 filename = '../WORKSPACE/MOD09CMG.005_%04d.%02d.%02d.hdf' % (date.year,date.month,date.day)

 #Check to see if file exists
 if os.path.exists(filename) == False:
  print "File does not exist"
  return
 '''
 if os.path.exists(filename) == False:
  print "File does not exist, placing an empty holder"
  old_file = '../DATA/MOD09_NDVI/DAILY/MOD09CMG_20030201_daily_0.250deg.nc'
  os.system("cp %s %s" % (old_file,netcdf_file))
  return
 '''

 #Open the hdf file
 hdf = sd.SD(filename)

 #Select and read the sds data
 sdsname = 'Coarse Resolution Surface Reflectance Band 1'
 sds = hdf.select(sdsname)
 band1= np.float32(sds.get())
 sdsname = 'Coarse Resolution Surface Reflectance Band 2'
 sds = hdf.select(sdsname)
 band2= np.float32(sds.get())
 #Calculate NDVI
 NDVI = (band2-band1)/(band2+band1)
 NDVI[NDVI>1] = float('NaN')
 NDVI[NDVI<=0] = float('NaN')
 sdsname = 'Coarse Resolution State QA'
 try:
  sds = hdf.select(sdsname)
  cqc = np.int32(sds.get())
  #Create the qc masks
  #Cloud shadow
  cloud_shadow = getBit_array(cqc,2)
  #land/water flag
  tmp1 = getBit_array(cqc,3)
  tmp2 = getBit_array(cqc,4)
  tmp3 = getBit_array(cqc,5)
  land_water_flag = np.ones(cqc.shape)
  land_water_flag[(tmp1 == 1) & (tmp2 == 0) & (tmp3 == 0)] = 0
  #cirrus flag
  tmp1 = getBit_array(cqc,8)
  tmp2 = getBit_array(cqc,9)
  cirrus_flag = np.ones(cqc.shape)
  cirrus_flag[(tmp1 == 0) & (tmp2 == 0)] = 0
  #internal cloud algorithm
  internal_cloud_flag = getBit_array(cqc,10)
  #pixel adjacent to cloud
  pixel_adjacent_cloud_flag = getBit_array(cqc,13)
  #cloud flag
  tmp1 = getBit_array(cqc,0)
  tmp2 = getBit_array(cqc,1)
  cqc[:] = float('NaN')
  cqc = np.ones(cqc.shape)
  cqc[(tmp1 == 0) & (tmp2 == 0)] = 0
  NDVI[cqc > 0]  = float('NaN')
  NDVI[cloud_shadow > 0] = float('NaN')
  NDVI[land_water_flag > 0] = float('NaN')
  NDVI[cirrus_flag > 0] = float('NaN')
  NDVI[internal_cloud_flag > 0] = float('NaN')
  NDVI[pixel_adjacent_cloud_flag > 0] = float('NaN')

 except:
  print "We do not have state QA data"
  NDVI[:] = float('NaN')

 #Terminate access to the data set
 sds.endaccess()

 #Close the file
 hdf.end()

 #Define MODIS dimensions
 mdims = {}
 mdims['nlat'] = NDVI.shape[0] #720
 mdims['nlon'] = NDVI.shape[1] #1440
 mdims['res'] = 360.0/mdims['nlon']#i0.250
 mdims['minlat'] = -90.0 + mdims['res']/2
 mdims['minlon'] = -180.0 + mdims['res']/2
 mdims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
 mdims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)

 #Create and open access to netcdf file
 file = '%s/modis_ndvi.nc' % workspace
 vars = ['ndvi']
 vars_info = ['NDVI (MOD09CMG.005)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(mdims,file,vars,vars_info,date,tstep,nt)

 #Write to file
 fp.variables['ndvi'][0] = np.flipud(NDVI)

 #Close file
 fp.close()

 #Open access to file
 ga("sdfopen %s" % file)

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Regrid to 1/4 degree
 Grads_Regrid("ndvi","data",dims)

 #Create and open access to netcdf file
 nt = 1
 vars = ['ndvi']
 vars_info = ['ndvi']
 tstep = 'days'
 fp = Create_NETCDF_File(dims,netcdf_file,vars,vars_info,date,tstep,nt)

 #Add data to file
 fp.variables['ndvi'][0] = np.ma.getdata(ga.exp("data"))

 #Close access to grads file
 ga("close 1")

 #Remove unnecessary files
 os.system("rm ../WORKSPACE/*.nc")
 os.system("rm ../WORKSPACE/*.hdf")
 
 #Close the output file
 fp.close()

 #Update the control file
 ctl_out = '../DATA/MOD09_NDVI/DAILY/MOD09CMG_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('MOD09CMG','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 return

def getBit_array(x, p):
 
 #Determine the value of the pth bit position (0,1) 
 #array cell_clim)
 array = np.zeros(x.shape,dtype=int)
 array[:] = 2**p
 array = np.bitwise_and(x,array)
 array[array!=0] = 1
 return array

def Compute_NDVI_moving_average(date,dims,Reprocess_Flag):

 dt_array = [datetime.timedelta(days=60),datetime.timedelta(days=30),datetime.timedelta(days=20),datetime.timedelta(days=10),datetime.timedelta(days=5),datetime.timedelta(days=1)]
 modis_root = '../DATA/MOD09_NDVI_MA/DAILY'
 netcdf_file = 'MOD09CMG_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 file = '%s/%s' % (modis_root,netcdf_file)
 file_org = "../DATA/MOD09_NDVI/DAILY/%s" % netcdf_file
 if os.path.exists(file_org) == False:
  return

 #If the date is before the product's start date:
 itime = datetime.datetime(2003,1,1)
 if date < itime:
  return

 #Update the control file
 ctl_out = '../DATA/MOD09_NDVI_MA/DAILY/MOD09CMG_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('MOD09CMG','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the file exists then exit (unless we are reprocessing the data)
 if os.path.exists(netcdf_file) == True and Reprocess_Flag == False:
  return

 #If file exists exit
 if os.path.exists(file) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Computing the MODIS MOD09CMG product n-day moving averages (NDVI)")

 #Open access to the NDVI
 ctl_file = "../DATA/MOD09_NDVI/DAILY/MOD09CMG_daily_0.25deg.ctl"

 #Open access to files
 ga("xdfopen %s" % ctl_file)

 #Create file
 nt = 1
 file = '%s/%s' % (modis_root,netcdf_file)
 vars = []
 vars_info = []
 for dt in dt_array:
  vars.append('ndvi%d' % dt.days)
  vars_info.append('ndvi %d-day moving window' % dt.days)
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file,vars,vars_info,date,tstep,nt)

 #Compute averages for all the moving averages
 for dt in dt_array:

  var = 'ndvi%d' % dt.days
  #Compute the average between the two intervals
  time1 = datetime2gradstime(date - dt/2)
  time2 = datetime2gradstime(date + dt/2)
  ga("data = ave(ndvi,time=%s,time=%s)" % (time1,time2))

  #Add data to file
  fp.variables[var][0] = np.ma.getdata(ga.exp("data"))

 #Close access to grads file
 ga("close 1")

 #Close access to new file
 fp.close()

def Extract_Data_Period_Average(idate_out,fdate_out,dt_down,dt_up,dt,ctl_in,var,type,open_type):

 #Open access to the control file
 ga("%s %s" % (open_type,ctl_in))
 #ga.open(ctl_in)

 #Determine initial and final time step
 ga('set t 1')
 idate_all = gradstime2datetime(ga.exp(var).grid.time[0])
 ga('set t last')
 fdate_all = gradstime2datetime(ga.exp(var).grid.time[0])

 date = idate_out
 count = 0
 while date <= fdate_out:

  #print date
  date1 = date - dt_down
  date2 = date + dt_up
  if date1 < idate_all and date2 < idate_all:
   date = date + dt
   continue
  if date1 < idate_all and date2 >= idate_all:
   date1 = idate_all
  if date1 > fdate_all and date2 > fdate_all:
   date = date +dt
   continue
  if date1 <= fdate_all and date2 > fdate_all:
   date2 = fdate_all
  t1 = datetime2gradstime(date1)
  t2 = datetime2gradstime(date2)

  #Extract data
  if type == "ave":
   ga("data = ave(%s,time=%s,time=%s)" % (var,t1,t2))
   tmp = np.ma.getdata(ga.exp("maskout(data,data)"))
 
  if type == "all":
   ga("set time %s %s" % (t1,t2))
   tmp = np.ma.getdata(ga.exp("maskout(%s,%s)" % (var,var)))

  #Convert 2-d arrays to 3-d for stacking
  if len(tmp.shape) == 2:
   tmp = np.reshape(tmp,(1,tmp.shape[0],tmp.shape[1]))

  #Append data
  if count == 0:#date == idate_out:
   count = 1
   data = tmp
  else:
   data = np.vstack((data,tmp))

  #Update time step
  date = date + dt 
 
 #Close access to the grads data
 ga("close 1")

 return data

def Calculate_Percentiles(data_clim,data):

 pct = np.zeros(data.shape)
 pct[0,:,:] = float('NaN')

 #Iterate through all cells and calculate percentiles
 for i in xrange(0,data.shape[1]):
  #print i
  for j in xrange(0,data.shape[2]):
   if data[0,i,j] == -9.99e+08:
    #pct[0,i,j] = float('NaN')
    continue
   data_cell_clim = data_clim[:,i,j]
   data_cell_clim = data_cell_clim[np.isnan(data_cell_clim) == 0]
   pct[:,i,j] = library_f90.percentileofscore(data_cell_clim,data[0,i,j],random.random())

 return pct

def Calculate_SPI(data_clim,data):

 #spi = np.zeros(data.shape)
 pct = np.zeros(data.shape)

 #Iterate through all cells and calculate the SPI values
 for i in range(data.shape[1]):
  #print i
  for j in range(data.shape[2]):
   data_clim_cell = data_clim[:,i,j]
   data_clim_cell = data_clim_cell[np.isnan(data_clim_cell) == 0]
   score = data[0,i,j]
   pct[:,i,j] = library_f90.percentileofscore(data_clim_cell,score,random.random())/100.0
    #pct.append(ss.percentileofscore(data_clim_cell,score,kind='mean')/100.0)
   #pct = np.array(pct)
   #pct[pct > 0.9999] = 0.9999
   #pct[pct < 0.0001] = 0.0001
   #spi[:,i,j] = ss.norm.ppf(pct)

 pct[pct > 0.9999] = 0.9999
 pct[pct < 0.0001] = 0.0001
 spi = ss.norm.ppf(pct)

 return spi

def CDF_Match(data_baseline_clim,data_biased_clim,data_biased):

 data_corrected = np.zeros(data_biased.shape)

 #Iterate through all cells and cdf match the data 
 for i in range(data_baseline_clim.shape[1]):
  for j in range(data_baseline_clim.shape[2]):
   #data_baseline_clim_cell = data_baseline_clim[:,i,j]
   #data_biased_clim_cell = data_biased_clim[:,i,j]
   data_biased_cell = data_biased[:,i,j]
   pct = []
   for score in data_biased_cell:
    #pct.append(ss.percentileofscore(data_biased_clim[:,i,j],score,kind='mean'))
    pct.append(library_f90.percentileofscore(data_biased_clim[:,i,j],score,random.random()))
   #score_unbiased = ss.scoreatpercentile(data_baseline_clim[:,i,j],pct)
   score_unbiased = library_f90.scoreatpercentile(data_baseline_clim[:,i,j],pct)
   data_corrected[:,i,j] = score_unbiased

 return data_corrected

def Regrid_and_Output_3B42rt(date,dims,Reprocess_Flag):

 print date
 file_out = '../DATA/3B42RT_RG/DAILY/3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return
 itime = datetime.datetime(2000,3,1)
 if date < itime:
  return

 file_in = '../DATA/3B42RT/DAILY/3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 ga("sdfopen %s" % file_in)

 #Regrid data
 ga("data = re(re(maskout(prec,prec),1.0),0.25)")
 data = ga.exp("data")

 #Output data
 nt = 1
 tstep = 'days'
 vars = ['prec']
 vars_info = ['daily total precip (mm) [Regridded to 1.0 deg and bilinear interpolation back to 0.25 deg]']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)
 #Write to file
 fp.variables['prec'][0] = data
 #Close files
 fp.close()
 ga("close 1")

 #Update the control file
 ctl_out = '../DATA/3B42RT_RG/DAILY/3B42RT_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('3B42RT','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 return

def BiasCorrect_and_Output_Forcing_GFS_Daily(date,dims,Reprocess_Flag):

 date = date + datetime.timedelta(days=1)
 dir = '../DATA/GFS_BC/%04d%02d%02d' % (date.year,date.month,date.day)
 
 #Reinitialize grads
 ga('reinit')

 #If the date is before the product's start date:
 if date < datetime.datetime(2012,8,14):
  return

 print_info_to_command_line("Bias-correcting the gfs 7-day forecast")

 #If the directory already exists exit:
 if os.path.exists(dir) == False:
  os.system("mkdir %s" % dir)

 dt = relativedelta.relativedelta(years=1)
 type = "all"

 #baseline forcing
 ctl_pgf = "../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl"

 vars = ["prec","tmax","tmin","wind"]
 vars_info = ['daily total precip (mm)','daily maximum temperature (K)','daily minimum temperature (K)','daily mean wind speed (m/s)']
 idate = date
 for t in range(7):
 
  #Set info for current day
  idate_pgf = date - relativedelta.relativedelta(years=date.year - 2001)#datetime.datetime(1950,date.month,date.day)
  fdate_pgf = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2008,date.month,date.day)
  idate_gfs = idate - relativedelta.relativedelta(years=date.year - 2012)#datetime.datetime(2001,date.month,date.day)
  fdate_gfs = idate - relativedelta.relativedelta(years=date.year - 2013)#datetime.datetime(2012,date.month,date.day)

  #original forcing
  ctl_gfs = "../DATA/GFS/gfs_daily_0.250deg_day%d.ctl" % (t+1)
  file = dir + '/gfs_%04d%02d%02d_daily_%.3fdeg_day%d.nc' % (idate.year,idate.month,idate.day,dims['res'],t+1)

  #Update control file
  date_ctl = datetime.datetime(2012,8,14)
  ndays = (idate - date_ctl).days + 1
  ctl_new = '../DATA/GFS_BC/gfs_daily_0.250deg_day%d.ctl' % (t+1)
  fp = open(ctl_new,'w')
  fp.write('dset ^%s%s%s/gfs_%s%s%s_daily_0.250deg_day%d.nc\n' % ('%y4','%m2','%d2','%y4','%m2','%d2',t+1))
  fp.write('options template\n')
  fp.write('dtype netcdf\n')
  fp.write('tdef t %d linear 14aug2012 1dy\n' % ndays)
  fp.close()

  #Determine if we skip the time step
  if os.path.exists(file) == True and Reprocess_Flag == False:
   date = date + datetime.timedelta(days=1)
   continue

  #Create file
  fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,'days',1)

  for var in vars:

   print var

   #Extract the required data
   print "Extracing gfs data (Climatology)"
   dt_up = relativedelta.relativedelta(days=20)
   dt_down = relativedelta.relativedelta(days=20)
   data_gfs_clim = Extract_Data_Period_Average(idate_gfs,fdate_gfs,dt_down,dt_up,dt,ctl_gfs,var,type,'xdfopen')
   print "Extracting pgf data (Climatology)"
   dt_up = relativedelta.relativedelta(days=3)
   dt_down = relativedelta.relativedelta(days=3)
   data_pgf_clim = Extract_Data_Period_Average(idate_pgf,fdate_pgf,dt_down,dt_up,dt,ctl_pgf,var,type,'xdfopen')
   print "Extracing gfs data (To correct)"
   dt_up = relativedelta.relativedelta(days=0)
   dt_down = relativedelta.relativedelta(days=0)
   data_gfs = Extract_Data_Period_Average(idate,idate,dt_down,dt_up,dt,ctl_gfs,var,type,'xdfopen')

   #CDF match the data
   print "Matching the daily gfs to the pgf"
   data_gfs_corrected = CDF_Match(data_pgf_clim,data_gfs_clim,data_gfs)

   #Write to file
   fp.variables[var][0] = data_gfs_corrected

  #Close file
  fp.close()

  #Update control file
  #date_ctl = datetime.datetime(2012,8,14)
  #ndays = (idate - date_ctl).days + 1
  #ctl_new = '../DATA/GFS_BC/gfs_daily_0.250deg_day%d.ctl' % (t+1)
  #fp = open(ctl_new,'w')
  #fp.write('dset ^%s%s%s/gfs_%s%s%s_daily_0.250deg_day%d.nc\n' % ('%y4','%m2','%d2','%y4','%m2','%d2',t+1))
  #fp.write('options template\n')
  #fp.write('dtype netcdf\n')
  #fp.write('tdef t %d linear 14aug2012 1dy\n' % ndays)
  #fp.close()

  #Update the time step
  date = date + datetime.timedelta(days=1)

 return

def BiasCorrect_and_Output_GFSANL_Daily(date,dims,Reprocess_Flag):

 #define parameters
 file_out = '../DATA/GFS_ANL_BC/DAILY/GFS_ANL_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])

 #Update the control file
 itime = datetime.datetime(2008,1,1)
 ctl_out = '../DATA/GFS_ANL_BC/DAILY/gfsanl_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('GFS_ANL','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the file exists then exit (unless we are reprocessing the data)
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return
 
 if date < itime:
  return

 dt = relativedelta.relativedelta(years=1)
 idate_pgf = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(1950,date.month,date.day)
 fdate_pgf = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2008,date.month,date.day)
 idate_anl = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2001,date.month,date.day)
 fdate_anl = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2012,date.month,date.day)
 type = "all"

 #original forcing
 ctl_anl = "../DATA/GFS_ANL/DAILY/gfsanl_daily_0.25deg.ctl"
 #baseline forcing
 ctl_pgf = "../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl"

 #Output data
 nt = 1
 tstep = 'days'
 #vars = ['prec']
 vars = ["tmax","tmin","wind"]
 vars_info = ['daily maximum temperature (K)','daily minimum temperature (K)','daily mean wind speed (m/s)']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)

 for var in vars:
 
  print var

  #Extract the required data
  print "Extracting pgf data (Climatology)"
  dt_up = relativedelta.relativedelta(days=30)
  dt_down = relativedelta.relativedelta(days=30)
  data_pgf_clim = Extract_Data_Period_Average(idate_pgf,fdate_pgf,dt_down,dt_up,dt,ctl_pgf,var,type,'xdfopen')
  print "Extracing anl data (Climatology)"
  dt_up = relativedelta.relativedelta(days=30)
  dt_down = relativedelta.relativedelta(days=30)
  data_anl_clim = Extract_Data_Period_Average(idate_anl,fdate_anl,dt_down,dt_up,dt,ctl_anl,var,type,'xdfopen')
  print "Extracing anl data (To correct)"
  dt_up = relativedelta.relativedelta(days=0)
  dt_down = relativedelta.relativedelta(days=0)
  data_anl = Extract_Data_Period_Average(date,date,dt_down,dt_up,dt,ctl_anl,var,type,'xdfopen')


  #CDF match the data
  print "Matching the daily anl to the pgf"
  data_anl_corrected = CDF_Match(data_pgf_clim,data_anl_clim,data_anl)

  #Write to file
  fp.variables[var][0] = data_anl_corrected

 #Close file
 fp.close()

 return

def BiasCorrect_and_Output_Forcing_FNL_Daily(date,dims,Reprocess_Flag):

 #define parameters
 file_out = '../DATA/FNL_ANALYSIS_BC/DAILY/fnlanl_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])

 #If the file exists then exit (unless we are reprocessing the data)
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return

 if date < datetime.datetime(2010,1,1):
  return

 dt = relativedelta.relativedelta(years=1)
 idate_pgf = date - relativedelta.relativedelta(years=date.year - 2006)#datetime.datetime(1950,date.month,date.day)
 fdate_pgf = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2008,date.month,date.day)
 idate_fnl = date - relativedelta.relativedelta(years=date.year - 2010)#datetime.datetime(2001,date.month,date.day)
 fdate_fnl = date - relativedelta.relativedelta(years=date.year - 2012)#datetime.datetime(2012,date.month,date.day)
 type = "all"

 #original forcing
 ctl_fnl = "../DATA/FNL_ANALYSIS/DAILY/fnlanl_daily_0.25deg.ctl"
 #baseline forcing
 ctl_pgf = "../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl"

 #Output data
 nt = 1
 tstep = 'days'
 #vars = ['prec']
 vars = ["prec","tmax","tmin","wind"]
 vars_info = ['daily total precip (mm)','daily maximum temperature (K)','daily minimum temperature (K)','daily mean wind speed (m/s)']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)

 for var in vars:
 
  print var

  #Extract the required data
  print "Extracting pgf data (Climatology)"
  dt_up = relativedelta.relativedelta(days=10)
  dt_down = relativedelta.relativedelta(days=10)
  data_pgf_clim = Extract_Data_Period_Average(idate_pgf,fdate_pgf,dt_down,dt_up,dt,ctl_pgf,var,type,'xdfopen')
  print "Extracing fnl data (Climatology)"
  dt_up = relativedelta.relativedelta(days=10)
  dt_down = relativedelta.relativedelta(days=10)
  data_fnl_clim = Extract_Data_Period_Average(idate_fnl,fdate_fnl,dt_down,dt_up,dt,ctl_fnl,var,type,'xdfopen')
  print "Extracing fnl data (To correct)"
  dt_up = relativedelta.relativedelta(days=0)
  dt_down = relativedelta.relativedelta(days=0)
  data_fnl = Extract_Data_Period_Average(date,date,dt_down,dt_up,dt,ctl_fnl,var,type,'xdfopen')


  #CDF match the data
  print "Matching the daily fnl to the pgf"
  data_fnl_corrected = CDF_Match(data_pgf_clim,data_fnl_clim,data_fnl)

  #Write to file
  fp.variables[var][0] = data_fnl_corrected

 #Close file
 fp.close()

 return
   
def BiasCorrect_and_Output_Forcing_3B42RT_Daily(date,dims,Reprocess_Flag):

 #define parameters
 file_out = '../DATA/3B42RT_BC/DAILY/3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])

 #Update the control file
 itime = datetime.datetime(2000,3,1)
 ctl_out = '../DATA/3B42RT_BC/DAILY/3B42RT_daily_0.25deg.ctl'
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('3B42RT','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the file exists then exit (unless we are reprocessing the data)
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return
 
 if date < itime:
  return

 dt = relativedelta.relativedelta(years=1)
 idate_pgf = date - relativedelta.relativedelta(years=date.year - 2001)#datetime.datetime(1950,date.month,date.day)
 fdate_pgf = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2008,date.month,date.day)
 idate_3b42rt = date - relativedelta.relativedelta(years=date.year - 2001)#datetime.datetime(2001,date.month,date.day)
 fdate_3b42rt = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2012,date.month,date.day)
 var = "prec"
 type = "all"

 #original precipitation
 ctl_3b42rt = "../DATA/3B42RT_RG/DAILY/3B42RT_daily_0.25deg.ctl"
 #baseline precipitation
 ctl_pgf = "../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl"

 #Extract the required data
 print "Extracting pgf data (Climatology)"
 dt_up = relativedelta.relativedelta(days=5)
 dt_down = relativedelta.relativedelta(days=5)
 data_pgf_clim = Extract_Data_Period_Average(idate_pgf,fdate_pgf,dt_down,dt_up,dt,ctl_pgf,var,type,'xdfopen')
 print "Extracing 3b42rt data (Climatology)"
 dt_up = relativedelta.relativedelta(days=5)
 dt_down = relativedelta.relativedelta(days=5)
 data_3b42rt_clim = Extract_Data_Period_Average(idate_3b42rt,fdate_3b42rt,dt_down,dt_up,dt,ctl_3b42rt,var,type,'xdfopen')
 print "Extracing 3b42rt data (To correct)"
 dt_up = relativedelta.relativedelta(days=0)
 dt_down = relativedelta.relativedelta(days=0)
 data_3b42rt = Extract_Data_Period_Average(date,date,dt_down,dt_up,dt,ctl_3b42rt,var,type,'xdfopen')

 #CDF match the data
 print "Matching the daily 3b42rt to the pgf"
 data_3b42rt_corrected = CDF_Match(data_pgf_clim,data_3b42rt_clim,data_3b42rt)

 #Output data
 nt = 1
 tstep = 'days'
 vars = ['prec']
 vars_info = ['daily total precip (mm) [cdf-matched against pgf]']
 #Create file
 #file_out = '../DATA/3B42RT_BC_SS/DAILY/3B42RT_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)
 #Write to file
 fp.variables['prec'][0] = data_3b42rt_corrected
 #Close file
 fp.close()

 return

def Calculate_and_Output_SPI(date,dims,dataset,idate,Reprocess_Flag):

 #define parameters
 itime = datetime.datetime(1950,1,1)
 if date >= itime and date <= datetime.datetime(2008,12,31):
  ctl_in = "../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl"
 elif date > datetime.datetime(2008,12,31):
  ctl_in = "../DATA/3B42RT_BC/DAILY/3B42RT_daily_0.25deg.ctl"
 else:
  return
 if dataset == 'forecast':
  ctl_in_forecast = "../DATA/GFS_DERIVED/%04d%02d%02d/output_grid_%04d%02d%02d00.ctl" % (idate.year,idate.month,idate.day,idate.year,idate.month,idate.day)

 #Define output and ctl files
 if dataset == 'monitor':
  file_out = '../DATA/SPI/DAILY/SPI_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/SPI/DAILY/spi_daily_0.25deg.ctl'
 elif dataset == 'forecast':
  itime = idate
  file_out = '../DATA/GFS_DERIVED/%04d%02d%02d/SPI_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/GFS_DERIVED/%04d%02d%02d/spi_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)

 #Precipitation dataset
 dt = relativedelta.relativedelta(years=1)
 idate_clim = date - relativedelta.relativedelta(years=date.year - 1950)#datetime.datetime(1950,date.month,date.day)
 fdate_clim = date - relativedelta.relativedelta(years=date.year - 2008)#datetime.datetime(2008,date.month,date.day)
 idate_tstep = date
 fdate_tstep = date
 var = "prec"
 dt_up = relativedelta.relativedelta(days=0)
 type = "ave"

 #Update the control file
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('SPI','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the product for today exists then exit
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return

 spi = []
 spi_months = [1,3,6,12]
 for spi_month in spi_months:
  
  print "SPI months: %d" % spi_month 
  dt_down = relativedelta.relativedelta(months=spi_month)

  #Extract the climatology and the desired data
  #If the climatology has already been saved open and if not save it
  file_climatology = '../DATA/SPI/DAILY_CLIMATOLOGY/%dmonth/%02d%02d_daily_%.3fdeg.pck' % (spi_month,date.month,date.day,dims['res'])
  if os.path.exists(file_climatology) == True: 
   print "Loading the climatology data from storage"
   data_clim = pickle.load(open(file_climatology,"rb"))
  else:
   print "Calculating the climatology"
   data_clim = Extract_Data_Period_Average(idate_clim,fdate_clim,dt_down,dt_up,dt,ctl_in,var,type,'xdfopen')
   print "Saving the climatology data to storage"
   pickle.dump(data_clim,open(file_climatology,'wb'))
  #Extract the period we are interested in 
  if dataset == 'monitor':
   data = Extract_Data_Period_Average(idate_tstep,fdate_tstep,dt_down,dt_up,dt,ctl_in,var,type,'xdfopen')
  if dataset == 'forecast':
   data_obs = Extract_Data_Period_Average(idate_tstep,fdate_tstep,dt_down,dt_up,dt,ctl_in,var,type,'xdfopen')
   data = Extract_Data_Period_Average(idate_tstep,fdate_tstep,dt_down,dt_up,dt,ctl_in_forecast,var,type,'open')
   data[-nt:,:,:] = data_obs[-nt:,:,:]

  #Calculate SPI
  print "Calculating the SPI"
  spi.append(Calculate_SPI(data_clim,data))

 #Create files and output all the data
 nt = 1
 tstep = 'days'
 vars = []
 vars_info = []
 for spi_month in spi_months:
  vars.append('spi%d' % spi_month)
  vars_info.append('%d month SPI' % spi_month)

 #Create file
 #file_out = '../DATA/SPI/DAILY/SPI_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)

 #Write to file
 for j in xrange(0,len(spi_months)):
  fp.variables[vars[j]][0] = spi[j]
  j = j + 1

 #Close file
 fp.close()
 
 return

def Calculate_and_Output_Streamflow_Percentiles(date,dims,Reprocess_Flag,dataset,idate):

 #define parameters
 root = '../DATA/ROUTING_VIC_DERIVED'
 itime = datetime.datetime(1950,1,1)
 if date >= itime and date <= datetime.datetime(2008,12,31):
  ctl_in = "../DATA/ROUTING_VIC_PGF/DAILY/Streamflow.ctl"
 elif date > datetime.datetime(2008,12,31):
  ctl_in = "../DATA/ROUTING_VIC_3B42RT/DAILY/Streamflow.ctl"
 else:
  return
 if dataset == 'forecast':
  ctl_in = "../DATA/GFS_DERIVED/%04d%02d%02d/Streamflow.ctl" % (idate.year,idate.month,idate.day)
 
 if dataset == 'monitor':
  file_out = '../DATA/ROUTING_VIC_DERIVED/DAILY/routing_vic_derived_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/ROUTING_VIC_DERIVED/DAILY/routing_vic_derived_daily_0.25deg.ctl'
 elif dataset == 'forecast':
  itime = idate
  file_out = '../DATA/GFS_DERIVED/%04d%02d%02d/routing_vic_derived_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/GFS_DERIVED/%04d%02d%02d/routing_vic_derived_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)
 dt = relativedelta.relativedelta(years=1)
 idate_clim = date - relativedelta.relativedelta(years=date.year - 1950)
 fdate_clim = date - relativedelta.relativedelta(years=date.year - 2008)
 type = "all"

 #Update the control file
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('routing_vic_derived','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the product for today exists then exit
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Computing streamflow percentiles")

 #Extract the required data
 var = 'flw'
 dt_up = datetime.timedelta(days=1)
 dt_down = datetime.timedelta(days=1)
 print "Extracing streamflow data (Climatology)"
 #If the climatology has already been saved open and if not save it
 file_climatology = '../DATA/ROUTING_VIC_DERIVED/DAILY_CLIMATOLOGY/%02d%02d_daily_%.3fdeg.pck' % (date.month,date.day,dims['res'])
 if os.path.exists(file_climatology) == True:
  print "Loading the climatology data from storage"
  data_clim = pickle.load(open(file_climatology,"rb"))
 else:
  print "Calculating the climatology"
  data_clim = Extract_Data_Period_Average(idate_clim,fdate_clim,dt_down,dt_up,dt,ctl_in,var,type,'open')
  print "Saving the climatology data to storage"
  pickle.dump(data_clim,open(file_climatology,'wb'))

 print "Extracing streamflow data (To calculate percentile)"
 dt_up = relativedelta.relativedelta(days=0)
 dt_down = relativedelta.relativedelta(days=0)
 data = Extract_Data_Period_Average(date,date,dt_down,dt_up,dt,ctl_in,var,type,'open')

 #Calculate the percentiles
 print "Calculating the percentiles"
 pct = Calculate_Percentiles(data_clim,data)

 #Output data
 nt = 1
 tstep = 'days'
 vars = ['flw_pct']
 vars_info = ['Streamflow percentiles']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)
 #Write to file
 fp.variables['flw_pct'][0] = pct
 #Close file
 fp.close()
 
 return

def Calculate_and_Output_SM_Percentiles(date,dims,dataset,idate,Reprocess_Flag):

 #define parameters
 root = '../DATA/VIC_DERIVED'
 itime = datetime.datetime(1950,1,1)

 #Define control files
 if date >= itime and date <= datetime.datetime(2008,12,31):
  ctl_in = "../DATA/VIC_PGF/DAILY/vic_daily_0.25deg.ctl"
 elif date > datetime.datetime(2008,12,31):
  ctl_in = "../DATA/VIC_3B42RT/DAILY/vic_daily_0.25deg.ctl"
 else:
  return
 if dataset == 'forecast':
  ctl_in = "../DATA/GFS_DERIVED/%04d%02d%02d/output_grid_%04d%02d%02d00.ctl" % (idate.year,idate.month,idate.day,idate.year,idate.month,idate.day)

 #Determine output files
 if dataset == 'monitor':
  file_out = '../DATA/VIC_DERIVED/DAILY/vic_derived_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/VIC_DERIVED/DAILY/vic_derived_daily_0.25deg.ctl'
 elif dataset == 'forecast':
  itime = idate
  file_out = '../DATA/GFS_DERIVED/%04d%02d%02d/vic_derived_%04d%02d%02d_daily_%.3fdeg.nc' % (idate.year,idate.month,idate.day,date.year,date.month,date.day,dims['res'])
  ctl_out = '../DATA/GFS_DERIVED/%04d%02d%02d/vic_derived_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)
 dt = relativedelta.relativedelta(years=1)

 #Define other parameters and files
 soil_ctl = "../DATA/VIC/INPUT/soil_Africa_0.25deg_calibrated_final.ctl"
 dt = relativedelta.relativedelta(years=1)
 idate_clim = date - relativedelta.relativedelta(years=date.year - 1950)
 fdate_clim = date - relativedelta.relativedelta(years=date.year - 2008)
 idate_tstep = date
 fdate_tstep = date
 dt_up = relativedelta.relativedelta(days=1)
 dt_down = relativedelta.relativedelta(days=1)
 type = "all"

 #Update the control file
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('vic_derived','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the product for today exists then exit
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Computing soil moisture percentiles")

 #Extract necessary soil information
 ga("open %s" % soil_ctl)
 depth1 = ga.exp("depth1")
 depth2 = ga.exp("depth2")
 porosity1 = 1 - ga.exp("bulk1")/ga.exp("density1")
 porosity2 = 1 - ga.exp("bulk2")/ga.exp("density2")
 satsm1 = 1000.0*depth1*porosity1
 satsm2 = 1000.0*depth2*porosity1
 satsm = satsm1 + satsm2
 ga("close 1")

 #Extract the climatology and the desired data
 #If the climatology has already been saved open and if not save it
 file_climatology = '../DATA/VIC_DERIVED/DAILY_CLIMATOLOGY/%02d%02d_daily_%.3fdeg.pck' % (date.month,date.day,dims['res'])
 if os.path.exists(file_climatology) == True:
  print "Loading the climatology data from storage"
  vc_clim = pickle.load(open(file_climatology,"rb"))
 else:
  print "Calculating the climatology"
  sm1_clim = Extract_Data_Period_Average(idate_clim,fdate_clim,dt_down,dt_up,dt,ctl_in,"sm1",type,"open")
  sm2_clim = Extract_Data_Period_Average(idate_clim,fdate_clim,dt_down,dt_up,dt,ctl_in,"sm2",type,"open")
  vc_clim = 100.0*(sm1_clim + sm2_clim)/(satsm)

  print "Saving the climatology data to storage"
  pickle.dump(vc_clim,open(file_climatology,'wb'))

 #Extract the period we are interested in 
 dt_up = relativedelta.relativedelta(days=0)
 dt_down = relativedelta.relativedelta(days=0)
 sm1 = Extract_Data_Period_Average(idate_tstep,fdate_tstep,dt_down,dt_up,dt,ctl_in,"sm1",type,"open")
 sm2 = Extract_Data_Period_Average(idate_tstep,fdate_tstep,dt_down,dt_up,dt,ctl_in,"sm2",type,"open")
 vc1 = 100.0*sm1/satsm1
 vc2 = 100.0*sm2/satsm2
 vc = 100.0*(sm1 + sm2)/(satsm)
 vc[sm1 < 0] = -9.99e+08

 #Calculate the percentiles
 print "Calculating the percentiles"
 vcpct = Calculate_Percentiles(np.ma.getdata(vc_clim),np.ma.getdata(vc))

 #Output data
 nt = 1
 tstep = 'days'
 vars = ['vc1','vc2','vc','vcpct']
 vars_info = ['Volumetric water content (layer 1)','Volumetric water content (layer 2)','Volumetric water content (layers 1+2)','Volumetric water content percentile (layers 1+2)']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)
 #Write to file
 fp.variables['vc1'][0] = vc1
 fp.variables['vc2'][0] = vc2
 fp.variables['vc'][0] = vc
 fp.variables['vcpct'][0] = vcpct
 #Close file
 fp.close()

 return

def Update_Control_File(type,idate,dims,nt,tstep,file_template,ctl_file):

 if type == 'nc':
  fp = open(ctl_file,'w')
  fp.write('dset %s\n' % file_template)
  fp.write('options template\n')
  fp.write('dtype netcdf\n')
  fp.write('tdef t %d linear %s %s\n' % (nt,datetime2gradstime(idate),tstep))
  fp.close()

def Calculate_and_Output_NDVI_Percentiles(date,dims,Reprocess_Flag):

 file_org = '../DATA/MOD09_NDVI_MA/DAILY/MOD09CMG_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 if os.path.exists(file_org) == False:
  return

 #define parameters
 file_out = '../DATA/MOD09_NDVI_MA_DERIVED/DAILY/MOD09CMG_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])

 itime = datetime.datetime(2003,1,1)
 if date < itime:
  return

 #Update the control file
 ctl_out = '../DATA/MOD09_NDVI_MA_DERIVED/DAILY/MOD09CMG_daily_0.25deg.ctl' 
 nt = (date - itime).days + 1
 file_template = '^%s_%s%s%s_daily_0.250deg.nc' % ('MOD09CMG','%y4','%m2','%d2')
 Update_Control_File('nc',itime,dims,nt,'1dy',file_template,ctl_out)

 #If the file exists then exit (unless we are reprocessing the data)
 if os.path.exists(file_out) == True and Reprocess_Flag == False:
  return

 print_info_to_command_line("Computing the NDVI percentiles")

 dt = relativedelta.relativedelta(years=1)
 idate_modis = date - relativedelta.relativedelta(years=date.year - 2003)
 fdate_modis = date - relativedelta.relativedelta(years=date.year - 2012)
 var = "ndvi30"
 type = "all"

 #moving average ndvi product
 ctl_modis = "../DATA/MOD09_NDVI_MA/DAILY/MOD09CMG_daily_0.25deg.ctl"

 #Extract the required data
 print "Extracing modis data (Climatology)"
 dt_up = relativedelta.relativedelta(days=5)
 dt_down = relativedelta.relativedelta(days=5)
 data_modis_clim = Extract_Data_Period_Average(idate_modis,fdate_modis,dt_down,dt_up,dt,ctl_modis,var,type,'xdfopen')
 print "Extracing modis data (To calculate percentile)"
 dt_up = relativedelta.relativedelta(days=0)
 dt_down = relativedelta.relativedelta(days=0)
 data_modis = Extract_Data_Period_Average(date,date,dt_down,dt_up,dt,ctl_modis,var,type,'xdfopen')

 #Calculate the percentiles
 print "Calculating the percentiles"
 pct = Calculate_Percentiles(data_modis_clim,data_modis)

 #Output data
 nt = 1
 tstep = 'days'
 vars = ['pct30day']
 vars_info = ['NDVI percentiles']
 #Create file
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,date,tstep,nt)
 #Write to file
 fp.variables['pct30day'][0] = pct
 #Close file
 fp.close()

 return

def datetime2gradstime(date):

 #Convert datetime to grads time
 str = date.strftime('%HZ%d%b%Y')

 return str

def gradstime2datetime(str):

 #Convert grads time to datetime
 date = datetime.datetime.strptime(str,'%HZ%d%b%Y')

 return date

def Compute_Time_Averages(ctl_in,file_out,idate,fdate,dims,open_type):

 #Open file in grads
 ga("%s %s" % (open_type,ctl_in))

 #What variables exist in the dataset?
 qh = ga.query("file")

 #Create and open access to netcdf file
 nt = 1
 tstep = 'hours'
 fp = Create_NETCDF_File(dims,file_out,qh.vars,qh.var_info,idate,tstep,nt) #Create and open access to netcdf file

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
   #ga("data = %s(maskout(%s,%s),time=%s,time=%s)" % (type,var,var,t1,t2))
   ga("data = %s(%s,time=%s,time=%s)" % (type,var,t1,t2))
   ga("data = const(data,-9.99e+08,-u)")

   #Write to file
   data = np.ma.getdata(ga.exp("data"))
   fp.variables[var][0] = data

 #Close access to the grads file
 ga("close 1")

def Download_and_Process_Seasonal_Forecast(date,Reprocess_Flag):

 dir0 = '../DATA/SEASONAL_FORECAST/%04d%02d' % (date.year,date.month)

 #If the date is before the product's start date:
 if date < datetime.datetime(2013,1,1):
  return

 #If the directory already exists exit
 if os.path.exists(dir0) == True and Reprocess_Flag == False:
  return

 #If it before the 15th do not attempt to download
 if date.day < 15:
  return

 print_info_to_command_line("Downloading the monthly seasonal forecast")

 models = ('CMC1-CanCM3','CMC2-CanCM4','COLA-RSMAS-CCSM3','GFDL-CM2p1-aer04','NASA-GMAO-062012')
 nensembles = (10,10,6,10,12)
 type = ('.FORECAST/','.FORECAST/','','','')
 vars = ('prec','tref')
 month = date.strftime('%b') 
 year = date.year
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
    os.system('wget -nv -O {0} {1}'.format(file_out,http_file))
   i = i + 1

 #Download CFS-V2
 models = ('CFSv2',)
 nensembles = (24,)
 vars = {'PRATE':'surface','TMP':'2m_above_ground'}
 idate = datetime.datetime(date.year,date.month,1) - datetime.timedelta(days=1)
 if not os.path.exists(dir0):
  os.makedirs(dir0)
 for var in vars:
  dir = dir0 + '/' + var
  if not os.path.exists(dir):
   os.makedirs(dir)
  i = 0 
  for model in models:
   for ens in xrange(1,nensembles[i]+1):
    date = idate - (ens-1)*datetime.timedelta(days=1)
    print ens,date
    root = ''
    http_file = 'http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.EMC/.CFSv2/.MONTHLY_REALTIME/.FLXF/.{3}/.{4}/M/%280%29%280%29RANGEEDGES/S/%280000%20{0}%20{1}%20{2}%29VALUES/-9.9900000E08/setmissing_value/data.cdf'.format(date.day,date.strftime('%b'),date.year,vars[var],var)
    file_out = '{0}/{1}_{2}.nc'.format(dir,model,ens)
    os.system('wget -nv -O {0} {1}'.format(file_out,http_file))
    os.system('ncwa -O -a S,M %s %s' % (file_out,file_out))
   i = i + 1

 #Regrid and extract only months 2-7 for each ensemble member
 dims_globe = {}
 dims_globe['minlat'] = -90.0 #-89.8750
 dims_globe['minlon'] = 0.0 #0.1250
 dims_globe['nlat'] = 181 #720
 dims_globe['nlon'] = 360 #1440
 dims_globe['res'] = 1.0
 dims_globe['maxlat'] = dims_globe['minlat'] + dims_globe['res']*(dims_globe['nlat']-1)
 dims_globe['maxlon'] = dims_globe['minlon'] + dims_globe['res']*(dims_globe['nlon']-1)
 nt = 6
 tstep = 'months'
 models = {'CFSv2':'NCEP-CFSv2',}
 vars = {'PRATE':'prec','TMP':'tref'}
 for var in vars:
  dir_old = dir0 + '/' + var
  dir_new = dir0 + '/' + vars[var]
  for model in models:
   for ens in xrange(1,24+1):
    print model,var,ens
    file_out_org = '{0}/{1}_{2}.nc'.format(dir_old,model,ens)
    file_out_new = '{0}/{1}_{2}.nc'.format(dir_new,models[model],ens)
    variables = [vars[var],]
    vars_info = [vars[var],]
    #Create file
    fp = Create_NETCDF_File(dims_globe,file_out_new,variables,vars_info,date,tstep,nt)

    #Open file in grads
    ga("sdfopen %s" % file_out_org)
    #Set to new region
    ga("set lat %f %f" % (dims_globe['minlat'],dims_globe['maxlat']))
    ga("set lon %f %f" % (dims_globe['minlon'],dims_globe['maxlon']))
    #Regrid and export to new file
    for t in xrange(2,7+1):
     ga("set t %d" % t)
     if var == 'PRATE':
      ga("data = 24*3600*PRATE")
      Grads_Regrid("data","prec",dims_globe)
     elif var == 'TMP':
      ga("data = TMP")
      Grads_Regrid("data","tref",dims_globe)
     data = ga.exp(vars[var])
     fp.variables[vars[var]][t-2] = data
    #Close file
    fp.close()
    #Close grads file 
    ga("close 1")

 return

def BiasCorrect_and_Compute_Seasonal_Forecast_Products(date,dims,Reprocess_Flag):

 models = {'CMC1-CanCM3':10,'CMC2-CanCM4':10,'COLA-RSMAS-CCSM3':6,'GFDL-CM2p1-aer04':10,'NASA-GMAO-062012':11,'NCEP-CFSv2':24}
 #models = {'CMC1-CanCM3':1,'CMC2-CanCM4':1,'COLA-RSMAS-CCSM3':1,'GFDL-CM2p1-aer04':1,'NASA-GMAO-062012':1}
 file = '../DATA/SEASONAL_FORECAST/%04d%02d/CMC1-CanCM3_10_sf.nc' % (date.year,date.month)

 #If the date is before the product's start date:
 if date < datetime.datetime(2013,1,1):
  return

 #If the directory already exists exit
 if os.path.exists(file) == True and Reprocess_Flag == False:
  return

 #If it before the 15th do not attempt to download
 if date.day < 15:
  return

 #Compute anomalies and SPI
 print_info_to_command_line("Bias-correcting and computing temperature anomalies for the Seasonal forecast")

 #Change directory to run the program
 owd = os.getcwd()
 os.chdir('../DATA/SEASONAL_FORECAST/code')

 #Write the namelist
 fp = open('namelist','w')
 fp.write('&nmmepar\n')
 fp.write('wlon = %d\n' % dims['nlon'])
 fp.write('wlat = %d\n' % dims['nlat'])
 fp.write('slon = %.3f\n' % dims['minlon'])
 fp.write('slat = %.3f\n' % dims['minlat'])
 fp.write('res  = %.2f\n' % dims['res'])
 fp.write('yr   = %d\n' % date.year)
 fp.write('nm   = %d\n' % date.month)
 fp.write("ext  = '../'\n")
 fp.write("ext0 = '../../3B42RT_BC/MONTHLY//'\n")
 fp.write('/\n')
 fp.close()

 #Run program
 os.system('./NMME-SPI.sh >& ../WORKSPACE/NMME-SPI_log.txt')

 #Change back to original directory
 os.chdir(owd)
 
 OUTPUT = {}
 vars = {'pr.1':'prec','t2m.1':'t2m','spi1.2':'spi1','spi3.2':'spi3','spi6.2':'spi6','spi12.2':'spi12','t2ano.2':'t2ano'}
 #Extract the mean of all the ensembles and models
 for model in models:
  OUTPUT[model] = {}
  nens = models[model]
  for ensemble in xrange(1,nens+1):
   print model,ensemble
   #Open bias corrected data
   ctl_bc = '../DATA/SEASONAL_FORECAST/%04d%02d/%s_%d_BC.nc' % (date.year,date.month,model,ensemble)
   ctl_sf = '../DATA/SEASONAL_FORECAST/%04d%02d/%s_%d_sf.nc' % (date.year,date.month,model,ensemble)
   ga("sdfopen %s" % ctl_bc)
   ga("sdfopen %s" % ctl_sf)
   ga("set t 1 6")
   for var in vars:
    if ensemble == 1:
     OUTPUT[model][var] = ga.exp(var)
    else:
     OUTPUT[model][var] = OUTPUT[model][var] + ga.exp(var)
   #Close file
   ga("close 2")
   ga("close 1")
  #Calculate the ensemble mean
  for var in vars:
   OUTPUT[model][var] = OUTPUT[model][var]/nens

 #Calculate the multi model mean
 OUTPUT['MultiModel'] = {}
 nmodels = 0
 for model in models:
  nmodels = nmodels + 1
  for var in vars:
   if nmodels == 1:
    OUTPUT['MultiModel'][var] = OUTPUT[model][var]
   else:
    OUTPUT['MultiModel'][var] = OUTPUT['MultiModel'][var] + OUTPUT[model][var]
 for var in vars:
  OUTPUT['MultiModel'][var] = OUTPUT['MultiModel'][var]/nmodels
 
 #Output data
 variables = ['prec','t2m','spi1','spi3','spi6','spi12','t2ano']
 vars_info = ['Precipitation (mm)','Temperature (C)','SPI 1-month','SPI 3-month','SPI 6-month','SPI 12-month','Temperature anomaly (C)']
 for model in OUTPUT:
  file = '../DATA/SEASONAL_FORECAST/%04d%02d/%s_%04d%02d_monthly.nc' % (date.year,date.month,model,date.year,date.month)
  nt = 6
  tstep = 'months'
  #Create file
  fp = Create_NETCDF_File(dims,file,variables,vars_info,date,tstep,nt)
  #Write to file
  for t in xrange(0,6):
   for var in vars:
    fp.variables[vars[var]][t] = OUTPUT[model][var][t,:,:]
  #Close file
  fp.close()

  #Create control file
  idate = datetime.datetime(2013,1,1)
  fdate = date
  nmonths = 12*(fdate.year - idate.year) + max(fdate.month-idate.month,0) + 1
  ctl = '../DATA/SEASONAL_FORECAST/CTL/%s.ctl' % model
  fp = open(ctl,'w')
  fp.write('dset ^../%s/%s_%s_monthly.nc\n' % ('%e',model,'%e'))
  fp.write('options template\n')
  fp.write('dtype netcdf\n')
  fp.write('title Seasonal Forecast %s\n' % model)
  fp.write('undef -9.99e+08\n')
  fp.write('xdef %d  linear %f %f\n' % (dims['nlon'],dims['minlon'],dims['res']))
  fp.write('ydef %d  linear %f %f\n' % (dims['nlat'],dims['minlat'],dims['res']))
  fp.write('tdef %d linear 00Z01%s%d 1mo\n' % (nmonths+5,idate.strftime('%b'),idate.year))
  fp.write('zdef 1 linear 1 1\n')
  fp.write('edef %d\n' % nmonths)
  date_tmp = idate
  while date_tmp <= fdate:
   fp.write('%04d%02d 6 00Z01%s%d 1mo\n' % (date_tmp.year,date_tmp.month,date_tmp.strftime('%b'),date_tmp.year))
   date_tmp = date_tmp + relativedelta.relativedelta(months=1)
  fp.write('endedef\n')
  fp.write('vars 7\n')
  fp.write('spi1 1 t,y,x data\n')
  fp.write('spi3 1 t,y,x data\n')
  fp.write('spi6 1 t,y,x data\n')
  fp.write('spi12 1 t,y,x data\n')
  fp.write('prec 1 t,y,x data\n')
  fp.write('t2ano 1 t,y,x data\n')
  fp.write('t2m 1 t,y,x data\n')
  fp.write('endvars\n')
  fp.close()

 return

def print_info_to_command_line(line):

 print "\n"
 print "#######################################################################################" 
 print "%s" % line
 print "#######################################################################################"
 print "\n"

 return

def Reprocess_PGF(date,dims): 

 pgf_new_root = '../DATA/PGF/DAILY'

 #If the date is after the product's final date return:
 if date > datetime.datetime(2008,12,31):
  return

 #If the file already exists exit:
 file = 'pgf_%04d%02d%02d_daily_%.3fdeg.nc' % (date.year,date.month,date.day,dims['res'])
 netcdf_file = '%s/%s' % (pgf_new_root,file)
 if os.path.exists(netcdf_file) == True:
  return

 #Load datasets into grads
 prec_ctl = '../DATA/PGF/ORIGINAL/prec/prec_1948-2008.ctl'
 tmax_ctl = '../DATA/PGF/ORIGINAL/tmax/tmax_1948-2008.ctl'
 tmin_ctl = '../DATA/PGF/ORIGINAL/tmin/tmin_1948-2008.ctl'
 wind_ctl = '../DATA/PGF/ORIGINAL/wind/wind_1948-2008.ctl'
 ga("open %s" % prec_ctl)
 ga("open %s" % tmax_ctl)
 ga("open %s" % tmin_ctl)
 ga("open %s" % wind_ctl)

 #Set time
 time = datetime2gradstime(date)
 ga("set time %s" % time)

 #Set to new region
 ga("set lat %f %f" % (dims['minlat'],dims['maxlat']))
 ga("set lon %f %f" % (dims['minlon'],dims['maxlon']))

 #Create new data
 Grads_Regrid("tmax.2","tmax",dims)
 Grads_Regrid("tmin.3","tmin",dims)
 Grads_Regrid("prec.1","prec",dims)
 Grads_Regrid("wind.4","wind",dims)

 #Create and open access to netcdf file
 vars = ['tmax','tmin','prec','wind']
 vars_info = ['daily tmax (K)','daily tmin (K)','daily total precip (mm)','daily mean wind speed (m/s)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,netcdf_file,vars,vars_info,date,tstep,nt)

 #Write to file
 for var in vars:
  data = np.ma.getdata(ga.exp(var))
  fp.variables[var][0] = data

 #Finalize and close NETCDF file
 fp.close()

 #Close access to all files in grads
 ga("close 4")
 ga("close 3")
 ga("close 2")
 ga("close 1")

 return

def Collect_Dataset_Bounds(dataset,ctl_in,info,open_type):

 #Determine dataset bounds
 ga("%s %s" % (open_type,ctl_in))
 qh = ga.query("file")
 vars = qh.vars

 #Iterate through the variables
 ga("set t 1")
 itime = gradstime2datetime(ga.exp(vars[0]).grid.time[0])
 ga("set t last")
 ftime = gradstime2datetime(ga.exp(vars[0]).grid.time[0])
 ga("close 1")

 #Place info
 info['itime'] = itime
 info['ftime'] = ftime

 return info

def Compute_Monthly_Yearly_Averages(date,dims,dt,dataset,ctl_in,open_type,reprocess_flag,itime,ftime):

 #Update the control files (MONTHLY)
 ctl_out = '../DATA/%s/MONTHLY/%s_monthly_%.3fdeg.ctl' % (dataset,dataset,dims['res'])
 nt = 12*(ftime.year - itime.year) + max(ftime.month - itime.month,0)# + 1
 file_template = '^%s_%s%s_monthly_%.3fdeg.nc' % (dataset,'%y4','%m2',dims['res'])
 Update_Control_File('nc',itime,dims,nt,'1mo',file_template,ctl_out)

 #Update the control files (YEARLY)
 ctl_out = '../DATA/%s/YEARLY/%s_yearly_%.3fdeg.ctl' % (dataset,dataset,dims['res'])
 nt = max(ftime.year - itime.year,0)# + 1
 file_template = '^%s_%s_yearly_%.3fdeg.nc' % (dataset,'%y4',dims['res'])
 Update_Control_File('nc',itime,dims,nt,'1yr',file_template,ctl_out)

 #Check for new month or new year
 ndate = date + dt
 if date.month == ndate.month:
  return

 #If we are before or after the datasets last time step, exit
 if date < itime or date > ftime:
  return

 #Determine dataset bounds
 #ga("%s %s" % (open_type,ctl_in))
 #qh = ga.query("file")
 #vars = qh.vars

 #Iterate through the variables
 #ga("set t 1")
 #itime = gradstime2datetime(ga.exp(vars[0]).grid.time[0])
 #ga("set t last")
 #ftime = gradstime2datetime(ga.exp(vars[0]).grid.time[0])
 #ga("close 1")
 
 yearly_dir = "../DATA/%s/YEARLY" % dataset
 monthly_dir = "../DATA/%s/MONTHLY" % dataset

 #Check for new month
 ndate = date + dt
 idate = datetime.datetime(date.year,date.month,1)
 fdate = date
 if date.month == ndate.month:
  return

 #If the file does not exist compute the mean
 file_out = "../DATA/%s/MONTHLY/%s_%04d%02d_monthly_%.3fdeg.nc" % (dataset,dataset,idate.year,idate.month,dims['res'])
 if os.path.exists(file_out) == False or reprocess_flag == True:

  if os.path.exists(monthly_dir) == False:
   os.mkdir(monthly_dir)
  print_info_to_command_line("Computing %s monthly time averages" % dataset)

  #Compute Monthly Average
  Compute_Time_Averages(ctl_in,file_out,idate,fdate,dims,open_type)
  
 #Update the control file
 ctl_out = '../DATA/%s/MONTHLY/%s_monthly_%.3fdeg.ctl' % (dataset,dataset,dims['res'])
 nt = 12*(date.year - itime.year) + max(date.month - itime.month,0) + 1
 file_template = '^%s_%s%s_monthly_%.3fdeg.nc' % (dataset,'%y4','%m2',dims['res'])
 Update_Control_File('nc',itime,dims,nt,'1mo',file_template,ctl_out)

 #Check for new year
 if date.year == ndate.year:
  return

 idate = datetime.datetime(date.year,1,1)
 fdate = date

 #If the file already exists exit:
 file_out = "../DATA/%s/YEARLY/%s_%04d_yearly_%.3fdeg.nc" % (dataset,dataset,idate.year,dims['res'])
 if os.path.exists(file_out) == False or reprocess_flag == True:

  if os.path.exists(yearly_dir) == False:
   os.mkdir(yearly_dir)
  print_info_to_command_line("Computing %s yearly time averages" % dataset)

  #Compute Annual Average
  Compute_Time_Averages(ctl_in,file_out,idate,fdate,dims,open_type)

 #Update the control file
 ctl_out = '../DATA/%s/YEARLY/%s_yearly_%.3fdeg.ctl' % (dataset,dataset,dims['res'])
 nt = max(date.year - itime.year,0) + 1
 file_template = '^%s_%s_yearly_%.3fdeg.nc' % (dataset,'%y4',dims['res'])
 Update_Control_File('nc',itime,dims,nt,'1yr',file_template,ctl_out)

 return

def Prepare_VIC_Global_Parameter_File(idate,fdate,dims,dataset):

 file = '../WORKSPACE/Global_Parameter.txt'
 fp = open(file,'w')
 fdate = fdate + datetime.timedelta(days=1)

 #Write the VIC parameters to file

 # Define Global Parameters
 fp.write('NLAYER          3       # number of layers\n')
 fp.write('TIME_STEP       24       # model time step in hours (= 24 for water balance)\n')
 fp.write('STARTYEAR       %d      # year model simulation starts\n' % idate.year)
 fp.write('STARTMONTH      %d      # month model simulation starts\n' % idate.month)
 fp.write('STARTDAY        %d      # day model simulation starts\n' % idate.day)
 fp.write('STARTHOUR       0       # hour model simulation starts\n')
 fp.write('ENDYEAR         %d      # year model simulation ends\n' % fdate.year)
 fp.write('ENDMONTH        %d      # month model simulation ends\n' % fdate.month) 
 fp.write('ENDDAY          %d      # day model simulation ends\n' % fdate.day) 
 fp.write('SKIPYEAR        0       # no. of startup yrs to skip before writing output\n')
 fp.write('WIND_H          10.0    # height of wind speed measurement\n')
 fp.write('MEASURE_H       2.0     # height of humidity measurement\n')
 fp.write('NODES           10       # number of soil thermal nodes\n')
 fp.write('MAX_SNOW_TEMP   0.5     # maximum temperature at which snow can fall\n')
 fp.write('MIN_RAIN_TEMP   -0.5    # minimum temperature at which rain can fall\n')

 # Define Global Parameters
 fp.write('FULL_ENERGY     FALSE    # calculate full energy balance\n')
 #fp.write('FROZEN_SOIL     TRUE    # calculate frozen soils\n')
 fp.write('DIST_PRCP       TRUE        # use distributed precipitation\n')
 fp.write('COMPRESS        FALSE       # compress input and output files when done\n')
 fp.write('CORRPREC        FALSE       # correct precipitation for gauge undercatch\n')
 fp.write('GRID_DECIMAL    3           # number of decimals to use in gridded file names\n')
 fp.write('PRT_SNOW_BAND   FALSE   # print snow variables\n')
 fp.write('SNOW_STEP       1        # time step in hours to solve snow bands\n')
 fp.write('ROOT_ZONES      2               # number of root zones in veg parameter file\n')
 fp.write('BINARY_OUTPUT   TRUE   # default is ASCII, unless LDAS format\n')
 fp.write('MIN_WIND_SPEED  0.1     # minimum allowable wind speed\n')
 fp.write('PREC_EXPT       0.6             # fraction of grid cell receiving\n')
 fp.write('GRND_FLUX       FALSE # true for full energy, false for water balance\n')
 fp.write('QUICK_FLUX      FALSE   # true uses Liang (1999), false uses finite diff.\n')
 fp.write('NOFLUX          FALSE  # false uses const. T at damping depth\n')

 # Define (Meteorological) Forcing Files
 fp.write('FORCING1        ../DATA/VIC/FORCING/DAILY/forcing_daily_\n')
 fp.write('N_TYPES         4\n')
 fp.write('FORCE_TYPE      TMAX    SIGNED  10\n')
 fp.write('FORCE_TYPE      TMIN    SIGNED  10\n')
 fp.write('FORCE_TYPE      WIND    UNSIGNED 10\n')
 fp.write('FORCE_TYPE      PREC    UNSIGNED 10\n')
 fp.write('FORCE_FORMAT    BINARY \n')
 fp.write('FORCE_ENDIAN    LITTLE      # LITTLE for PC arch., BIG for Sun or HP-UX\n')
 fp.write('FORCE_DT        24            # time step of two input met files\n')
 fp.write('FORCEYEAR   %d   # year model meteorological forcing files start\n' % idate.year)
 fp.write('FORCEMONTH   %d   # month model meteorological forcing files start\n' % idate.month)
 fp.write('FORCEDAY        %d                   # day meteorological forcing files start\n' % idate.day)
 fp.write('FORCEHOUR       00                  # hour meteorological forcing files start\n')


 # INPUT and OUTPUT TYPE from PRINCETON  (mpan and lluo)
 fp.write('INPUT_GRID_DEF %d %d %.3f %.3f %.3f %.3f\n' % (dims['nlon'],dims['nlat'],dims['minlon'],dims['minlat'],dims['res'],dims['res']))
 fp.write('OUTPUT_GRID_DEF %d %d %.3f %.3f %.3f %.3f\n' % (dims['nlon'],dims['nlat'],dims['minlon'],dims['minlat'],dims['res'],dims['res']))
 #fp.write('INPUT_STEP_PER_FILE     1        # number of timesteps per input file\n')
 fp.write('OUTPUT_PER_STEP TRUE # number of timesteps per output file\n')
 #fp.write('INPUT_GZIP              FALSE       # true if forcing file gzipped\n')
 #fp.write('OUTPUT_GZIP             FALSE        # true for writing gzipped output\n')
 fp.write('GRID_INPUT          TRUE #true for reading the input in GrADS binary, default is false\n')
 fp.write('GRID_OUTPUT         TRUE  #true for writing the output in GrADS binary,default is false\n')
 fp.write('REGULAR_OUTPUT      FALSE  #true for writing the output in standard version, default is false\n')



 # Define Input and Output Data Files
 fp.write('SNOW_BAND       ../DATA/VIC/INPUT/global_lai_0.25deg.txt\n')
 fp.write('ARC_SOIL        FALSE   # read soil parameters from ARC/INFO ASCII grids\n')
 #fp.write('SOIL            ../DATA/VIC/INPUT/tmp.txt\n')
 fp.write('SOIL            ../DATA/VIC/INPUT/soil_Africa_0.25deg_calibrated_final.txt\n')
 fp.write('VEGPARAM        ../DATA/VIC/INPUT/global_lai_0.25deg.txt\n')
 fp.write('VEGLIB          ../DATA/VIC/INPUT/veglib.dat\n')
 fp.write('GLOBAL_LAI      TRUE      # true if veg param file has monthly LAI\n')
 if dataset == 'pgf':
  fp.write('RESULT_DIR      ../DATA/VIC_PGF/DAILY/\n')
 if dataset == '3b42rt':
  fp.write('RESULT_DIR      ../DATA/VIC_3B42RT/DAILY/\n')
 if dataset == 'gfsanl':
  fp.write('RESULT_DIR      ../DATA/VIC_GFSANL/DAILY/\n')
 if dataset == 'gfs_forecast':
  fp.write('RESULT_DIR      ../DATA/GFS_DERIVED/%04d%02d%02d/\n' % (idate.year,idate.month,idate.day))

 # Define the state file
 if dataset == 'pgf':
  fp.write('STATENAME ../DATA/VIC/STATE_PGF/state\n')
 if dataset == '3b42rt':
  fp.write('STATENAME ../DATA/VIC/STATE_3B42RT/state\n')
 if dataset == 'gfsanl':
  fp.write('STATENAME ../DATA/VIC/STATE_GFSANL/state\n')
 if dataset == 'gfsforecast':
  fp.write('STATENAME ../DATA/GFS_DERIVED/%04d%02d%02d/\n' % (idate.year,idate.month,idate.day))
 fp.write('STATEYEAR %d\n' % fdate.year)
 fp.write('STATEMONTH %d\n' % fdate.month)
 fp.write('STATEDAY %d\n' % fdate.day)
 if dataset == 'pgf':
  file_state = '../DATA/VIC/STATE_PGF/state_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 if dataset == '3b42rt' or dataset == 'gfs_forecast':
  file_state = '../DATA/VIC/STATE_3B42RT/state_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 if dataset == 'gfsanl':
  file_state = '../DATA/VIC/STATE_GFSANL/state_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 if os.path.exists(file_state) == True:
  fp.write('INIT_STATE %s\n' % file_state)

 #Close the file
 fp.close()

def Prepare_VIC_Forcings_GFS_forecast(idate,fdate,dims):

 #Forcing_Filename
 forcing_file = '../DATA/VIC/FORCING/DAILY/forcing_daily_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 fp = open(forcing_file,'wb')

 #Extract and print data
 date = idate
 itime = datetime2gradstime(date)
 dt = datetime.timedelta(days=1)
 iday = 0
 while date <= fdate:

  print date
  #Open the control file
  iday = iday + 1
  ctl_file = '../DATA/GFS_BC/gfs_daily_0.250deg_day%d.ctl' % iday
  ga('xdfopen %s' % ctl_file)

  #Set time
  ga("set time %s" % itime)
  prec = np.ma.getdata(ga.exp("prec"))
  tmax = np.ma.getdata(ga.exp("tmax-273.15"))
  tmin = np.ma.getdata(ga.exp("tmin-273.15"))
  wind = np.ma.getdata(ga.exp("wind"))

  #Append to the outgoing file
  prec.tofile(fp)
  tmax.tofile(fp)
  tmin.tofile(fp)
  wind.tofile(fp)

  #Update the time step
  date = date + dt

  #Close current gfs forecast day
  ga("close 1")

 #Close the outgoing file
 fp.close()

 return forcing_file

 return

def Prepare_VIC_Forcings_3B42RT(idate,fdate,dims):
  
 #Open files
 pgf_ctl = '../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl' 
 tmpabc_ctl = '../DATA/3B42RT_BC/DAILY/3B42RT_daily_0.25deg.ctl'
 gfsanl_ctl = '../DATA/GFS_ANL_BC/DAILY/gfsanl_daily_0.25deg.ctl'
 ga("xdfopen %s" % pgf_ctl)
 ga("xdfopen %s" % gfsanl_ctl)
 ga("xdfopen %s" % tmpabc_ctl)

 #Forcing_Filename
 forcing_file = '../DATA/VIC/FORCING/DAILY/forcing_daily_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 fp = open(forcing_file,'wb')

 #Extract and print data
 date = idate
 dt = datetime.timedelta(days=1)
 while date <= fdate:
  
  print date
  time = datetime2gradstime(date)
  ga("set time %s" % time)
  if date < datetime.datetime(2009,1,1):
   met_id = 1
  else:
   met_id = 2
  prec = np.ma.getdata(ga.exp("prec.3"))
  tmax = np.ma.getdata(ga.exp("tmax.%d-273.15" % met_id))
  tmin = np.ma.getdata(ga.exp("tmin.%d-273.15" % met_id))
  wind = np.ma.getdata(ga.exp("wind.%d" % met_id))

  #Append to the outgoing file
  prec.tofile(fp)
  tmax.tofile(fp)
  tmin.tofile(fp)
  wind.tofile(fp)

  #Update the time step
  date = date + dt

 #Close the outgoing file
 fp.close()
 ga("close 3")
 ga("close 2")
 ga("close 1")

 return forcing_file

def Prepare_VIC_Forcings_Historical(idate,fdate,dims):

 #Open files
 pgf_ctl = '../DATA/PGF/DAILY/pgf_daily_0.25deg.ctl' 
 ga("xdfopen %s" % pgf_ctl)

 #Forcing_Filename
 forcing_file = '../DATA/VIC/FORCING/DAILY/forcing_daily_%04d%02d%02d' % (idate.year,idate.month,idate.day)
 fp = open(forcing_file,'wb')

 #Extract and print data
 date = idate
 dt = datetime.timedelta(days=1)
 while date <= fdate:
  
  print date
  time = datetime2gradstime(date)
  ga("set time %s" % time)
  prec = np.ma.getdata(ga.exp("prec"))
  tmax = np.ma.getdata(ga.exp("tmax-273.15"))
  tmin = np.ma.getdata(ga.exp("tmin-273.15"))
  wind = np.ma.getdata(ga.exp("wind"))

  #Append to the outgoing file
  prec.tofile(fp)
  tmax.tofile(fp)
  tmin.tofile(fp)
  wind.tofile(fp)
  
  #Update the time step
  date = date + dt

 #Close the outgoing file
 fp.close()
 ga("close 1")

 return forcing_file

def Run_VIC(idate,fdate,dims,dataset,Reprocess_Flag):

 dt = relativedelta.relativedelta(years=7)
 VIC_exe = '../SOURCE/VIC_4.0.5_image_mode/VIC_dev.exe'
 #VIC_exe = '../SOURCE/VIC_4.0.5_image_mode_openmp/VIC_dev.exe'
 VIC_global = '../WORKSPACE/Global_Parameter.txt'

 #If it is pgf
 if dataset == 'pgf':
  idate_all = datetime.datetime(1948,1,1)
  if fdate > datetime.datetime(2008,12,31):
   fdate = datetime.datetime(2008,12,31)
  if idate < datetime.datetime(1948,1,1):
   idate = datetime.datetime(1948,1,1)
  if (fdate - idate).days <= 0:
   return
  file = '../DATA/VIC_PGF/DAILY/output_grid_%04d%02d%02d00' % (fdate.year,fdate.month,fdate.day)
  if os.path.exists(file) and Reprocess_Flag == False:
   return
 #If it is 3b42rt
 if dataset == '3b42rt':
  idate_all = datetime.datetime(2003,1,1)
  if idate < datetime.datetime(2003,1,1):
   idate = datetime.datetime(2003,1,1)
  if (fdate-idate).days <= 0:
   return
  file = '../DATA/VIC_3B42RT/DAILY/output_grid_%04d%02d%02d00' % (fdate.year,fdate.month,fdate.day)
  if os.path.exists(file) and Reprocess_Flag == False:
   return
 #If it is gfsanl
 if dataset == 'gfsanl':
  if idate < datetime.datetime(2010,1,1):
   idate = datetime.datetime(2010,1,1)
  if (fdate-idate).days <= 0:
   return
 if dataset == 'gfs_forecast':
  idate = fdate + datetime.timedelta(days=1)
  fdate = fdate + datetime.timedelta(days=7)
  if idate < datetime.datetime(2013,8,1):
   return
  file = '../DATA/GFS_DERIVED/%04d%02d%02d/output_grid_%04d%02d%02d00' % (idate.year,idate.month,idate.day,fdate.year,fdate.month,fdate.day)
  if os.path.exists(file) and Reprocess_Flag == False:
   return

 #Run model until completed 
 idate_tmp = idate
 while idate_tmp <= fdate:

  fdate_tmp = idate_tmp + dt - datetime.timedelta(days=1)
  if fdate_tmp > fdate:
   fdate_tmp = fdate
  print idate_tmp,fdate_tmp

  #Prepare the VIC global parameter file
  print "Preparing the global parameter file"
  Prepare_VIC_Global_Parameter_File(idate_tmp,fdate_tmp,dims,dataset)

  #Prepare the VIC forcings
  print "Preparing the VIC forcings"
  if dataset == 'pgf':
   forcing_file = Prepare_VIC_Forcings_Historical(idate_tmp,fdate_tmp,dims) #Historical
  if dataset == '3b42rt':
   forcing_file = Prepare_VIC_Forcings_3B42RT(idate_tmp,fdate_tmp,dims) #Historical
  if dataset == 'gfsanl': 
   forcing_file = Prepare_VIC_Forcings_GFSANL(idate_tmp,fdate_tmp,dims) #GFS analysis
  if dataset == 'gfs_forecast':
   forcing_file = Prepare_VIC_Forcings_GFS_forecast(idate_tmp,fdate_tmp,dims) #GFS forecast

  #Run the model
  print "Running VIC"
  os.system('%s -g %s >& ../WORKSPACE/VIC_log.txt' % (VIC_exe,VIC_global))

  #Remove the forcing file
  os.system('rm %s' % forcing_file)

  #Update initial time step
  idate_tmp = fdate_tmp + datetime.timedelta(days=1)

 #Update the control file
 if dataset not in ['pgf','3b42rt']:
  return
 if dataset == 'pgf':
  ctl = '../DATA/VIC_PGF/DAILY/vic_daily_0.25deg.ctl'
 if dataset == '3b42rt':
  ctl = '../DATA/VIC_3B42RT/DAILY/vic_daily_0.25deg.ctl'
 nt = (fdate - idate_all).days + 1
 fp = open(ctl,'w')
 fp.write('dset ^output_grid_%y4%m2%d2%h2\n')
 fp.write('options template\n')
 fp.write('title Gridded VIC Output\n')
 fp.write('undef -999.900000\n')
 fp.write('xdef 296 linear -18.875000 0.250000\n')
 fp.write('ydef 292 linear -34.875000 0.250000\n')
 fp.write('zdef 1 linear 1 1\n')
 fp.write('tdef %d linear %s 24hr\n' % (nt,idate_all.strftime('%d%b%Y')))
 fp.write('vars 18\n')
 fp.write('prec            0  61 prec (mm)\n')
 fp.write('evap            0  57 evap (mm)\n')
 fp.write('runoff          0 235 runoff (mm)\n')
 fp.write('baseflow        0 234 baseflow (mm)\n')
 fp.write('Wdew            0  99 Wdew (mm)\n')
 fp.write('sm1             0 151 soil moisture in layer 1 (mm)\n')
 fp.write('sm2             0 151 soil moisture in layer 2 (mm)\n')
 fp.write('sm3             0 151 soil moisture in layer 3 (mm)\n')
 fp.write('evap_canop      0 200 evap_canop (mm)\n')
 fp.write('evap_veg        0 210 evap_veg (mm)\n')
 fp.write('evap_bare       0 199 evap_bare (mm)\n')
 fp.write('net_short       0 111 net_short (W/m^2)\n')
 fp.write('net_long       0 111 net_long (W/m^2)\n')
 fp.write('r_net           0  99 r_net (W/m^2)\n')
 fp.write('surf_temp           0  99 surf_temp (C)\n')
 fp.write('swq             0  65 swq (mm)\n')
 fp.write('snow_depth      0  66 snow_depth (cm)\n')
 fp.write('snow_canop      0  99 snow_canop (mm)\n')
 fp.write('endvars\n')
 fp.close()

 return

def Extract_VIC_Baseflow_and_Runoff(date,dims):

 ctl_file = '../DATA/VIC/OUTPUT/DAILY/vic_daily_0.25deg.ctl'
 idate = datetime.datetime(1950,1,1)
 if date < idate:
  return
 file = '../DATA/VIC/OUTPUT/DAILY/RUNOFF/runoff_%04d%02d%02d_daily_0.250deg.nc' % (date.year,date.month,date.day)
 if os.path.exists(file):
  return

 #Open file 
 ga("open %s" % ctl_file)
 
 #Set time
 ga("set time %s" % datetime2gradstime(date))

 #Extract data
 data = ga.exp("baseflow") + ga.exp("runoff")

 #Write data
 vars = ['runoff']
 vars_info = ['runoff (mm/day)']
 nt = 1
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file,vars,vars_info,idate,tstep,nt)

 #Write to file
 fp.variables['runoff'][0] = data

 #Close access to all files in grads
 ga("close 1")

 #Finalize and close NETCDF file
 fp.close()

 return

def Run_VDSC(idate,fdate,dims,dataset):

 #If beginning of record start with no initial state
 ist = 1
 if idate == datetime.datetime(1948,1,1):
  ist = -1
 if dataset == 'pgf':
  idate_ctl = datetime.datetime(1948,1,1)
  fdate_ctl = datetime.datetime(2008,12,31)
  if (idate < idate_ctl or fdate > fdate_ctl):
   return
  input_dir = '../DATA/VIC_PGF/DAILY/output_grid_'
  output_dir = '../DATA/ROUTING_VIC_PGF/DAILY'
  state_dir = '../DATA/ROUTING_VIC_PGF/STATE'
  ctl = '../DATA/ROUTING_VIC_PGF/DAILY/Streamflow.ctl'
 if dataset == '3b42rt':
  idate_ctl = datetime.datetime(2003,1,1)
  input_dir = '../DATA/VIC_3B42RT/DAILY/output_grid_'
  output_dir = '../DATA/ROUTING_VIC_3B42RT/DAILY'
  state_dir = '../DATA/ROUTING_VIC_3B42RT/STATE'
  ctl = '../DATA/ROUTING_VIC_3B42RT/DAILY/Streamflow.ctl'
 if dataset == 'gfsanl':
  input_dir = '../DATA/VIC_GFSANL/DAILY/output_grid_'
  output_dir = '../DATA/ROUTING_VIC_GFSANL/DAILY'
  state_dir = '../DATA/ROUTING_VIC_GFSANL/STATE'
  ctl = '../DATA/ROUTING_VIC_GFSANL/DAILY/Streamflow.ctl'
 if dataset == 'gfs_forecast':
  idate = fdate + datetime.timedelta(days=1)
  fdate = idate + datetime.timedelta(days=6)
  idate_ctl = idate
  input_dir = '../DATA/GFS_DERIVED/%04d%02d%02d/output_grid_' % (idate.year,idate.month,idate.day)
  output_dir = '../DATA/GFS_DERIVED/%04d%02d%02d' % (idate.year,idate.month,idate.day)
  state_dir = '../DATA/ROUTING_VIC_3B42RT/STATE'
  ctl = '../DATA/GFS_DERIVED/%04d%02d%02d/Streamflow.ctl' % (idate.year,idate.month,idate.day)

 #Add current directory info
 input_dir = os.getcwd() + '/' + input_dir
 output_dir = os.getcwd() + '/' + output_dir
 state_dir = os.getcwd() + '/' + state_dir

 #Define time info
 iyear = idate.year
 imonth = idate.month
 iday = idate.day
 nt = (fdate - idate).days + 1

 #Run the model
 os.system('perl ../SOURCE/Grid_Routing/scripts/Run_Grid_Routing.pl %d %d %d %d %d %s %s %s' % (iyear,imonth,iday,nt,ist,input_dir,output_dir,state_dir))

 #Update the control file
 nt = (fdate - idate_ctl).days + 1
 fp = open(ctl,'w')
 fp.write('dset ^Streamflow_%y4%m2%d2.bin\n')
 fp.write('options template little_endian\n')
 fp.write('title x,y direction for AFR at 900s resolution \n')
 fp.write('undef -9.99e+08\n')
 fp.write('xdef %d  linear %f %f\n' % (dims['nlon'],dims['minlon'],dims['res']))
 fp.write('ydef %d  linear %f %f\n' % (dims['nlat'],dims['minlat'],dims['res']))
 fp.write('tdef %d linear 00Z%s 1dy\n' % (nt,idate_ctl.strftime('%d%b%Y')))
 fp.write('zdef 1 linear 1 1\n')
 fp.write('vars 1\n')
 fp.write('flw 0 99 Routed Streamflow [m3/s]\n')
 fp.write('endvars\n')
 fp.close()

 return

def Finalize_GFS_forecast(idate,dims):

 print_info_to_command_line("Finalizing the GFS forecast products")

 DATA = {}
 for t in xrange(0,7):
  
  #Define the timestamp
  date = idate + t*datetime.timedelta(days=1)
  gradsdate = datetime2gradstime(date)
  gradsidate = datetime2gradstime(idate)

  #Define all the control files
  ctl_vic = "../DATA/GFS_DERIVED/%04d%02d%02d/output_grid_%04d%02d%02d00.ctl" % (idate.year,idate.month,idate.day,idate.year,idate.month,idate.day)
  ctl_vic_derived = '../DATA/GFS_DERIVED/%04d%02d%02d/vic_derived_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)
  ctl_routing = '../DATA/GFS_DERIVED/%04d%02d%02d/Streamflow.ctl' % (idate.year,idate.month,idate.day)
  ctl_forcing = '../DATA/GFS_BC/gfs_daily_0.250deg_day%d.ctl' % (t+1)
  ctl_spi = '../DATA/GFS_DERIVED/%04d%02d%02d/spi_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)
  ctl_routing_derived = '../DATA/GFS_DERIVED/%04d%02d%02d/routing_vic_derived_daily_0.25deg.ctl' % (idate.year,idate.month,idate.day)
 
  #Create dictionary of ctl info
  ctl_info = {ctl_vic:'open',ctl_vic_derived:'xdfopen',ctl_routing:'open',ctl_forcing:'xdfopen',ctl_spi:'xdfopen',ctl_routing_derived:'xdfopen'}

  #Iterate through all data files extracting the desired data
  for ctl in ctl_info:

   #Open all the control files
   ga("%s %s" % (ctl_info[ctl],ctl))
   
   #Set timestamp
   ga("set time %s" % gradsdate)
   if ctl == ctl_forcing:
    ga("set time %s" % gradsidate)

   #Extract all the variable names
   vars = ga.query("file").vars
   vars_info = ga.query("file").var_info

   #Extract all the info for the variable at that time step
   count = 0
   for var in vars:
    if var == 'prec' and ctl == ctl_vic:
     continue
    if t == 0:
     DATA[var] = {}
     DATA[var]['data'] = []
     DATA[var]['info'] = vars_info[count][3]
    count = count + 1
    DATA[var]['data'].append(ga.exp(var))

   #Close all the control files
   ga("close 1")

 #Write all the data to a file
 file_out = "../DATA/GFS_DERIVED/%04d%02d%02d/gfs_derived_%04d%02d%02d00.nc" % (idate.year,idate.month,idate.day,idate.year,idate.month,idate.day)
 vars = []
 vars_info = []
 for var in DATA:
  vars.append(var)
  vars_info.append(DATA[var]['info'])
 nt = 7
 tstep = 'days'
 fp = Create_NETCDF_File(dims,file_out,vars,vars_info,idate,tstep,nt)

 #Write to file
 for t in xrange(0,7):
  for var in vars:
   fp.variables[var][t] = DATA[var]['data'][t]

 #Close access to file
 fp.close()

 #Update the control file
 idate_all = datetime.datetime(2013,1,1)
 fdate_all = idate
 ndays = (fdate_all - idate_all).days + 1
 ctl = '../DATA/GFS_DERIVED/GFS_DERIVED.ctl'
 fp = open(ctl,'w')
 fp.write('dset ^%e/gfs_derived_%e00.nc\n')
 fp.write('options template\n')
 fp.write('dtype netcdf\n')
 fp.write('title GFS Forecast\n')
 fp.write('undef -9.99e+08\n')
 fp.write('xdef %d  linear %f %f\n' % (dims['nlon'],dims['minlon'],dims['res']))
 fp.write('ydef %d  linear %f %f\n' % (dims['nlat'],dims['minlat'],dims['res']))
 fp.write('tdef %d linear 00Z%s 1dy\n' % (ndays+7-1,idate_all.strftime('%d%b%Y')))
 fp.write('zdef 1 linear 1 1\n')
 fp.write('edef %d\n' % ndays)
 date_tmp = idate_all
 while date_tmp <= fdate_all:
  fp.write('%04d%02d%02d 7 00Z%s 1dy\n' % (date_tmp.year,date_tmp.month,date_tmp.day,date_tmp.strftime('%d%b%Y')))
  date_tmp = date_tmp + relativedelta.relativedelta(days=1)
 fp.write('endedef\n')
 fp.write('vars %d\n' % len(vars))
 for ivar in xrange(0,len(vars)):
  fp.write('%s 1 t,y,x %s\n' % (vars[ivar],vars_info[ivar]))
 fp.write('endvars\n')
 fp.close()

 return
